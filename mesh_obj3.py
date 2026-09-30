from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pygalmesh

ROOT = Path(__file__).resolve().parent
SOURCES = {
    "modal": ROOT / "modal_cantilever3.obj",
    "uniform": ROOT / "uniform_cantilever3.obj",
}
VOXEL = 0.25
PAD = 2
OFFSETS = np.array([0.419, 0.371, 0.583])  # Avoid rays through source edges.
FACE_BLOCK = 500_000


def mesh_case(case: str) -> None:
    source = SOURCES[case]
    output = ROOT / f"{source.stem}_volume.msh"

    with source.open() as obj:
        header = obj.readline()
        match = re.fullmatch(r"#\s*(\d+)\s+vertices,\s*(\d+)\s+faces\s*\n?", header)
        if not match:
            raise RuntimeError(f"Unexpected OBJ header: {header!r}")
        vertex_count, face_count = map(int, match.groups())
        points = np.loadtxt(
            obj, dtype=np.float32, usecols=(1, 2, 3), max_rows=vertex_count, ndmin=2
        )
        if len(points) != vertex_count:
            raise RuntimeError(f"Expected {vertex_count} vertices, found {len(points)}")

        lower = np.floor(points.min(axis=0) / VOXEL) * VOXEL
        shape = np.ceil((points.max(axis=0) - lower) / VOXEL).astype(int)
        intersected_rays: list[np.ndarray] = []
        intersection_x: list[np.ndarray] = []
        source_volume = 0.0
        loaded_faces = 0

        while loaded_faces < face_count:
            wanted = min(FACE_BLOCK, face_count - loaded_faces)
            faces = np.loadtxt(
                obj, dtype=np.int32, comments="#", usecols=(1, 2, 3), max_rows=wanted, ndmin=2
            ) - 1
            if len(faces) != wanted:
                raise RuntimeError(f"Expected {face_count} faces, found {loaded_faces + len(faces)}")
            loaded_faces += len(faces)

            a, b, c = (points[faces[:, i]] for i in range(3))
            source_volume += np.einsum(
                "ij,ij->", a.astype(np.float64), np.cross(b, c)
            ) / 6

            yz_min = np.minimum(np.minimum(a[:, 1:], b[:, 1:]), c[:, 1:])
            yz_max = np.maximum(np.maximum(a[:, 1:], b[:, 1:]), c[:, 1:])
            if np.any(yz_max - yz_min >= VOXEL):
                raise RuntimeError("An OBJ triangle spans a full voxel; reduce VOXEL or refine the surface")
            ij = np.ceil((yz_min - lower[1:]) / VOXEL - OFFSETS[1:]).astype(np.int32)
            q = lower[1:] + (ij + OFFSETS[1:]) * VOXEL
            den = (b[:, 1] - c[:, 1]) * (a[:, 2] - c[:, 2]) + (
                c[:, 2] - b[:, 2]
            ) * (a[:, 1] - c[:, 1])
            valid = (
                (ij[:, 0] >= 0)
                & (ij[:, 0] < shape[1])
                & (ij[:, 1] >= 0)
                & (ij[:, 1] < shape[2])
                & (q[:, 0] <= yz_max[:, 0])
                & (q[:, 1] <= yz_max[:, 1])
                & (np.abs(den) > 1e-12)
            )
            rows = np.flatnonzero(valid)
            wa = (
                (b[rows, 1] - c[rows, 1]) * (q[rows, 1] - c[rows, 2])
                + (c[rows, 2] - b[rows, 2]) * (q[rows, 0] - c[rows, 1])
            ) / den[rows]
            wb = (
                (c[rows, 1] - a[rows, 1]) * (q[rows, 1] - c[rows, 2])
                + (a[rows, 2] - c[rows, 2]) * (q[rows, 0] - c[rows, 1])
            ) / den[rows]
            inside = (wa >= 0) & (wb >= 0) & (wa + wb <= 1)
            rows, wa, wb = rows[inside], wa[inside], wb[inside]
            intersected_rays.append(ij[rows, 0] * shape[2] + ij[rows, 1])
            intersection_x.append(
                wa * a[rows, 0] + wb * b[rows, 0] + (1 - wa - wb) * c[rows, 0]
            )
            print(f"Read {loaded_faces:,}/{face_count:,} faces", flush=True)

    rays = np.concatenate(intersected_rays)
    xs = np.concatenate(intersection_x)
    order = np.lexsort((xs, rays))
    rays, xs = rays[order], xs[order]
    solid = np.zeros(tuple(shape), dtype=np.uint8)
    boundaries = np.flatnonzero(np.diff(rays)) + 1
    odd_rays = 0
    for begin, end in zip(np.r_[0, boundaries], np.r_[boundaries, len(rays)]):
        ray = rays[begin]
        ray_xs = xs[begin:end]
        ray_xs = ray_xs[np.r_[True, np.diff(ray_xs) > 1e-4]]
        if len(ray_xs) % 2:
            odd_rays += 1
            continue
        y, z = divmod(int(ray), int(shape[2]))
        for left, right in ray_xs.reshape(-1, 2):
            lo = max(0, int(np.ceil((left - lower[0]) / VOXEL - OFFSETS[0])))
            hi = min(shape[0], int(np.ceil((right - lower[0]) / VOXEL - OFFSETS[0])))
            solid[lo:hi, y, z] = 1

    assert odd_rays < max(10, len(boundaries) // 100)
    voxel_volume = int(solid.sum()) * VOXEL**3
    print(
        f"Source/voxel volume: {abs(source_volume):.1f}/{voxel_volume:.1f}; "
        f"skipped odd rays: {odd_rays}"
    )
    assert abs(voxel_volume - abs(source_volume)) / abs(source_volume) < 0.15

    mesh = pygalmesh.generate_from_array(
        np.pad(solid, PAD),
        (VOXEL,) * 3,
        min_facet_angle=25,
        max_radius_surface_delaunay_ball=0.75,
        max_facet_distance=0.1,
        max_circumradius_edge_ratio=3,
        max_cell_circumradius=1.0,
        verbose=True,
        seed=1,
    )
    mesh.points += lower + (OFFSETS - PAD) * VOXEL
    mesh.point_data.clear()
    mesh.cell_data.clear()

    tetra = mesh.get_cells_type("tetra")
    a, b, c, d = (mesh.points[tetra[:, i]] for i in range(4))
    mesh_volume = np.abs(
        np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))
    ).sum() / 6
    assert len(tetra) > 0
    assert abs(mesh_volume - abs(source_volume)) / abs(source_volume) < 0.1
    mesh.write(output, file_format="gmsh22", binary=True)
    print(f"Wrote {output.name}: {len(mesh.points)} points, {len(tetra)} tetra, volume {mesh_volume:.1f}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Voxelize and tetrahedralize a cantilever3 OBJ.")
    parser.add_argument("case", choices=SOURCES)
    args = parser.parse_args()
    mesh_case(args.case)


if __name__ == "__main__":
    main()
