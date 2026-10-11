from __future__ import annotations

import argparse
import re
import struct
from pathlib import Path

import numpy as np

from run import SOURCES

ROOT = Path(__file__).resolve().parent
SCALES = {
    "modal_cantilever3": 1.0,
    "uniform_cantilever3": 1.0,
    "modal_cantilever4": 0.5,
    "uniform_cantilever4": 0.5,
    "modal-8": 1.0,
    "uniform-8": 1.0,
}
VOXELS = {
    "modal_cantilever3": 0.25,
    "uniform_cantilever3": 0.25,
    "modal_cantilever4": 0.125,
    "uniform_cantilever4": 0.2,
    "modal-8": 0.2,
    "uniform-8": 0.2,
}
PAD = 2
OFFSETS = np.array([0.419, 0.371, 0.583])  # Avoid rays through source edges.
FACE_BLOCK = 500_000
STL_DTYPE = np.dtype([("normal", "<f4", (3,)), ("vertices", "<f4", (3, 3)), ("attribute", "<u2")])


def surface_triangles(source: Path):
    """Stream binary STL or the existing indexed, header-counted triangle OBJ."""
    if source.suffix.lower() == ".stl":
        with source.open("rb") as file:
            header = file.read(84)
        if len(header) != 84:
            raise ValueError(f"Incomplete binary STL: {source}")
        count = struct.unpack("<I", header[80:84])[0]
        if not count or source.stat().st_size != 84 + count * STL_DTYPE.itemsize:
            raise ValueError(f"Expected complete binary STL, not ASCII STL: {source}")
        records = np.memmap(source, dtype=STL_DTYPE, offset=84, mode="r", shape=(count,))
        for start in range(0, count, FACE_BLOCK):
            vertices = records["vertices"][start:start + FACE_BLOCK]
            if not np.isfinite(vertices).all():
                raise ValueError(f"Nonfinite STL coordinates: {source}")
            yield vertices
        return
    if source.suffix.lower() != ".obj":
        raise ValueError(f"Unsupported surface format: {source.suffix}")
    with source.open() as obj:
        header = obj.readline()
        match = re.fullmatch(r"#\s*(\d+)\s+vertices,\s*(\d+)\s+faces\s*\n?", header)
        if not match:
            raise ValueError(f"Unexpected OBJ header: {header!r}")
        vertex_count, face_count = map(int, match.groups())
        if not vertex_count or not face_count:
            raise ValueError("Empty OBJ")
        points = np.loadtxt(obj, dtype=np.float32, usecols=(1, 2, 3), max_rows=vertex_count, ndmin=2)
        if len(points) != vertex_count or not np.isfinite(points).all():
            raise ValueError("Incomplete or nonfinite OBJ vertices")
        for start in range(0, face_count, FACE_BLOCK):
            wanted = min(FACE_BLOCK, face_count - start)
            faces = np.loadtxt(obj, dtype=np.int32, comments="#", usecols=(1, 2, 3), max_rows=wanted, ndmin=2) - 1
            if len(faces) != wanted or faces.min() < 0 or faces.max() >= vertex_count:
                raise ValueError("Incomplete OBJ faces or invalid vertex index")
            yield points[faces]


def voxelize(source: Path, voxel: float):
    if not np.isfinite(voxel) or voxel <= 0:
        raise ValueError("Voxel size must be positive and finite")
    stamp = (source.stat().st_size, source.stat().st_mtime_ns)
    lower = np.full(3, np.inf)
    upper = -lower.copy()
    face_count = 0
    # ponytail: two streaming passes avoid a multi-GB vertex buffer; cache bounds if IO dominates.
    for vertices in surface_triangles(source):
        lower = np.minimum(lower, vertices.min(axis=(0, 1)))
        upper = np.maximum(upper, vertices.max(axis=(0, 1)))
        face_count += len(vertices)
    if not face_count:
        raise ValueError("Empty surface")
    lower = np.floor(lower / voxel) * voxel
    shape = np.ceil((upper - lower) / voxel).astype(int)
    if np.any(shape <= 0):
        raise ValueError("Surface does not enclose a 3D domain")
    intersected_rays: list[np.ndarray] = []
    intersection_x: list[np.ndarray] = []
    source_volume = 0.0
    loaded_faces = 0
    for vertices in surface_triangles(source):
        loaded_faces += len(vertices)
        a, b, c = (vertices[:, i].astype(np.float64) for i in range(3))
        source_volume += np.einsum("ij,ij->", a, np.cross(b, c)) / 6
        yz_min = np.minimum(np.minimum(a[:, 1:], b[:, 1:]), c[:, 1:])
        yz_max = np.maximum(np.maximum(a[:, 1:], b[:, 1:]), c[:, 1:])
        if np.any(yz_max - yz_min >= voxel):
            raise RuntimeError("A surface triangle spans a full voxel; increase VOXEL or refine the surface")
        ij = np.ceil((yz_min - lower[1:]) / voxel - OFFSETS[1:]).astype(np.int32)
        q = lower[1:] + (ij + OFFSETS[1:]) * voxel
        den = (b[:, 1] - c[:, 1]) * (a[:, 2] - c[:, 2]) + (c[:, 2] - b[:, 2]) * (a[:, 1] - c[:, 1])
        valid = (
            (ij[:, 0] >= 0) & (ij[:, 0] < shape[1])
            & (ij[:, 1] >= 0) & (ij[:, 1] < shape[2])
            & (q[:, 0] <= yz_max[:, 0]) & (q[:, 1] <= yz_max[:, 1])
            & (np.abs(den) > 1e-12)
        )
        rows = np.flatnonzero(valid)
        wa = ((b[rows, 1] - c[rows, 1]) * (q[rows, 1] - c[rows, 2])
              + (c[rows, 2] - b[rows, 2]) * (q[rows, 0] - c[rows, 1])) / den[rows]
        wb = ((c[rows, 1] - a[rows, 1]) * (q[rows, 1] - c[rows, 2])
              + (a[rows, 2] - c[rows, 2]) * (q[rows, 0] - c[rows, 1])) / den[rows]
        inside = (wa >= 0) & (wb >= 0) & (wa + wb <= 1)
        rows, wa, wb = rows[inside], wa[inside], wb[inside]
        intersected_rays.append(ij[rows, 0] * shape[2] + ij[rows, 1])
        intersection_x.append(wa * a[rows, 0] + wb * b[rows, 0] + (1 - wa - wb) * c[rows, 0])
        print(f"Read {loaded_faces:,}/{face_count:,} faces", flush=True)
    if loaded_faces != face_count or stamp != (source.stat().st_size, source.stat().st_mtime_ns):
        raise RuntimeError("Surface changed while reading")
    if not np.isfinite(source_volume) or source_volume == 0:
        raise ValueError("Surface has invalid or zero signed volume")
    rays = np.concatenate(intersected_rays)
    xs = np.concatenate(intersection_x)
    if not len(rays):
        raise ValueError("No surface intersections at this voxel size")
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
            lo = max(0, int(np.ceil((left - lower[0]) / voxel - OFFSETS[0])))
            hi = min(shape[0], int(np.ceil((right - lower[0]) / voxel - OFFSETS[0])))
            solid[lo:hi, y, z] = 1
    assert odd_rays < max(10, len(boundaries) // 100)
    voxel_volume = int(solid.sum()) * voxel**3
    print(f"Source/voxel volume: {abs(source_volume):.1f}/{voxel_volume:.1f}; skipped odd rays: {odd_rays}")
    assert abs(voxel_volume - abs(source_volume)) / abs(source_volume) < 0.15
    return solid, lower, abs(source_volume)


def mesh_case(case: str) -> None:
    import pygalmesh  # Needed only for tetrahedralization, not surface-reader checks.

    source = SOURCES[case]
    output = ROOT / f"{source.stem}_volume.msh"
    scale, voxel = SCALES[case], VOXELS[case]
    solid, lower, source_volume = voxelize(source, voxel)
    mesh = pygalmesh.generate_from_array(
        np.pad(solid, PAD), (voxel,) * 3,
        min_facet_angle=25,
        max_radius_surface_delaunay_ball=0.75 * scale,
        max_facet_distance=0.1 * scale,
        max_circumradius_edge_ratio=3,
        max_cell_circumradius=1.0 * scale,
        verbose=True, seed=1,
    )
    mesh.points += lower + (OFFSETS - PAD) * voxel
    mesh.point_data.clear()
    mesh.cell_data.clear()
    tetra = mesh.get_cells_type("tetra")
    a, b, c, d = (mesh.points[tetra[:, i]] for i in range(4))
    mesh_volume = np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))).sum() / 6
    assert len(tetra) > 0
    assert abs(mesh_volume - source_volume) / source_volume < 0.1
    mesh.write(output, file_format="gmsh22", binary=True)
    print(f"Wrote {output.name}: {len(mesh.points)} points, {len(tetra)} tetra, volume {mesh_volume:.1f}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Voxelize and tetrahedralize a binary STL or triangle OBJ.")
    parser.add_argument("case", choices=SOURCES)
    args = parser.parse_args()
    mesh_case(args.case)


if __name__ == "__main__":
    main()
