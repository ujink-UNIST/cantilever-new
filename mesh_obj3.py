from __future__ import annotations

import argparse
import re
import struct
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile

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
    "modal_cantilever3": 0.2,
    "uniform_cantilever3": 0.2,
    "modal_cantilever4": 0.125,
    "uniform_cantilever4": 0.1,
    "modal-8": 0.2,
    "uniform-8": 0.2,
}
PAD = 2
OFFSETS = np.array([0.419, 0.371, 0.583])  # Avoid rays through source edges.
GRAPH_OFFSETS = {
    "modal_cantilever3": np.array([0.547, 0.673, 0.293]),
    "uniform_cantilever3": np.array([0.547, 0.673, 0.293]),
    "modal_cantilever4": np.array([0.389, 0.811, 0.137]),
    "uniform_cantilever4": np.array([0.389, 0.811, 0.137]),
}
RADIUS_SCALES = {"uniform_cantilever3": 1.015, "uniform_cantilever4": 1.005}
FACE_BLOCK = 500_000
STL_DTYPE = np.dtype([("normal", "<f4", (3,)), ("vertices", "<f4", (3, 3)), ("attribute", "<u2")])
CORE = "http://schemas.microsoft.com/3dmanufacturing/core/2015/02"
BEAM = "http://schemas.microsoft.com/3dmanufacturing/beamlattice/2017/02"
BALL = "http://schemas.microsoft.com/3dmanufacturing/beamlattice/balls/2020/07"
REFERENCE_VOLUMES = {
    "modal_cantilever3": 86399.552,
    "uniform_cantilever3": 86399.552,
    "modal_cantilever4": 10805.796875,
    "uniform_cantilever4": 10805.796875,
}


def graph_primitives(source: Path):
    with ZipFile(source) as archive:
        models = [name for name in archive.namelist() if name.lower().endswith(".model")]
        if len(models) != 1:
            raise ValueError(f"Expected one 3MF model, found {len(models)}")
        model = ET.fromstring(archive.read(models[0]))
    if model.get("unit", "millimeter") != "millimeter":
        raise ValueError("Only millimeter 3MF models are supported")
    objects = model.findall(f".//{{{CORE}}}object")
    if len(objects) != 1:
        raise ValueError(f"Expected one 3MF object, found {len(objects)}")
    obj = objects[0]
    item = model.find(f".//{{{CORE}}}item")
    if item is None or item.get("transform") is not None:
        raise ValueError("Expected one untransformed 3MF build item")
    vertices = obj.findall(f".//{{{CORE}}}vertex")
    points = np.array(
        [[float(vertex.get(axis)) for axis in "xyz"] for vertex in vertices], dtype=np.float64
    )
    lattice = obj.find(f".//{{{BEAM}}}beamlattice")
    if lattice is None or lattice.get("cap", "sphere") != "butt":
        raise ValueError("Expected a butt-cap 3MF beam lattice")
    default_radius = float(lattice.get("radius"))
    beams = np.array(
        [
            [
                int(beam.get("v1")), int(beam.get("v2")),
                float(beam.get("r1", default_radius)), float(beam.get("r2", default_radius)),
            ]
            for beam in lattice.findall(f".//{{{BEAM}}}beam")
        ],
        dtype=np.float64,
    )
    default_ball_radius = float(lattice.get(f"{{{BALL}}}ballradius", default_radius))
    balls = np.array(
        [
            [int(ball.get("vindex")), float(ball.get("r", default_ball_radius))]
            for ball in lattice.findall(f".//{{{BALL}}}ball")
        ],
        dtype=np.float64,
    )
    if (
        points.ndim != 2 or points.shape[1] != 3 or not len(beams) or not len(balls)
        or not np.isfinite(points).all() or not np.isfinite(beams).all() or not np.isfinite(balls).all()
        or beams[:, :2].min() < 0 or beams[:, :2].max() >= len(points)
        or balls[:, 0].min() < 0 or balls[:, 0].max() >= len(points)
        or beams[:, 2:].min() <= 0 or balls[:, 1].min() <= 0
    ):
        raise ValueError("Invalid 3MF beam or ball data")
    return points, beams, balls


def voxelize_graph(
    source: Path, voxel: float, reference_volume: float,
    radius_scale: float = 1.0, offsets: np.ndarray = OFFSETS,
    *, implicit: bool = False,
):
    if not np.isfinite(voxel) or voxel <= 0 or not np.isfinite(radius_scale) or radius_scale <= 0:
        raise ValueError("Voxel size and radius scale must be positive and finite")
    points, beams, balls = graph_primitives(source)
    beams[:, 2:] *= radius_scale
    balls[:, 1] *= radius_scale
    max_radius = max(beams[:, 2:].max(), balls[:, 1].max())
    lower = np.floor((points.min(axis=0) - max_radius) / voxel) * voxel
    upper = points.max(axis=0) + max_radius
    band = 2 * voxel if implicit else 0.0
    lower -= band
    upper += band
    shape = np.ceil((upper - lower) / voxel).astype(int)
    # ponytail: a two-cell distance band covers zero crossings; widen for offset surfaces.
    solid = np.full(tuple(shape), band, dtype=np.float32) if implicit else np.zeros(tuple(shape), dtype=bool)

    def index_bounds(minimum, maximum):
        begin = np.ceil((minimum - lower) / voxel - offsets).astype(int)
        end = np.floor((maximum - lower) / voxel - offsets).astype(int) + 1
        return np.maximum(begin, 0), np.minimum(end, shape)

    def coordinates(begin, end):
        return [
            lower[axis] + (np.arange(begin[axis], end[axis]) + offsets[axis]) * voxel
            for axis in range(3)
        ]

    for number, (v1, v2, r1, r2) in enumerate(beams, 1):
        start, end = points[int(v1)], points[int(v2)]
        axis = end - start
        length_squared = axis @ axis
        if length_squared <= 0:
            raise ValueError("Zero-length 3MF beam")
        radius = max(r1, r2) + band
        begin, finish = index_bounds(np.minimum(start, end) - radius, np.maximum(start, end) + radius)
        x, y, z = coordinates(begin, finish)
        dx = x[:, None, None] - start[0]
        dy = y[None, :, None] - start[1]
        dz = z[None, None, :] - start[2]
        projection = dx * axis[0] + dy * axis[1] + dz * axis[2]
        fraction = projection / length_squared
        local_radius = r1 + (r2 - r1) * fraction
        distance_squared = dx * dx + dy * dy + dz * dz - projection * projection / length_squared
        slices = tuple(slice(begin[axis], finish[axis]) for axis in range(3))
        if implicit:
            length = np.sqrt(length_squared)
            radial = (np.sqrt(np.maximum(distance_squared, 0)) - local_radius) / np.sqrt(
                1 + (r2 - r1) ** 2 / length_squared
            )
            cap = np.maximum(-projection, projection - length_squared) / length
            # Exact capped-cylinder SDF for r1=r2; sign-correct field for tapered beams.
            distance = np.hypot(np.maximum(radial, 0), np.maximum(cap, 0)) + np.minimum(
                np.maximum(radial, cap), 0
            )
            np.minimum(solid[slices], distance, out=solid[slices])
        else:
            inside = (
                (projection >= 0) & (projection <= length_squared)
                & (distance_squared <= local_radius * local_radius)
            )
            solid[slices] |= inside
        if number % 1000 == 0 or number == len(beams):
            print(f"Rasterized {number:,}/{len(beams):,} beams", flush=True)

    for vindex, radius in balls:
        center = points[int(vindex)]
        begin, finish = index_bounds(center - radius - band, center + radius + band)
        x, y, z = coordinates(begin, finish)
        distance_squared = (
            (x[:, None, None] - center[0]) ** 2
            + (y[None, :, None] - center[1]) ** 2
            + (z[None, None, :] - center[2]) ** 2
        )
        slices = tuple(slice(begin[axis], finish[axis]) for axis in range(3))
        if implicit:
            np.minimum(solid[slices], np.sqrt(distance_squared) - radius, out=solid[slices])
        else:
            solid[slices] |= distance_squared <= radius * radius

    voxel_volume = int(np.count_nonzero(solid < 0) if implicit else solid.sum()) * voxel**3
    difference = abs(voxel_volume - reference_volume) / reference_volume
    print(
        f"3MF target/voxel volume: {reference_volume:.3f}/{voxel_volume:.3f} mm^3 "
        f"({difference:.3%}); radius scale {radius_scale:.6f}; "
        f"{len(points):,} vertices, {len(beams):,} beams, {len(balls):,} balls"
    )
    if difference > 0.01:
        raise RuntimeError(f"Sampled volume differs from the reference target by {difference:.3%}")
    return solid if implicit else solid.astype(np.uint8), lower, reference_volume


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


def mesh_case(case: str, voxel_only: bool = False) -> None:
    source = SOURCES[case]
    output = ROOT / f"{case}_volume.msh"
    scale, voxel = SCALES[case], VOXELS[case]
    offsets = GRAPH_OFFSETS.get(case, OFFSETS)
    if source.suffix.lower() == ".3mf":
        solid, lower, source_volume = voxelize_graph(
            source, voxel, REFERENCE_VOLUMES[case], RADIUS_SCALES.get(case, 1.0), offsets
        )
    else:
        solid, lower, source_volume = voxelize(source, voxel)
    if voxel_only:
        return

    import pygalmesh  # Needed only for tetrahedralization, not voxel validation.
    mesh = pygalmesh.generate_from_array(
        np.pad(solid, PAD), (voxel,) * 3,
        min_facet_angle=25,
        max_radius_surface_delaunay_ball=0.75 * scale,
        max_facet_distance=0.1 * scale,
        max_circumradius_edge_ratio=3,
        max_cell_circumradius=1.0 * scale,
        verbose=True, seed=1,
    )
    mesh.points += lower + (offsets - PAD) * voxel
    mesh.point_data.clear()
    mesh.cell_data.clear()
    tetra = mesh.get_cells_type("tetra")
    a, b, c, d = (mesh.points[tetra[:, i]] for i in range(4))
    mesh_volume = np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))).sum() / 6
    if not len(tetra):
        raise RuntimeError("CGAL generated no tetrahedra")
    volume_difference = abs(mesh_volume - source_volume) / source_volume
    if volume_difference > 0.03:
        raise RuntimeError(f"Tetra volume differs from the target by {volume_difference:.3%}")
    mesh.write(output, file_format="gmsh22", binary=True)
    print(
        f"Wrote {output.name}: {len(mesh.points)} points, {len(tetra)} tetra, "
        f"volume {mesh_volume:.1f} ({volume_difference:.3%} from target)"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Voxelize and tetrahedralize a 3MF beam graph or binary STL.")
    parser.add_argument("case", choices=SOURCES)
    parser.add_argument("--voxel-only", action="store_true", help="validate voxel volume without CGAL")
    args = parser.parse_args()
    mesh_case(args.case, args.voxel_only)


if __name__ == "__main__":
    main()
