"""Union the supplied STEP supports with the six lattice graphs; export print-only STL."""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
from pathlib import Path

import cadquery as cq
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

from export_voxel_surface import implicit_surface, validate_surface, write_stl
from mesh_obj3 import GRAPH_OFFSETS, VOXELS, graph_primitives, voxelize_graph, surface_triangles
from run import ROOT, NEW_CASES, SOURCES


def volume(surface):
    mass = vtk.vtkMassProperties()
    mass.SetInputData(surface)
    mass.Update()
    return mass.GetVolume()


def load_supports(step: Path):
    supports = {}
    solids = cq.importers.importStep(str(step)).solids().vals()
    if len(solids) != 2:
        raise ValueError("Expected exactly two STEP support solids")
    for solid in solids:
        if not solid.isValid() or solid.Volume() <= 0:
            raise ValueError("Invalid STEP support solid")
        front = max(
            (f for f in solid.Faces() if f.geomType() == "PLANE"
             and abs(f.Center().x) < 1e-5 and f.normalAt().x > 0.99),
            key=lambda f: f.Area(),
        )
        # The supplied supports have their 25/50 mm interface centered on (base/2, base/2).
        base = int(round(2 * front.Center().y))
        bounds = solid.BoundingBox()
        if (base not in (25, 50) or base in supports
            or abs(front.Center().z - base / 2) > 0.002
            or abs(bounds.xmin + 10) > 1e-5 or abs(bounds.xmax) > 1e-5):
            raise ValueError("Unexpected support placement; expected x=-10..0 and a 25/50 mm base")
        vertices, faces = solid.tessellate(0.02, 0.1)
        points = vtk.vtkPoints()
        points.SetData(numpy_to_vtk(np.array([v.toTuple() for v in vertices]), deep=True))
        cells = vtk.vtkCellArray()
        cells.SetData(
            numpy_to_vtk(np.arange(len(faces) + 1, dtype=np.int64) * 3, deep=True),
            numpy_to_vtk(np.array(faces, dtype=np.int64).ravel(), deep=True),
        )
        mesh = vtk.vtkPolyData()
        mesh.SetPoints(points)
        mesh.SetPolys(cells)
        clean = vtk.vtkCleanPolyData()
        clean.SetInputData(mesh)
        clean.ToleranceIsAbsoluteOn()
        clean.SetAbsoluteTolerance(1e-8)
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputConnection(clean.GetOutputPort())
        normals.SplittingOff()
        normals.ConsistencyOn()
        normals.AutoOrientNormalsOn()
        normals.Update()
        surface = vtk.vtkPolyData()
        surface.DeepCopy(normals.GetOutput())
        validate_surface(surface)
        mesh_volume = volume(surface)
        if abs(mesh_volume / solid.Volume() - 1) > 0.002:
            raise RuntimeError("Support tessellation changed STEP volume by more than 0.2%")
        supports[base] = (surface, solid.Volume())
        print(f"STEP support {base} mm: CAD {solid.Volume():.3f}, mesh {mesh_volume:.3f} mm^3; closed/manifold, no intersections", flush=True)
    return supports


def sample_support(surface, voxel, offsets):
    bounds = np.array(surface.GetBounds()).reshape(3, 2)
    lower = np.floor(bounds[:, 0] / voxel) * voxel - 2 * voxel
    upper = bounds[:, 1] + 2 * voxel
    shape = np.ceil((upper - lower) / voxel).astype(int) + 1
    origin = lower + offsets * voxel
    maximum = origin + (shape - 1) * voxel
    distance = vtk.vtkImplicitPolyDataDistance()
    distance.SetInput(surface)
    if distance.EvaluateFunction((-5, 0, 0)) >= 0:
        raise RuntimeError("Support distance sign is not negative inside the root material")
    sample = vtk.vtkSampleFunction()
    sample.SetImplicitFunction(distance)
    sample.SetModelBounds(*np.column_stack((origin, maximum)).ravel())
    sample.SetSampleDimensions(*shape)
    sample.SetOutputScalarTypeToFloat()
    sample.ComputeNormalsOff()
    sample.CappingOff()
    print(f"Sampling STEP support field: grid {voxel} mm, {tuple(shape)}", flush=True)
    sample.Update()
    field = vtk_to_numpy(sample.GetOutput().GetPointData().GetScalars()).reshape(tuple(shape), order="F").copy()
    np.clip(field, -2 * voxel, 2 * voxel, out=field)
    if not np.isfinite(field).all():
        raise RuntimeError("Nonfinite support distance field")
    skin = implicit_surface(field, lower, voxel, offsets)
    sampled_volume = volume(skin)
    print(f"Support field surface volume {sampled_volume:.3f} mm^3", flush=True)
    return field, lower, sampled_volume


def union_fields(lattice, lattice_lower, support, support_lower, voxel):
    lower = np.minimum(lattice_lower, support_lower)
    starts = [np.rint((item - lower) / voxel).astype(int) for item in (lattice_lower, support_lower)]
    for start, item in zip(starts, (lattice_lower, support_lower)):
        if not np.allclose(lower + start * voxel, item, atol=1e-8, rtol=0):
            raise ValueError("Field grids are not aligned")
    ends = [start + np.array(field.shape) for start, field in zip(starts, (lattice, support))]
    shared_start, shared_end = np.maximum(*starts), np.minimum(*ends)
    if np.any(shared_end <= shared_start):
        raise RuntimeError("Support and lattice have no overlapping grid region")
    shared = [tuple(slice(a, b) for a, b in zip(shared_start - start, shared_end - start)) for start in starts]
    overlap = (lattice[shared[0]] < 0) & (support[shared[1]] < 0)
    count = int(np.count_nonzero(overlap))
    if count < 8:
        raise RuntimeError("Insufficient solid overlap; touching surfaces alone are not a safe joint")
    joined = np.full(tuple(np.maximum(*ends)), 2 * voxel, dtype=np.float32)
    for start, end, field in zip(starts, ends, (lattice, support)):
        region = tuple(slice(a, b) for a, b in zip(start, end))
        np.minimum(joined[region], field, out=joined[region])
    return joined, lower, count * voxel**3


def audit_stl(path: Path, expected_volume: float):
    signed_volume = 0.0
    for vertices in surface_triangles(path):
        a, b, c = (vertices[:, i].astype(np.float64) for i in range(3))
        normal = np.cross(b - a, c - a)
        if np.any(np.einsum("ij,ij->i", normal, normal) == 0):
            raise RuntimeError(f"Written STL has zero-area triangles: {path}")
        signed_volume += np.einsum("ij,ij->", a, np.cross(b, c)) / 6
    if signed_volume <= 0 or abs(signed_volume / expected_volume - 1) > 1e-6:
        raise RuntimeError(f"Written STL volume/winding differs from the validated mesh: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for block in iter(lambda: file.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return {"stl_signed_volume_mm3": signed_volume, "stl_sha256": digest.hexdigest()}


def export(case, support, sampled, step_hash):
    voxel, offsets = VOXELS[case], GRAPH_OFFSETS[case]
    _, cad_volume = support
    support_field, support_lower, support_volume = sampled
    if abs(support_volume / cad_volume - 1) > 0.01:
        raise RuntimeError("Support implicit surface differs from STEP volume by more than 1%")
    lattice, lattice_lower, _ = voxelize_graph(SOURCES[case], voxel, None, offsets=offsets, implicit=True)
    joined, lower, overlap_volume = union_fields(lattice, lattice_lower, support_field, support_lower, voxel)
    del lattice
    surface = implicit_surface(joined, lower, voxel, offsets)
    joined_volume = volume(surface)
    if not np.isfinite(joined_volume) or joined_volume <= max(support_volume, 0):
        raise RuntimeError("Invalid joined volume")
    output = ROOT / "geometry" / f"{case}_with_support.stl"
    write_stl(surface, output)
    audit = audit_stl(output, joined_volume)
    print(f"Wrote {output.name}: {surface.GetNumberOfPolys():,} triangles, {joined_volume:.3f} mm^3, sampled overlap {overlap_volume:.3f} mm^3; all surface checks passed", flush=True)
    return {
        "file": output.name, "source_3mf": SOURCES[case].name, **audit,
        "source_3mf_sha256": hashlib.sha256(SOURCES[case].read_bytes()).hexdigest(),
        "support_step_sha256": step_hash, "support_base_mm": 50 if case.endswith("full") else 25,
        "support_transform": "identity; STEP x=-10..0, lattice root x=0; no radius changes",
        "grid_mm": voxel, "offsets": offsets.tolist(),
        "support_cad_volume_mm3": cad_volume, "support_surface_volume_mm3": support_volume,
        "sampled_overlap_volume_mm3": overlap_volume, "combined_surface_volume_mm3": joined_volume,
        "bounds_mm": list(surface.GetBounds()), "triangles": surface.GetNumberOfPolys(),
        "bytes": output.stat().st_size, "unit": "mm (STL is unitless)",
        "connected_surfaces": 1, "boundary_edges": 0, "nonmanifold_edges": 0, "self_intersecting_faces": 0,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", type=Path, default=ROOT / "geometry/Support.step")
    parser.add_argument("--cases", nargs="+", choices=NEW_CASES, default=NEW_CASES)
    args = parser.parse_args()
    step_hash = hashlib.sha256(args.step.read_bytes()).hexdigest()
    supports = load_supports(args.step)
    sampled_fields = {}
    report = ROOT / "geometry/supported_print_models.json"
    results = json.loads(report.read_text(encoding="utf-8")) if report.is_file() else {}
    for case in args.cases:
        points, _, _ = graph_primitives(SOURCES[case])
        base = 50 if case.endswith("full") else 25
        if not np.allclose(points.min(0), 0) or not np.allclose(points.max(0)[1:], base):
            raise ValueError(f"Unexpected lattice alignment: {case}")
        key = (base, VOXELS[case], tuple(GRAPH_OFFSETS[case]))
        if key not in sampled_fields:
            sampled_fields[key] = sample_support(supports[base][0], VOXELS[case], GRAPH_OFFSETS[case])
        results[case] = export(case, supports[base], sampled_fields[key], step_hash)
        staged = report.with_suffix(".json.tmp")
        staged.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
        staged.replace(report)
        gc.collect()
    if hashlib.sha256(args.step.read_bytes()).hexdigest() != step_hash:
        raise RuntimeError("Support STEP changed during export")


if __name__ == "__main__":
    main()
