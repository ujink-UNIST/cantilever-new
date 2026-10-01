from __future__ import annotations

import argparse
from pathlib import Path

import fast_simplification
import numpy as np
import pymeshfix
import pyvista as pv
import vtk

ROOT = Path(__file__).resolve().parent
CASES = ("modal", "uniform")


def edge_count(mesh: vtk.vtkPolyData, *, boundary: bool = False, nonmanifold: bool = False) -> int:
    edges = vtk.vtkFeatureEdges()
    edges.SetInputData(mesh)
    edges.FeatureEdgesOff()
    edges.ManifoldEdgesOff()
    edges.SetBoundaryEdges(boundary)
    edges.SetNonManifoldEdges(nonmanifold)
    edges.Update()
    return edges.GetOutput().GetNumberOfCells()


def reduce(case: str, reduction: float) -> None:
    source = ROOT / f"{case}_cantilever3.obj"
    output = ROOT / f"{case}_cantilever3_print.obj"

    reader = vtk.vtkOBJReader()
    reader.SetFileName(str(source))
    reader.Update()
    original = pv.wrap(reader.GetOutput())
    if not original.n_cells:
        raise RuntimeError(f"No faces in {source.name}")

    simplified = fast_simplification.simplify_mesh(
        original, target_reduction=reduction, preserve_border=True
    )
    faces = simplified.faces.reshape(-1, 4)[:, 1:].astype(np.int32, copy=False)
    points, faces = pymeshfix.clean_from_arrays(
        np.asarray(simplified.points, dtype=np.float64),
        faces,
        joincomp=False,
        remove_smallest_components=False,
    )
    repaired = pv.PolyData(
        points, np.column_stack((np.full(len(faces), 3, dtype=np.int32), faces)).ravel()
    )

    difference = abs(repaired.volume - original.volume) / original.volume
    if difference > 0.003:
        raise RuntimeError(f"Volume changed by {difference:.3%}; reduce target reduction")
    if edge_count(repaired, boundary=True) or edge_count(repaired, nonmanifold=True):
        raise RuntimeError("Reduced mesh is not closed and manifold")

    writer = vtk.vtkOBJWriter()
    writer.SetFileName(str(output))
    writer.SetInputData(repaired)
    if not writer.Write():
        raise RuntimeError(f"Failed to write {output.name}")

    print(
        f"{case}: {original.n_cells:,} -> {repaired.n_cells:,} faces; "
        f"volume {original.volume:.3f} -> {repaired.volume:.3f} mm^3 "
        f"({difference:.4%}); wrote {output.name}",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Create closed, repaired OBJ files for 3D printing.")
    parser.add_argument("case", nargs="?", choices=CASES, help="omit to process both")
    parser.add_argument("--reduction", type=float, default=0.90, help="requested face reduction (default: 0.90)")
    args = parser.parse_args()
    if not 0 < args.reduction < 1:
        parser.error("reduction must be in (0,1)")
    for case in (args.case,) if args.case else CASES:
        reduce(case, args.reduction)


if __name__ == "__main__":
    main()
