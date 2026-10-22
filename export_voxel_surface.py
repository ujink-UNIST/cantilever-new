from __future__ import annotations

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pymeshlab
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

from mesh_obj3 import (
    GRAPH_OFFSETS,
    RADIUS_SCALES,
    REFERENCE_VOLUMES,
    VOXELS,
    voxelize_graph,
)
from run import ROOT, SOURCES, NEW_CASES

CASES = tuple(case for case, source in SOURCES.items() if source.suffix.lower() == ".3mf")


def edge_count(surface: vtk.vtkPolyData, *, boundary: bool = False, nonmanifold: bool = False) -> int:
    edges = vtk.vtkFeatureEdges()
    edges.SetInputData(surface)
    edges.FeatureEdgesOff()
    edges.ManifoldEdgesOff()
    edges.SetBoundaryEdges(boundary)
    edges.SetNonManifoldEdges(nonmanifold)
    edges.Update()
    return edges.GetOutput().GetNumberOfCells()


def implicit_surface(
    field: np.ndarray, lower: np.ndarray, voxel: float, offsets: np.ndarray, *, validate: bool = True,
):
    image = vtk.vtkImageData()
    image.SetDimensions(*field.shape)
    image.SetSpacing(voxel, voxel, voxel)
    image.SetOrigin(*(lower + offsets * voxel))
    # VTK expects high values inside to orient triangle winding outward.
    image.GetPointData().SetScalars(numpy_to_vtk((-field).ravel(order="F")))

    contour = vtk.vtkFlyingEdges3D()
    contour.SetInputData(image)
    contour.SetValue(0, 0.0)
    contour.Update()
    surface = contour.GetOutput()
    if not surface.GetNumberOfPolys():
        raise RuntimeError("No exterior surface was generated")
    if not validate:
        return surface  # Volume-only iterations; validate the final surface before publishing.

    connectivity = vtk.vtkPolyDataConnectivityFilter()
    connectivity.SetInputData(surface)
    connectivity.SetExtractionModeToAllRegions()
    connectivity.Update()
    if connectivity.GetNumberOfExtractedRegions() != 1:
        raise RuntimeError(f"Expected one connected surface, found {connectivity.GetNumberOfExtractedRegions()}")
    if edge_count(surface, boundary=True) or edge_count(surface, nonmanifold=True):
        raise RuntimeError("Generated surface is not closed and manifold")

    meshset = pymeshlab.MeshSet()
    meshset.add_mesh(pymeshlab.Mesh(
        vertex_matrix=vtk_to_numpy(surface.GetPoints().GetData()).astype(np.float64),
        face_matrix=vtk_to_numpy(surface.GetPolys().GetConnectivityArray()).reshape(-1, 3).astype(np.int32),
    ))
    meshset.compute_selection_by_self_intersections_per_face()
    intersections = meshset.current_mesh().selected_face_number()
    if intersections:
        raise RuntimeError(f"Generated surface has {intersections} self-intersecting faces")
    return surface


def export(case: str, output: Path) -> None:
    voxel = VOXELS[case]
    offsets = GRAPH_OFFSETS[case]
    field, lower, target_volume = voxelize_graph(
        SOURCES[case], voxel, REFERENCE_VOLUMES[case], RADIUS_SCALES.get(case, 1.0), offsets,
        implicit=True,
    )
    surface = implicit_surface(field, lower, voxel, offsets)
    mass = vtk.vtkMassProperties()
    mass.SetInputData(surface)
    mass.Update()
    volume_difference = abs(mass.GetVolume() - target_volume) / target_volume
    if volume_difference > 0.03:
        raise RuntimeError(f"Surface volume differs from the reference target by {volume_difference:.3%}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent) as folder:
        staged = Path(folder) / output.name
        writer = vtk.vtkSTLWriter()
        writer.SetFileName(str(staged))
        writer.SetFileTypeToBinary()
        writer.SetInputData(surface)
        if not writer.Write() or not staged.is_file() or staged.stat().st_size != 84 + 50 * surface.GetNumberOfPolys():
            raise RuntimeError(f"Failed to write a complete binary STL: {output}")
        staged.replace(output)
    print(
        f"Wrote {output}: {surface.GetNumberOfPoints():,} points, "
        f"{surface.GetNumberOfPolys():,} triangles, volume {mass.GetVolume():.3f} mm^3 "
        f"({volume_difference:.3%} from reference target); boundary/nonmanifold edges 0/0, "
        "self-intersecting faces 0, one component"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract a closed STL from a continuous 3MF implicit field.")
    parser.add_argument("case", choices=CASES, nargs="?", default="modal_cantilever4")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--all-new", action="store_true", help="Export all six new matched modal/uniform graphs")
    args = parser.parse_args()
    if args.all_new and args.output:
        parser.error("--output is only available for a single model")
    for case in NEW_CASES if args.all_new else (args.case,):
        output = args.output or ROOT / "geometry" / f"{case}_implicit_surface.stl"
        export(case, output.resolve())


if __name__ == "__main__":
    main()
