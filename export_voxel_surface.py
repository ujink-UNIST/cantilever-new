from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk

from mesh_obj3 import (
    GRAPH_OFFSETS,
    RADIUS_SCALES,
    REFERENCE_VOLUMES,
    VOXELS,
    voxelize_graph,
)
from run import ROOT, SOURCES

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


def export(case: str, output: Path) -> None:
    voxel = VOXELS[case]
    offsets = GRAPH_OFFSETS[case]
    solid, lower, target_volume = voxelize_graph(
        SOURCES[case], voxel, REFERENCE_VOLUMES[case], RADIUS_SCALES.get(case, 1.0), offsets
    )
    padded = np.pad(solid, 1)
    image = vtk.vtkImageData()
    image.SetDimensions(*padded.shape)
    image.SetSpacing(voxel, voxel, voxel)
    image.SetOrigin(*(lower + (offsets - 1) * voxel))
    image.GetPointData().SetScalars(numpy_to_vtk(padded.ravel(order="F"), deep=True))

    contour = vtk.vtkFlyingEdges3D()
    contour.SetInputData(image)
    contour.SetValue(0, 0.5)
    contour.Update()
    surface = contour.GetOutput()
    if not surface.GetNumberOfPolys():
        raise RuntimeError("No exterior surface was generated")

    connectivity = vtk.vtkPolyDataConnectivityFilter()
    connectivity.SetInputData(surface)
    connectivity.SetExtractionModeToAllRegions()
    connectivity.Update()
    if connectivity.GetNumberOfExtractedRegions() != 1:
        raise RuntimeError(f"Expected one connected surface, found {connectivity.GetNumberOfExtractedRegions()}")
    if edge_count(surface, boundary=True) or edge_count(surface, nonmanifold=True):
        raise RuntimeError("Generated surface is not closed and manifold")

    mass = vtk.vtkMassProperties()
    mass.SetInputData(surface)
    mass.Update()
    output.parent.mkdir(parents=True, exist_ok=True)
    writer = vtk.vtkSTLWriter()
    writer.SetFileName(str(output))
    writer.SetFileTypeToBinary()
    writer.SetInputData(surface)
    if not writer.Write():
        raise RuntimeError(f"Failed to write {output}")
    print(
        f"Wrote {output}: {surface.GetNumberOfPoints():,} points, "
        f"{surface.GetNumberOfPolys():,} triangles, volume {mass.GetVolume():.3f} mm^3 "
        f"({abs(mass.GetVolume() - target_volume) / target_volume:.3%} from voxel target)"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract a closed exterior STL from a direct 3MF voxel union.")
    parser.add_argument("case", choices=CASES, nargs="?", default="modal_cantilever4")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or ROOT / "geometry" / f"{args.case}_direct_voxel_surface.stl"
    export(args.case, output.resolve())


if __name__ == "__main__":
    main()
