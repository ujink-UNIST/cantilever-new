"""Run: python test_supported_stl.py (small analytic field checks)."""
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import vtk

from export_supported_stl import union_fields, audit_stl
from export_voxel_surface import write_stl


def main():
    voxel = 0.125
    lattice = np.full((7, 5, 5), 2 * voxel, dtype=np.float32)
    support = np.full((5, 5, 5), 2 * voxel, dtype=np.float32)
    lattice[:4, 1:4, 1:4] = -0.1
    support[2:, 1:4, 1:4] = -0.1
    lattice_lower = np.zeros(3)
    support_lower = np.array([-2 * voxel, 0, 0])
    joined, lower, overlap_volume = union_fields(lattice, lattice_lower, support, support_lower, voxel)
    assert np.allclose(lower, support_lower)
    assert joined.shape == (9, 5, 5)
    assert np.all(joined[2:9] <= lattice), "Union erased lattice material"
    assert np.all(joined[:5] <= support), "Union erased support material"
    expected = np.count_nonzero(lattice < 0) + np.count_nonzero(support < 0) - overlap_volume / voxel**3
    assert np.count_nonzero(joined < 0) == expected
    assert overlap_volume == 27 * voxel**3
    for bad_lower in (np.array([-20.0, 0, 0]), np.array([-0.22, 0, 0])):
        try:
            union_fields(lattice, lattice_lower, support, bad_lower, voxel)
        except (ValueError, RuntimeError):
            pass
        else:
            raise AssertionError("Accepted an unaligned or disconnected support")
    cube = vtk.vtkCubeSource()
    cube.SetXLength(2)
    cube.SetYLength(2)
    cube.SetZLength(2)
    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputConnection(cube.GetOutputPort())
    triangles.Update()
    with TemporaryDirectory() as folder:
        output = Path(folder) / "cube.stl"
        write_stl(triangles.GetOutput(), output)
        audit = audit_stl(output, 8.0)
        assert np.isclose(audit["stl_signed_volume_mm3"], 8.0) and len(audit["stl_sha256"]) == 64
    print("PASS: aligned union, preserved materials, positive overlap, rejected bad placements, binary STL audit")


if __name__ == "__main__":
    main()
