"""Run: python test_implicit_surface.py (no ANSYS/CGAL needed)."""
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from export_voxel_surface import implicit_surface
from mesh_obj3 import CORE, BEAM, BALL, voxelize_graph


def main():
    # One cylinder plus endpoint spheres is a capsule with known volume.
    radius, length, voxel = 0.5, 2.0, 0.05
    expected = np.pi * radius**2 * length + 4 / 3 * np.pi * radius**3
    offsets = np.full(3, 0.5)
    xml = f'''<model xmlns="{CORE}" xmlns:b="{BEAM}" xmlns:b2="{BALL}" unit="millimeter">
      <resources><object id="1"><mesh><vertices>
        <vertex x="0" y="0" z="0"/><vertex x="2" y="0" z="0"/>
      </vertices><triangles/><b:beamlattice radius="0.5" cap="butt" b2:ballradius="0.5">
        <b:beams><b:beam v1="0" v2="1"/></b:beams>
        <b2:balls><b2:ball vindex="0"/><b2:ball vindex="1"/></b2:balls>
      </b:beamlattice></mesh></object></resources><build><item objectid="1"/></build></model>'''
    with TemporaryDirectory() as folder:
        source = Path(folder) / "capsule.3mf"
        with ZipFile(source, "w") as archive:
            archive.writestr("3D/3dmodel.model", xml)
        field, lower, _ = voxelize_graph(source, voxel, expected, offsets=offsets, implicit=True)
        x, y, z = np.meshgrid(*[
            lower[i] + (np.arange(field.shape[i]) + offsets[i]) * voxel for i in range(3)
        ], indexing="ij", sparse=True)
        exact = np.sqrt((x - np.clip(x, 0, length))**2 + y**2 + z**2) - radius
        assert np.array_equal(field < 0, exact < 0), "Caps/cylinder/sphere union changed"
        assert np.any((field > 0) & (field < voxel)), "Distance field was binarized"
        surface = implicit_surface(field, lower, voxel, offsets)
        mass = vtk.vtkMassProperties()
        mass.SetInputData(surface)
        mass.Update()
        assert abs(mass.GetVolume() / expected - 1) < 0.02
        points = vtk_to_numpy(surface.GetPoints().GetData()).astype(np.float64)
        faces = vtk_to_numpy(surface.GetPolys().GetConnectivityArray()).reshape(-1, 3)
        a, b, c = (points[faces[:, i]] for i in range(3))
        assert np.einsum("ij,ij->", a, np.cross(b, c)) > 0, "Triangle winding is inward"
        print(f"PASS: closed capsule, outward winding, volume {mass.GetVolume():.6f}/{expected:.6f}")


if __name__ == "__main__":
    main()
