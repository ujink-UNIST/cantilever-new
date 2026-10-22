"""Run: python test_new_modal_cases.py [--cgal for a real CGAL sphere smoke test]."""
import re
import sys
from tempfile import TemporaryDirectory
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

import run
import run_all_modal as modal
from mesh_obj3 import (
    GRAPH_OFFSETS, RADIUS_SCALES, REFERENCE_VOLUMES, SCALES, VOXELS,
    interpolate_field, mesh_implicit_field,
)


def main():
    assert len(run.NEW_CASES) == 6 and len(set(run.STEMS)) == 12
    template = modal.TEMPLATE.read_text(encoding="utf-8")
    modes = int(re.search(r"MODOPT,LANPCG,(\d+)", template).group(1))
    routing = dict((int(i), stem) for i, stem in re.findall(
        r"\*(?:ELSEIF|IF),MESHCASE,EQ,(\d+),THEN\s+/INPUT,\.\./([^,]+),cdb", template,
    ))
    assert routing == {i: f"{stem}_ansys_tet10" for i, stem in enumerate(run.STEMS, 1)}
    selected = [case for case in modal.CASES if run.STEMS[case[0] - 1] in run.NEW_CASES]
    assert len(selected) == 12
    physical_rows = []
    with TemporaryDirectory() as folder, patch.object(modal, "RUNS", Path(folder)):
        for meshcase, matcase, job, mesh, _ in selected:
            stem = run.STEMS[meshcase - 1]
            assert run.SOURCES[stem].is_file()
            assert VOXELS[stem] in (0.2, 0.125) and SCALES[stem] in (1.0, 0.5)
            assert GRAPH_OFFSETS[stem].shape == (3,) and REFERENCE_VOLUMES[stem] > 0
            assert RADIUS_SCALES.get(stem, 1.0) == 1.0, "Uniform radii were scaled twice"
            deck = modal.render(meshcase, matcase, job).read_text(encoding="ascii")
            assert f"MESHCASE={meshcase}\n" in deck and f"MATCASE={matcase}\n" in deck
            assert f"/FILNAME,{job},1" in deck
            physical_rows.extend(modal.scaled_rows(mesh, matcase, [1.0] * modes))
    assert len(physical_rows) == 24 and len({row[0] for row in physical_rows}) == 24
    assert all(len(row) == 3 + modes for row in physical_rows)

    x, y, z = np.indices((4, 5, 6), dtype=float)
    field = (x + 2 * y + 3 * z).astype(np.float32)
    origin, spacing = (-1.0, 2.0, -3.0), 0.2
    for index in ((1.25, 2.5, 3.75), (0.0, 0.0, 0.0), (2.7, 3.2, 4.1)):
        point = [origin[i] + index[i] * spacing for i in range(3)]
        assert np.isclose(interpolate_field(field, point, origin, spacing), index[0] + 2*index[1] + 3*index[2])
    assert interpolate_field(field, (-100, 0, 0), origin, spacing) > 0

    captured = {}
    def generate(domain, **kwargs):
        captured.update(domain=domain, kwargs=kwargs)
        return "mesh"
    fake = SimpleNamespace(DomainBase=object, generate_mesh=generate)
    with patch.dict(sys.modules, pygalmesh=fake):
        mesh, centered_origin = mesh_implicit_field(field, spacing, 0.5)
        point = centered_origin + np.array([1.25, 2.5, 3.75]) * spacing
        assert mesh == "mesh" and np.isclose(captured["domain"].eval(point), 17.5)
        assert captured["domain"].get_bounding_sphere_squared_radius() > centered_origin @ centered_origin
        assert captured["kwargs"]["max_cell_circumradius"] == 0.5
    print("PASS: six models, 12 decks, four scaled materials, no double radius scaling, trilinear CGAL adapter")

    if "--cgal" in sys.argv:
        axis = np.linspace(-1.6, 1.6, 33)
        x, y, z = np.meshgrid(axis, axis, axis, indexing="ij", sparse=True)
        field = (np.sqrt(x*x + y*y + z*z) - 1.0).astype(np.float32)
        mesh, _ = mesh_implicit_field(field, 0.1, 0.2)
        tetra = mesh.get_cells_type("tetra")
        a, b, c, d = (mesh.points[tetra[:, i]] for i in range(4))
        volume = np.abs(np.einsum("ij,ij->i", b-a, np.cross(c-a, d-a))).sum() / 6
        assert len(tetra) > 0 and abs(volume / (4*np.pi/3) - 1) < 0.08
        print(f"PASS: real CGAL sphere, {len(tetra)} tetrahedra, volume {volume:.6f}")


if __name__ == "__main__":
    main()
