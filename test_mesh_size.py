"""Run: python test_mesh_size.py; fake CGAL/MAPDL, never writes baseline project results."""
import argparse
import csv
import importlib.util
import os
from pathlib import Path
import re
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

import mesh_obj3 as mesher
import run
import run_all_modal as modal


def main():
    assert run.size_suffix(1) == run.size_suffix(1.0) == ""
    assert run.size_suffix(0.5) == "_size0p5"
    assert run.size_suffix(0.25) != run.size_suffix(0.5)
    for bad in (0, -1, "nan", "inf", "-inf", "invalid"):
        try:
            run.positive_size(bad)
        except argparse.ArgumentTypeError:
            pass
        else:
            raise AssertionError(f"Accepted invalid size: {bad}")

    with TemporaryDirectory(prefix="mesh size ") as folder:
        root = Path(folder)
        seen = []
        def fake_mesh():
            return SimpleNamespace(
                points=np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]),
                point_data={}, cell_data={},
                get_cells_type=lambda _: np.array([[0, 1, 2, 3]]),
                write=lambda path, **_: Path(path).write_text("sized mesh"),
            )
        def generate(*args, **kwargs):
            seen.append((args, kwargs))
            return fake_mesh()
        def sample(*args, **kwargs):
            field = np.ones((3, 3, 3), dtype=np.float32 if kwargs.get("implicit") else np.uint8)
            return field, np.zeros(3), 1 / 6
        fake_cgal = SimpleNamespace(DomainBase=object, generate_mesh=generate, generate_from_array=generate)
        with patch.dict(sys.modules, pygalmesh=fake_cgal), patch.object(mesher, "ROOT", root), patch.object(mesher, "voxelize_graph", sample):
            for case in ("modal16_full", "modal16_half", "modal_cantilever4"):
                old = root / f"{case}_volume.msh"
                old.write_text("baseline mesh")
                mesher.mesh_case(case, size=0.5)
                args, criteria = seen[-1]
                scale = mesher.SCALES[case]
                assert criteria["max_cell_circumradius"] == scale * 0.5
                assert criteria["max_radius_surface_delaunay_ball"] == 0.75 * scale * 0.5
                assert criteria["max_facet_distance"] == 0.1 * scale * 0.5
                assert criteria["max_circumradius_edge_ratio"] == 3
                assert (root / f"{case}_size0p5_volume.msh").is_file()
                assert old.read_text() == "baseline mesh"
                if case == "modal_cantilever4":
                    assert args[1] == (0.125,) * 3, "Size changed the field sampling grid"

        reads = []
        def read_mesh(path, **kwargs):
            reads.append(path)
            return fake_mesh()
        spec = importlib.util.spec_from_file_location("sized_exporter_test", Path(__file__).with_name("export_ansys_cdb.py"))
        exporter = importlib.util.module_from_spec(spec)
        with patch.dict(sys.modules, meshio=SimpleNamespace(read=read_mesh)):
            spec.loader.exec_module(exporter)
        exporter.ROOT = root
        old_cdb = root / "modal16_full_ansys_tet10.cdb"
        old_cdb.write_text("baseline CDB")
        exporter.export("modal16_full", size=0.5)
        assert reads[-1].name == "modal16_full_size0p5_volume.msh"
        assert "ET,1,187" in (root / "modal16_full_size0p5_ansys_tet10.cdb").read_text()
        assert old_cdb.read_text() == "baseline CDB"

        source = root / "source.3mf"
        source.write_text("test geometry")
        report = root / "matches.json"
        report.write_text("{}")
        script = root / "mesh_obj3.py"
        script.write_text("new script version")
        os.utime(source, (1000, 1000)); os.utime(report, (1000, 1000))
        os.utime(script, (3000, 3000))
        case = "uniform16_full"
        for suffix in ("_volume.msh", "_ansys_tet10.cdb"):
            p = root / f"{case}{suffix}"
            p.write_text("preserve baseline")
            os.utime(p, (2000, 2000))
        commands = []
        with patch.object(run, "ROOT", root), patch.object(run, "MATCH_REPORT", report), patch.object(run, "ensure_surface", return_value=source), patch.object(run, "execute", side_effect=lambda *args: commands.append(args)):
            with patch.object(sys, "argv", ["run.py", "--size", "0.5", "--cases", case]):
                run.main()
            assert [cmd[1] for cmd in commands] == ["mesh_obj3.py", "export_ansys_cdb.py", "run_all_modal.py"]
            assert all(cmd[cmd.index("--size") + 1] == "0.5" for cmd in commands)
            commands.clear()
            with patch.object(sys, "argv", ["run.py", "--cases", case]):
                run.main()
            assert [cmd[1] for cmd in commands] == ["run_all_modal.py"], "New code alone invalidated the baseline"

        baseline = root / "runs"
        baseline.mkdir()
        selected = [item for item in modal.CASES if run.STEMS[item[0] - 1] == case]
        modes = int(re.search(r"MODOPT,LANPCG,(\d+)", modal.TEMPLATE.read_text()).group(1))
        protected = []
        for item in selected:
            for suffix in (".inp", ".rst", ".out", "_frequencies.csv"):
                p = baseline / f"{item[2]}{suffix}"
                p.write_text("baseline results")
                protected.append((p, p.read_bytes(), p.stat().st_mtime_ns))
        summary_name = "fixed_xle0_modal_frequencies_all.csv" if "fixed" in modal.TEMPLATE.name else "free_free_modal_frequencies_all.csv"
        p = baseline / summary_name
        p.write_text("baseline summary")
        protected.append((p, p.read_bytes(), p.stat().st_mtime_ns))
        with patch.object(run, "ROOT", root), patch.object(modal, "ROOT", root), patch.object(modal, "RUNS", baseline):
            with patch.object(sys, "argv", ["run_all_modal.py", "--size", "0.5", "--cases", case, "--prepare-only"]):
                modal.main()
            sized_runs = run.run_directory(0.5)
            assert sized_runs == baseline / "size0p5"
            cdb = root / f"{case}_size0p5_ansys_tet10.cdb"
            cdb.write_text("sized CDB")
            for item in selected:
                text = (sized_runs / f"{item[2]}.inp").read_text()
                assert f"/INPUT,'{cdb.with_suffix('').as_posix()}',cdb" in text
                (sized_runs / f"{item[2]}.rst").write_text("sized RST")
                frequency_path = sized_runs / f"{item[2]}_frequencies.csv"
                frequency_path.write_text("mode,frequency_hz\n" + "".join(f"{i},{i * 1e-3}\n" for i in range(1, modes + 1)))
            with patch.object(modal, "find_ansys", return_value=root / "ANSYS.exe"), patch.object(modal, "run_mapdl", side_effect=AssertionError("Current sized results should be reused")):
                with patch.object(sys, "argv", ["run_all_modal.py", "--size", "0.5", "--cases", case]):
                    modal.main()
            with (sized_runs / summary_name).open(newline="") as f:
                assert len(list(csv.reader(f))) == 5, "Expected four material rows"
        for p, data, stamp in protected:
            assert p.read_bytes() == data and p.stat().st_mtime_ns == stamp, f"Touched baseline: {p}"
    print("PASS: sized CGAL criteria, unchanged field grid, isolated mesh/CDB/decks/results/summary, sized reuse, baseline bytes and mtimes preserved")


if __name__ == "__main__":
    main()
