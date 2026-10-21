from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parent
STEMS = (
    "modal_cantilever3",
    "uniform_cantilever3",
    "modal_cantilever4",
    "uniform_cantilever4",
    "modal-8",
    "uniform-8",
)
SOURCES = {
    "modal_cantilever3": ROOT / "geometry/modal_cantilever3_corrected.3mf",
    "uniform_cantilever3": ROOT / "geometry/uniform_cantilever3_volume_matched.3mf",
    "modal_cantilever4": ROOT / "geometry/cantilever16_mo_w0.50_dia1.0_noisland_from_modal_s0.5_thick_lattice_graph.3mf",
    "uniform_cantilever4": ROOT / "geometry/cantilever16_uniform_modal_s0.5_thick_lattice_graph.3mf",
    "modal-8": ROOT / "modal-8.stl",
    "uniform-8": ROOT / "uniform-8.stl",
}


def stale(output: Path, *sources: Path) -> bool:
    return not output.is_file() or any(output.stat().st_mtime < source.stat().st_mtime for source in sources)


def execute(*args: str) -> None:
    print(f"\n> {' '.join(args)}", flush=True)
    subprocess.run(args, cwd=ROOT, check=True)


def ensure_surface(stem: str) -> Path:
    source = SOURCES[stem]
    if source.is_file():
        return source
    if source.suffix.lower() == ".3mf":
        raise FileNotFoundError(source)
    archive = ROOT / f"{source.name}.zip"
    if not archive.is_file():
        url = f"https://github.com/ujink-UNIST/cantilever-new/releases/download/v1.0/{archive.name}"
        print(f"Downloading {archive.name} ...", flush=True)
        partial = archive.with_suffix(archive.suffix + ".part")
        try:
            urllib.request.urlretrieve(url, partial)
            partial.replace(archive)
        finally:
            partial.unlink(missing_ok=True)
    with ZipFile(archive) as zipped:
        if zipped.namelist() != [source.name]:
            raise RuntimeError(f"Unexpected contents in {archive.name}")
        with tempfile.TemporaryDirectory(dir=ROOT) as folder:
            staged = Path(folder) / source.name
            with zipped.open(source.name) as incoming, staged.open("xb") as outgoing:
                shutil.copyfileobj(incoming, outgoing, 8 * 1024 * 1024)
            os.link(staged, source)  # Publish a complete file without replacing user geometry.
    archive.unlink()
    return source


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create TET10 meshes directly from 3MF graphs/STL and run free-free modal cases."
    )
    parser.add_argument("--exe", type=Path, help="Path to ANSYSxxx.exe; auto-detected if omitted")
    parser.add_argument("--np", type=int, default=4, help="MAPDL CPU count (default: 4)")
    parser.add_argument("--cases", nargs="+", choices=STEMS, default=STEMS, help="Models to run; defaults to all")
    args = parser.parse_args()

    for stem in args.cases:
        source = ensure_surface(stem)
        volume = ROOT / f"{stem}_volume.msh"
        cdb = ROOT / f"{stem}_ansys_tet10.cdb"
        mesh_inputs = (source, ROOT / "mesh_obj3.py") if source.suffix.lower() == ".3mf" else (source,)
        if stale(volume, *mesh_inputs):
            execute(sys.executable, "mesh_obj3.py", stem)
        else:
            print(f"Reusing {volume.name}")
        if stale(cdb, volume):
            execute(sys.executable, "export_ansys_cdb.py", stem)
        else:
            print(f"Reusing {cdb.name}")

    command = [sys.executable, "run_all_modal.py", "--np", str(args.np), "--cases", *args.cases]
    if args.exe:
        command += ["--exe", str(args.exe.resolve())]
    execute(*command)
    print(f"\nResults: {ROOT / 'runs' / 'free_free_modal_frequencies_all.csv'}")


if __name__ == "__main__":
    main()
