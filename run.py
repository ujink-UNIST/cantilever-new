from __future__ import annotations

import argparse
import subprocess
import sys
import urllib.request
from pathlib import Path
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parent
STEMS = (
    "modal_cantilever3",
    "uniform_cantilever3",
    "modal_cantilever4",
    "uniform_cantilever4",
)


def stale(output: Path, source: Path) -> bool:
    return not output.is_file() or output.stat().st_mtime < source.stat().st_mtime


def execute(*args: str) -> None:
    print(f"\n> {' '.join(args)}", flush=True)
    subprocess.run(args, cwd=ROOT, check=True)


def ensure_obj(stem: str) -> Path:
    source = ROOT / f"{stem}.obj"
    if source.is_file():
        return source
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
        zipped.extract(source.name, ROOT)
    archive.unlink()
    return source


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create cantilever3/4 TET10 meshes and run x<=0 fixed modal cases."
    )
    parser.add_argument("--exe", type=Path, help="Path to ANSYSxxx.exe; auto-detected if omitted")
    parser.add_argument("--np", type=int, default=4, help="MAPDL CPU count (default: 4)")
    args = parser.parse_args()

    for stem in STEMS:
        source = ensure_obj(stem)
        volume = ROOT / f"{stem}_volume.msh"
        cdb = ROOT / f"{stem}_ansys_tet10.cdb"
        if stale(volume, source):
            execute(sys.executable, "mesh_obj3.py", stem)
        else:
            print(f"Reusing {volume.name}")
        if stale(cdb, volume):
            execute(sys.executable, "export_ansys_cdb.py", stem)
        else:
            print(f"Reusing {cdb.name}")

    command = [sys.executable, "run_all_modal.py", "--np", str(args.np)]
    if args.exe:
        command += ["--exe", str(args.exe.resolve())]
    execute(*command)
    print(f"\nResults: {ROOT / 'runs' / 'fixed_xle0_modal_frequencies_all.csv'}")


if __name__ == "__main__":
    main()
