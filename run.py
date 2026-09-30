from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CASES = ("modal", "uniform")


def stale(output: Path, source: Path) -> bool:
    return not output.is_file() or output.stat().st_mtime < source.stat().st_mtime


def execute(*args: str) -> None:
    print(f"\n> {' '.join(args)}", flush=True)
    subprocess.run(args, cwd=ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create modal3/uniform3 TET10 meshes and run unconstrained free-free modal cases."
    )
    parser.add_argument("--exe", type=Path, help="Path to ANSYSxxx.exe; auto-detected if omitted")
    parser.add_argument("--np", type=int, default=4, help="MAPDL CPU count (default: 4)")
    args = parser.parse_args()

    for case in CASES:
        source = ROOT / f"{case}_cantilever3.obj"
        volume = ROOT / f"{case}_cantilever3_volume.msh"
        cdb = ROOT / f"{case}_cantilever3_ansys_tet10.cdb"
        if not source.is_file():
            raise FileNotFoundError(source)
        if stale(volume, source):
            execute(sys.executable, "mesh_obj3.py", case)
        else:
            print(f"Reusing {volume.name}")
        if stale(cdb, volume):
            execute(sys.executable, "export_ansys_cdb.py", case)
        else:
            print(f"Reusing {cdb.name}")

    command = [sys.executable, "run_all_modal.py", "--np", str(args.np)]
    if args.exe:
        command += ["--exe", str(args.exe.resolve())]
    execute(*command)
    print(f"\nResults: {ROOT / 'runs' / 'free_free_modal_frequencies_all.csv'}")


if __name__ == "__main__":
    main()
