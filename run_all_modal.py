from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
from pathlib import Path

from run import STEMS, NEW_CASES, positive_size, size_suffix, run_directory

ROOT = Path(__file__).resolve().parent
RUNS = ROOT / "runs"
TEMPLATE = ROOT / "run_fixed_xle0_modal.inp"
MESHES = ("modal3", "uniform3", "modal4", "uniform4", "modal8", "uniform8") + NEW_CASES
CASES = [
    (meshcase, matcase, f"{mesh}_xle0_nu{nu}", mesh, f"whitened nu={poisson:.2f}")
    for meshcase, mesh in enumerate(MESHES, 1)
    for matcase, nu, poisson in ((201, "030", 0.30), (202, "039", 0.39))
]
MATERIALS = {
    201: (
        ("inconel718", "Inconel 718", 165000.0, 8.220e-9),
        ("structural_steel", "Structural Steel", 200000.0, 7.850e-9),
        ("ms1", "EOS MaragingSteel MS1 (heat treated)", 190000.0, 8.100e-9),
    ),
    202: (("nylon12", "Formlabs Nylon 12", 1850.0, 9.900e-10),),
}


def find_ansys(given: Path | None) -> Path | None:
    candidates = [given] if given else []
    candidates += [Path(v) / "ansys/bin/winx64" for k, v in os.environ.items() if k.startswith("AWP_ROOT")]
    for candidate in candidates:
        if not candidate:
            continue
        if candidate.is_file():
            return candidate.resolve()
        matches = sorted(candidate.glob("ANSYS*.exe"), reverse=True)
        if matches:
            return matches[0].resolve()
    for pattern in (
        "v*/ansys/bin/winx64/ANSYS*.exe",
        "ANSYS Student/v*/ansys/bin/winx64/ANSYS*.exe",
    ):
        matches = sorted(Path("C:/Program Files/ANSYS Inc").glob(pattern), reverse=True)
        if matches:
            return matches[0].resolve()
    return None


def render(meshcase: int, matcase: int, name: str, size: float = 1.0) -> Path:
    text = TEMPLATE.read_text(encoding="utf-8")
    suffix = size_suffix(size)
    if suffix:
        for stem in STEMS:
            marker = f"/INPUT,../{stem}_ansys_tet10,cdb"
            cdb_stem = (ROOT / f"{stem}{suffix}_ansys_tet10").as_posix()
            if text.count(marker) != 1:
                raise RuntimeError(f"Missing or repeated CDB template marker: {marker}")
            text = text.replace(marker, f"/INPUT,'{cdb_stem}',cdb")
    replacements = {
        "/FILNAME,fixed_xle0_modal,1": f"/FILNAME,{name},1",
        "\nMESHCASE=1\n": f"\nMESHCASE={meshcase}\n",
        "\nMATCASE=201\n": f"\nMATCASE={matcase}\n",
        "*CFOPEN,modal_frequencies,csv": f"*CFOPEN,{name}_frequencies,csv",
        "PARSAV,ALL,modal_frequency_parameters,parm": f"PARSAV,ALL,{name}_frequencies,parm",
    }
    for old, new in replacements.items():
        if text.count(old) != 1:
            raise RuntimeError(f"Template marker not found exactly once: {old!r}")
        text = text.replace(old, new)
    assert "MODOPT,LANPCG,6" in text and "PCGOPT,1,OFF,NO,OFF,OFF" in text and "MSAVE,ON" in text
    path = RUNS / f"{name}.inp"
    try:
        if path.is_file() and path.read_text(encoding="ascii") == text:
            return path
        path.write_text(text, encoding="ascii")
    except PermissionError as error:
        raise RuntimeError(f"{path.name} is locked by a running MAPDL process") from error
    return path


def read_frequencies(path: Path) -> list[float]:
    by_mode: dict[int, float] = {}
    with path.open(newline="") as file:
        for row in csv.reader(file):
            try:
                mode, frequency = int(float(row[0])), float(row[1])
            except (ValueError, IndexError):
                continue
            if 1 <= mode <= 6:
                by_mode.setdefault(mode, frequency)
    if set(by_mode) != set(range(1, 7)):
        raise RuntimeError(f"Expected modes 1-6 in {path}, found {sorted(by_mode)}")
    return [by_mode[mode] for mode in range(1, 7)]


def scaled_rows(mesh: str, matcase: int, frequencies: list[float]) -> list[list[str | float]]:
    return [
        [
            f"{mesh}_xle0_{suffix}", mesh, material,
            *(frequency * math.sqrt(elastic_modulus / density) for frequency in frequencies),
        ]
        for suffix, material, elastic_modulus, density in MATERIALS[matcase]
    ]


def cleanup_solver_files(job: str) -> None:
    for suffix in (
        "mode", "full", "full_mat2", "emat", "esav", "db", "rdb", "tri", "mntr",
        "mlv", "PCG", "PCS", "sscr", "page", "stat", "bat", "lock", "err", "log",
    ):
        (RUNS / f"{job}.{suffix}").unlink(missing_ok=True)
    for path in RUNS.glob(f"{job}_[0-9]*.*"):
        path.unlink(missing_ok=True)


def run_mapdl(command: list[str], job: str) -> int:
    process = subprocess.Popen(command, cwd=RUNS)
    try:
        return process.wait()
    except KeyboardInterrupt:
        print(f"\nStopping {job} and its child processes ...", flush=True)
        if os.name == "nt":
            subprocess.run(
                ["taskkill", "/PID", str(process.pid), "/T", "/F"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        else:
            process.terminate()
        process.wait()
        cleanup_solver_files(job)
        raise SystemExit(130)


def main() -> None:
    global RUNS
    parser = argparse.ArgumentParser(description="Run selected x<=0 fixed models using two whitened material solves.")
    parser.add_argument("--exe", type=Path, help="Path to ANSYSxxx.exe")
    parser.add_argument("--np", type=int, default=4, help="MAPDL CPU count (default: 4)")
    parser.add_argument("--prepare-only", action="store_true", help="Generate input decks without MAPDL")
    parser.add_argument("--cases", nargs="+", choices=STEMS, default=NEW_CASES, help="Models to run; defaults to the six new matched graphs")
    parser.add_argument("--size", type=positive_size, default=1.0, help="Use the sized CDBs and isolate results under runs/size.../")
    args = parser.parse_args()
    suffix = size_suffix(args.size)
    RUNS = run_directory(args.size)
    RUNS.mkdir(parents=True, exist_ok=True)

    selected = [case for case in CASES if STEMS[case[0] - 1] in args.cases]
    inputs = [(case, render(*case[:3], size=args.size)) for case in selected]
    if args.prepare_only:
        print("Generated:", *(path.name for _, path in inputs), sep="\n  ")
        return

    for stem in args.cases:
        cdb = ROOT / f"{stem}{suffix}_ansys_tet10.cdb"
        if not cdb.is_file():
            raise FileNotFoundError(cdb)

    exe = find_ansys(args.exe)
    if not exe:
        raise SystemExit("ANSYS executable not found. Use --exe C:/.../ANSYS252.exe")

    results = []
    for (meshcase, matcase, name, mesh, material), input_path in inputs:
        output_path = RUNS / f"{name}.out"
        frequency_path = RUNS / f"{name}_frequencies.csv"
        result_path = RUNS / f"{name}.rst"
        lock_path = RUNS / f"{name}.lock"
        cdb_path = ROOT / f"{STEMS[meshcase - 1]}{suffix}_ansys_tet10.cdb"
        if lock_path.exists():
            raise RuntimeError(f"{lock_path.name} exists; stop MAPDL and remove the stale lock")
        result_is_current = (
            frequency_path.is_file() and result_path.is_file()
            and min(frequency_path.stat().st_mtime, result_path.stat().st_mtime)
            >= max(input_path.stat().st_mtime, cdb_path.stat().st_mtime)
        )
        if result_is_current:
            try:
                frequencies = read_frequencies(frequency_path)
            except RuntimeError:
                pass
            else:
                print(f"Reusing completed {name} ...", flush=True)
                results.extend(scaled_rows(mesh, matcase, frequencies))
                cleanup_solver_files(name)
                continue

        cleanup_solver_files(name)
        frequency_path.unlink(missing_ok=True)
        result_path.unlink(missing_ok=True)
        print(f"Running {name} ... progress is in {output_path.name}", flush=True)
        returncode = run_mapdl(
            [str(exe), "-b", "-smp", "-np", str(args.np), "-j", name, "-i", input_path.name, "-o", output_path.name],
            name,
        )
        if returncode or not frequency_path.is_file():
            cleanup_solver_files(name)
            result_path.unlink(missing_ok=True)
            raise RuntimeError(f"{name} failed; inspect {output_path.name}")
        frequencies = read_frequencies(frequency_path)
        results.extend(scaled_rows(mesh, matcase, frequencies))
        cleanup_solver_files(name)

    summary = RUNS / "fixed_xle0_modal_frequencies_all.csv"
    with summary.open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["case", "mesh", "material", *(f"f{i}_hz" for i in range(1, 7))])
        writer.writerows(results)
    print(f"Saved {summary.name}")


if __name__ == "__main__":
    main()
