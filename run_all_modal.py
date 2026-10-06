from __future__ import annotations

import argparse
import csv
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUNS = ROOT / "runs"
TEMPLATE = ROOT / "run_free_free_modal.inp"
CASES = [
    (1, 101, "modal3_inconel718", "modal3", "Inconel 718"),
    (1, 102, "modal3_structural_steel", "modal3", "Structural Steel"),
    (1, 103, "modal3_nylon12", "modal3", "Formlabs Nylon 12"),
    (2, 101, "uniform3_inconel718", "uniform3", "Inconel 718"),
    (2, 102, "uniform3_structural_steel", "uniform3", "Structural Steel"),
    (2, 103, "uniform3_nylon12", "uniform3", "Formlabs Nylon 12"),
    (3, 101, "modal4_inconel718", "modal4", "Inconel 718"),
    (3, 102, "modal4_structural_steel", "modal4", "Structural Steel"),
    (3, 103, "modal4_nylon12", "modal4", "Formlabs Nylon 12"),
    (4, 101, "uniform4_inconel718", "uniform4", "Inconel 718"),
    (4, 102, "uniform4_structural_steel", "uniform4", "Structural Steel"),
    (4, 103, "uniform4_nylon12", "uniform4", "Formlabs Nylon 12"),
]


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


def render(meshcase: int, matcase: int, name: str) -> Path:
    text = TEMPLATE.read_text(encoding="utf-8")
    replacements = {
        "/FILNAME,free_free_modal,1": f"/FILNAME,{name},1",
        "\nMESHCASE=1\n": f"\nMESHCASE={meshcase}\n",
        "\nMATCASE=101\n": f"\nMATCASE={matcase}\n",
        "*CFOPEN,modal_frequencies,csv": f"*CFOPEN,{name}_frequencies,csv",
        "PARSAV,ALL,modal_frequency_parameters,parm": f"PARSAV,ALL,{name}_frequencies,parm",
    }
    for old, new in replacements.items():
        if text.count(old) != 1:
            raise RuntimeError(f"Template marker not found exactly once: {old!r}")
        text = text.replace(old, new)
    assert "MODOPT,LANPCG,12" in text and "PCGOPT,1,OFF,NO,OFF,OFF" in text and "MSAVE,ON" in text
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
            if 1 <= mode <= 12:
                by_mode.setdefault(mode, frequency)
    if set(by_mode) != set(range(1, 13)):
        raise RuntimeError(f"Expected modes 1-12 in {path}, found {sorted(by_mode)}")
    return [by_mode[mode] for mode in range(1, 13)]


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
    parser = argparse.ArgumentParser(description="Run unconstrained free-free cantilever3/4 cases with three materials.")
    parser.add_argument("--exe", type=Path, help="Path to ANSYSxxx.exe")
    parser.add_argument("--np", type=int, default=4, help="MAPDL CPU count (default: 4)")
    parser.add_argument("--prepare-only", action="store_true", help="Generate input decks without MAPDL")
    args = parser.parse_args()
    RUNS.mkdir(exist_ok=True)

    inputs = [(case, render(*case[:3])) for case in CASES]
    if args.prepare_only:
        print("Generated:", *(path.name for _, path in inputs), sep="\n  ")
        return

    for stem in ("modal_cantilever3", "uniform_cantilever3", "modal_cantilever4", "uniform_cantilever4"):
        cdb = ROOT / f"{stem}_ansys_tet10.cdb"
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
        if lock_path.exists():
            raise RuntimeError(f"{lock_path.name} exists; stop MAPDL and remove the stale lock")
        if frequency_path.is_file() and result_path.is_file():
            try:
                frequencies = read_frequencies(frequency_path)
            except RuntimeError:
                pass
            else:
                print(f"Reusing completed {name} ...", flush=True)
                results.append([name, mesh, material, *frequencies])
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
        if frequencies[6] and max(map(abs, frequencies[:6])) > abs(frequencies[6]) * 1e-3:
            print(f"WARNING: {name} modes 1-6 are not near-zero relative to mode 7", flush=True)
        results.append([name, mesh, material, *frequencies])
        cleanup_solver_files(name)

    summary = RUNS / "free_free_modal_frequencies_all.csv"
    with summary.open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow([
            "case", "mesh", "material",
            *(f"rigid_body_{i}_hz" for i in range(1, 7)),
            *(f"f{i}_hz" for i in range(1, 7)),
        ])
        writer.writerows(results)
    print(f"Saved {summary.name}")


if __name__ == "__main__":
    main()
