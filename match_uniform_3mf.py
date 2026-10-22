"""Match the new uniform graphs to modal surface volumes without nTop."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile
import xml.etree.ElementTree as ET

import numpy as np
import vtk

from csv_to_thick_graph_3mf import package
from export_voxel_surface import implicit_surface
from mesh_obj3 import CORE, BEAM, BALL, graph_primitives, voxelize_graph
from run import ROOT

GEOMETRY = ROOT / "geometry"
UNIFORM = GEOMETRY / "cantilever24_uniform5449_from_modal_s0.5.3mf"
PAIRS = {
    "16_full": ("cantilever16_mo_w0.50_dia1.0_sym_iface5_from_modal.3mf", 16, 2.0, 0.2),
    "16_half": ("cantilever16_mo_w0.50_dia1.0_sym_island_iface5_from_modal_s0.5.3mf", 16, 1.0, 0.125),
    "24_half": ("cantilever24_mo_w0.50_dia1.0_sym_island_iface5_from_modal_s0.5.3mf", 24, 1.0, 0.125),
}
OFFSETS = np.array([0.389, 0.811, 0.137])
TOLERANCE = 0.0001  # 0.01% of the modal surface volume at the recorded grid spacing.


def uniform_model(cells: int, scale: float):
    with ZipFile(UNIFORM) as archive:
        model = ET.fromstring(archive.read("3D/3dmodel.model"))
    vertices = model.find(f".//{{{CORE}}}vertices")
    beams = model.find(f".//{{{BEAM}}}beams")
    balls = model.find(f".//{{{BALL}}}balls")
    index = {}
    for old, vertex in enumerate(list(vertices)):
        if float(vertex.get("x")) > cells * 5.0 + 1e-9:
            vertices.remove(vertex)
            continue
        index[old] = len(index)
        for axis in "xyz":
            vertex.set(axis, f"{float(vertex.get(axis)) * scale:.12g}")
    for beam in list(beams):
        v1, v2 = int(beam.get("v1")), int(beam.get("v2"))
        if v1 not in index or v2 not in index:
            beams.remove(beam)
            continue
        beam.set("v1", str(index[v1]))
        beam.set("v2", str(index[v2]))
    for ball in list(balls):
        old = int(ball.get("vindex"))
        if old not in index:
            balls.remove(ball)
        else:
            ball.set("vindex", str(index[old]))
    lattice = model.find(f".//{{{BEAM}}}beamlattice")
    lattice.set("minlength", f"{float(lattice.get('minlength')) * scale:.12g}")
    return model


def set_radius(model, radius: float):
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("Uniform radius must be positive and finite")
    lattice = model.find(f".//{{{BEAM}}}beamlattice")
    value = f"{radius:.9f}"
    lattice.set("radius", value)
    lattice.set(f"{{{BALL}}}ballradius", value)
    for beam in lattice.findall(f".//{{{BEAM}}}beam"):
        beam.set("r1", value)
        beam.set("r2", value)
    for ball in lattice.findall(f".//{{{BALL}}}ball"):
        ball.set("r", value)


def measure(source: Path, voxel: float, radius_scale=1.0, *, validate=False):
    with redirect_stdout(io.StringIO()):
        field, lower, _ = voxelize_graph(source, voxel, None, radius_scale, OFFSETS, implicit=True)
    surface = implicit_surface(field, lower, voxel, OFFSETS, validate=validate)
    mass = vtk.vtkMassProperties()
    mass.SetInputData(surface)
    mass.Update()
    volume = mass.GetVolume()
    if not np.isfinite(volume) or volume <= 0:
        raise RuntimeError(f"Invalid surface volume for {source}")
    return volume, surface.GetNumberOfPolys()


def match(name: str):
    modal_name, cells, scale, voxel = PAIRS[name]
    modal = GEOMETRY / modal_name
    output = GEOMETRY / modal_name.replace("_mo_w0.50_dia1.0_", "_uniform5449_volume_matched_")
    modal_hash = hashlib.sha256(modal.read_bytes()).hexdigest()
    target, modal_faces = measure(modal, voxel, validate=True)
    print(f"{name}: modal surface volume {target:.6f} mm^3, grid {voxel} mm", flush=True)
    model = uniform_model(cells, scale)
    initial_radius = 0.37798 * scale
    set_radius(model, initial_radius)
    with TemporaryDirectory() as folder:
        trial = Path(folder) / "uniform.3mf"
        package(model, trial)
        points, beams, balls = graph_primitives(trial)
        modal_points, _, _ = graph_primitives(modal)
        assert np.allclose(points.min(0), modal_points.min(0))
        assert np.allclose(points.max(0), modal_points.max(0))
        assert len(points) == (3168 if cells == 16 else 4704)
        assert len(beams) == (8952 if cells == 16 else 13368)
        assert len(balls) == len(points)
        lo, hi = initial_radius * 0.5, initial_radius * 2
        vlo, _ = measure(trial, voxel, lo / initial_radius)
        vhi, _ = measure(trial, voxel, hi / initial_radius)
        if not vlo < target < vhi:
            raise RuntimeError("Target volume is outside the radius bracket")
        for iteration in range(25):
            # Safeguarded secant interpolation; bisection bounds retain a valid bracket.
            weight = np.clip((target - vlo) / (vhi - vlo), 0.1, 0.9)
            radius = float(lo + weight * (hi - lo))
            volume, _ = measure(trial, voxel, radius / initial_radius)
            print(f"  {iteration + 1}: radius {radius:.9f} mm, volume {volume:.6f}, delta {(volume/target-1):+.5%}", flush=True)
            if abs(volume / target - 1) <= TOLERANCE:
                break
            if volume < target:
                lo, vlo = radius, volume
            else:
                hi, vhi = radius, volume
        else:
            raise RuntimeError("Uniform radius did not converge")
        set_radius(model, radius)
        package(model, trial)
        volume, uniform_faces = measure(trial, voxel, validate=True)
        if abs(volume / target - 1) > TOLERANCE:
            raise RuntimeError("Serialized uniform did not meet the volume tolerance")
        checked_points, checked_beams, checked_balls = graph_primitives(trial)
        assert np.array_equal(checked_points, points), "Radius fitting changed cell geometry"
        assert np.allclose(checked_beams[:, 2:], round(radius, 9))
        assert np.allclose(checked_balls[:, 1], round(radius, 9))
        if hashlib.sha256(modal.read_bytes()).hexdigest() != modal_hash:
            raise RuntimeError("Modal source changed while matching")
        package(model, output)
    return {
        "modal": modal.name, "modal_sha256": modal_hash,
        "uniform": output.name, "uniform_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "cells": [cells, 5, 5], "cell_size_mm": 5.0 * scale,
        "node_bounds_mm": [points.min(0).tolist(), points.max(0).tolist()],
        "radius_mm": round(radius, 9), "diameter_mm": round(2 * round(radius, 9), 9),
        "grid_mm": voxel, "offsets": OFFSETS.tolist(),
        "modal_surface_volume_mm3": target, "uniform_surface_volume_mm3": volume,
        "relative_volume_difference": volume / target - 1,
        "modal_faces": modal_faces, "uniform_faces": uniform_faces,
        "final_surface_checks": "one component; closed edge-manifold; zero self-intersecting faces",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", choices=PAIRS, nargs="*", default=list(PAIRS))
    args = parser.parse_args()
    report = GEOMETRY / "new_uniform_volume_matches.json"
    results = json.loads(report.read_text(encoding="utf-8")) if report.is_file() else {}
    for name in args.cases:
        results[name] = match(name)
        report.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {report}")


if __name__ == "__main__":
    main()
