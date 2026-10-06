from __future__ import annotations

import argparse
import csv
import math
import tempfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

ROOT = Path(__file__).resolve().parent
DEFAULTS = (
    "cantilever16_mo_w0.50_dia1.0_noisland_from_modal_s0.5",
    "cantilever16_uniform_modal_s0.5",
)
RADIUS_OVERRIDES = {"cantilever16_uniform_modal_s0.5": 0.411625}
CORE = "http://schemas.microsoft.com/3dmanufacturing/core/2015/02"
BEAM = "http://schemas.microsoft.com/3dmanufacturing/beamlattice/2017/02"
BALL = "http://schemas.microsoft.com/3dmanufacturing/beamlattice/balls/2020/07"
CONTENT_TYPES = "http://schemas.openxmlformats.org/package/2006/content-types"
RELATIONSHIPS = "http://schemas.openxmlformats.org/package/2006/relationships"

ET.register_namespace("", CORE)
ET.register_namespace("b", BEAM)
ET.register_namespace("b2", BALL)


def rows(folder: Path, name: str) -> list[dict[str, str]]:
    with (folder / name).open(encoding="utf-8-sig", newline="") as file:
        return list(csv.DictReader(file))


def package(model: ET.Element, output: Path) -> None:
    types = ET.Element(f"{{{CONTENT_TYPES}}}Types")
    ET.SubElement(types, f"{{{CONTENT_TYPES}}}Default", {
        "Extension": "rels", "ContentType": "application/vnd.openxmlformats-package.relationships+xml"
    })
    ET.SubElement(types, f"{{{CONTENT_TYPES}}}Default", {
        "Extension": "model", "ContentType": "application/vnd.ms-package.3dmanufacturing-3dmodel+xml"
    })
    relationships = ET.Element(f"{{{RELATIONSHIPS}}}Relationships")
    ET.SubElement(relationships, f"{{{RELATIONSHIPS}}}Relationship", {
        "Target": "/3D/3dmodel.model", "Id": "rel0",
        "Type": "http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel",
    })
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=output.parent, delete=False, suffix=".3mf") as file:
        temporary = Path(file.name)
    try:
        with ZipFile(temporary, "w", ZIP_DEFLATED) as archive:
            archive.writestr("[Content_Types].xml", ET.tostring(types, encoding="utf-8", xml_declaration=True))
            archive.writestr("_rels/.rels", ET.tostring(relationships, encoding="utf-8", xml_declaration=True))
            archive.writestr("3D/3dmodel.model", ET.tostring(model, encoding="utf-8", xml_declaration=True))
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)


def convert(folder: Path, output: Path, radius_override: float | None = None) -> None:
    node_rows = rows(folder, "nodes.csv")
    points = {
        int(row.get("id", row.get("node", ""))): tuple(float(row[axis]) for axis in "xyz")
        for row in node_rows
    }
    vertex = {node: index for index, node in enumerate(points)}
    beams: list[tuple[int, int, float]] = []
    incident: defaultdict[int, list[float]] = defaultdict(list)
    min_length = math.inf
    for row in rows(folder, "struts.csv"):
        n1, n2 = int(row["n1"]), int(row["n2"])
        radius = radius_override if radius_override is not None else float(row["radius_mm"])
        length = math.dist(points[n1], points[n2])
        if radius <= 0 or length <= 0:
            raise ValueError(f"Invalid strut {row.get('id', row.get('strut'))}")
        beams.append((vertex[n1], vertex[n2], radius))
        incident[n1].append(radius)
        incident[n2].append(radius)
        min_length = min(min_length, length)
    node_radius = {
        node: math.sqrt(sum(radius * radius for radius in radii) / len(radii))
        for node, radii in incident.items()
    }

    model = ET.Element(f"{{{CORE}}}model", {
        "unit": "millimeter", "xml:lang": "en-US", "requiredextensions": "b b2"
    })
    ET.SubElement(model, f"{{{CORE}}}metadata", {"name": "Title"}).text = output.stem
    resources = ET.SubElement(model, f"{{{CORE}}}resources")
    obj = ET.SubElement(resources, f"{{{CORE}}}object", {"id": "1", "type": "model"})
    mesh = ET.SubElement(obj, f"{{{CORE}}}mesh")
    vertices = ET.SubElement(mesh, f"{{{CORE}}}vertices")
    for x, y, z in points.values():
        ET.SubElement(vertices, f"{{{CORE}}}vertex", {"x": repr(x), "y": repr(y), "z": repr(z)})
    ET.SubElement(mesh, f"{{{CORE}}}triangles")
    lattice = ET.SubElement(mesh, f"{{{BEAM}}}beamlattice", {
        "minlength": f"{min_length:.15g}", "radius": f"{beams[0][2]:.15g}", "cap": "butt",
        f"{{{BALL}}}ballmode": "mixed",
        f"{{{BALL}}}ballradius": f"{next(iter(node_radius.values())):.15g}",
    })
    beam_elements = ET.SubElement(lattice, f"{{{BEAM}}}beams")
    for v1, v2, radius in beams:
        value = f"{radius:.15g}"
        ET.SubElement(beam_elements, f"{{{BEAM}}}beam", {
            "v1": str(v1), "v2": str(v2), "r1": value, "r2": value,
        })
    balls = ET.SubElement(lattice, f"{{{BALL}}}balls")
    for node, radius in node_radius.items():
        ET.SubElement(balls, f"{{{BALL}}}ball", {
            "vindex": str(vertex[node]), "r": f"{radius:.15g}",
        })
    build = ET.SubElement(model, f"{{{CORE}}}build")
    ET.SubElement(build, f"{{{CORE}}}item", {"objectid": "1"})
    package(model, output)

    with ZipFile(output) as archive:
        checked = ET.fromstring(archive.read("3D/3dmodel.model"))
    assert len(checked.findall(f".//{{{BEAM}}}beam")) == len(beams)
    assert len(checked.findall(f".//{{{BALL}}}ball")) == len(node_radius)
    print(
        f"Wrote {output}: {len(points)} vertices, {len(beams)} beams, "
        f"{len(node_radius)} RMS balls, radius {min(r for _, _, r in beams):.5f}.."
        f"{max(r for _, _, r in beams):.5f} mm"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert nodes/struts CSV to an nTop-importable thick-graph 3MF.")
    parser.add_argument("folders", nargs="*", type=Path, help="CSV folders; defaults to both s0.5 datasets")
    args = parser.parse_args()
    for folder in args.folders or [ROOT / name for name in DEFAULTS]:
        folder = folder.resolve()
        convert(
            folder,
            ROOT / "geometry" / f"{folder.name}_thick_lattice_graph.3mf",
            RADIUS_OVERRIDES.get(folder.name),
        )


if __name__ == "__main__":
    main()
