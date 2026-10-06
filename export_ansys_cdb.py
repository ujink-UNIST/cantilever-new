from __future__ import annotations

import argparse
from pathlib import Path

import meshio
import numpy as np

ROOT = Path(__file__).resolve().parent
STEMS = (
    "modal_cantilever3",
    "uniform_cantilever3",
    "modal_cantilever4",
    "uniform_cantilever4",
)
EDGE_ORDER = np.array(((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)))


def export(case: str) -> None:
    stem = case
    source = ROOT / f"{stem}_volume.msh"
    output = ROOT / f"{stem}_ansys_tet10.cdb"
    mesh = meshio.read(source, file_format="gmsh")
    points = mesh.points
    corners = mesh.get_cells_type("tetra")
    if not len(corners):
        raise RuntimeError(f"No tetrahedra in {source.name}")

    edges = corners[:, EDGE_ORDER].reshape(-1, 2)
    edges.sort(axis=1)
    unique_edges, inverse = np.unique(edges, axis=0, return_inverse=True)
    midpoints = (points[unique_edges[:, 0]] + points[unique_edges[:, 1]]) / 2
    midpoint_ids = inverse.reshape(-1, 6) + len(points) + 1
    tet10 = np.hstack((corners + 1, midpoint_ids))
    total_nodes = len(points) + len(midpoints)
    assert tet10.shape[1] == 10 and tet10.max() == total_nodes

    with output.open("w", newline="\n", buffering=1024 * 1024) as out:
        out.write(
            "/COM,ANSYS RELEASE 2023 R2 BUILD 23.2\n"
            f"/COM,{stem} TET10; coordinates in mm\n"
            "/UNITS,MPA\n/PREP7\nET,1,187\n"
        )
        out.write(f"NBLOCK,6,SOLID,{total_nodes:9d},{total_nodes:9d}\n(3i9,6e21.13e3)\n")
        node_id = 1
        for node_block in (points, midpoints):
            for start in range(0, len(node_block), 10_000):
                block = node_block[start : start + 10_000]
                out.writelines(
                    f"{i:9d}{0:9d}{0:9d}{x:21.13E}{y:21.13E}{z:21.13E}\n"
                    for i, (x, y, z) in enumerate(block, node_id)
                )
                node_id += len(block)
        out.write("N,R5.3,LOC,       -1,\n")

        out.write(f"EBLOCK,19,SOLID,{len(tet10):9d},{len(tet10):9d}\n(19i9)\n")
        for start in range(0, len(tet10), 10_000):
            lines = []
            for element_id, nodes in enumerate(tet10[start : start + 10_000], start + 1):
                lines.append(
                    f"{1:9d}{1:9d}{1:9d}{1:9d}{0:9d}{0:9d}{0:9d}{0:9d}{10:9d}{0:9d}{element_id:9d}"
                    + "".join(f"{node:9d}" for node in nodes[:8])
                    + "\n"
                )
                lines.append("".join(f"{node:9d}" for node in nodes[8:]) + "\n")
            out.writelines(lines)
        out.write(f"{-1:9d}\nFINISH\n")

    print(f"Wrote {output.name}: {total_nodes} nodes, {len(tet10)} SOLID187 TET10 elements")


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert a cantilever tetra mesh to ANSYS SOLID187 CDB.")
    parser.add_argument("case", choices=STEMS)
    args = parser.parse_args()
    export(args.case)


if __name__ == "__main__":
    main()
