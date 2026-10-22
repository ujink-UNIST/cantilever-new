"""Run: python test_match_uniform_3mf.py (checks graph sizing and saved matches)."""
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from csv_to_thick_graph_3mf import package
from match_uniform_3mf import GEOMETRY, PAIRS, TOLERANCE, uniform_model, set_radius
from mesh_obj3 import graph_primitives


def main():
    report = json.loads((GEOMETRY / "new_uniform_volume_matches.json").read_text(encoding="utf-8"))
    with TemporaryDirectory() as folder:
        for name, (modal_name, cells, scale, _) in PAIRS.items():
            model = uniform_model(cells, scale)
            set_radius(model, 0.6)
            path = Path(folder) / "trial.3mf"
            package(model, path)
            points, beams, balls = graph_primitives(path)
            assert np.allclose(points.min(0), 0)
            assert np.allclose(points.max(0), [cells * 5 * scale, 25 * scale, 25 * scale])
            assert len(points) == (3168 if cells == 16 else 4704)
            assert len(beams) == (8952 if cells == 16 else 13368)
            assert np.all(beams[:, 2:] == 0.6) and np.all(balls[:, 1] == 0.6)
            matched = report[name]
            assert abs(matched["relative_volume_difference"]) <= TOLERANCE
            assert hashlib.sha256((GEOMETRY / modal_name).read_bytes()).hexdigest() == matched["modal_sha256"]
            assert hashlib.sha256((GEOMETRY / matched["uniform"]).read_bytes()).hexdigest() == matched["uniform_sha256"]
    print("PASS: three matching cell grids, constant beam/ball radii, unchanged modal sources, output hashes")


if __name__ == "__main__":
    main()
