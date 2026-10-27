"""Run with python test_surface_inputs.py; no CGAL or ANSYS execution needed."""
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile
import struct

import numpy as np

import mesh_obj3 as reader
import run
import run_all_modal as modal


def check():
    triangles = []
    for axis in range(3):
        others = [i for i in range(3) if i != axis]
        for side in (0, 1):
            normal = np.zeros(3); normal[axis] = 2 * side - 1
            for i in range(8):
                for j in range(8):
                    quad = np.zeros((4, 3)); quad[:, axis] = side * .4
                    quad[:, others[0]] = np.array([i, i + 1, i + 1, i]) * .05
                    quad[:, others[1]] = np.array([j, j, j + 1, j + 1]) * .05
                    for indices in ([0, 1, 2], [0, 2, 3]):
                        t = quad[indices]
                        if np.dot(np.cross(t[1] - t[0], t[2] - t[0]), normal) < 0:
                            t = t[[0, 2, 1]]
                        triangles.append(t)
    triangles = np.array(triangles, dtype=np.float32)
    with TemporaryDirectory() as folder:
        root = Path(folder); stl = root / 'cube.stl'; obj = root / 'cube.obj'
        records = np.zeros(len(triangles), dtype=reader.STL_DTYPE)
        records['vertices'] = triangles
        stl.write_bytes(bytes(80) + struct.pack('<I', len(records)) + records.tobytes())
        with obj.open('w') as f:
            f.write(f'# {len(triangles) * 3} vertices, {len(triangles)} faces\n')
            np.savetxt(f, triangles.reshape(-1, 3), fmt='v %.9g %.9g %.9g')
            np.savetxt(f, np.arange(1, len(triangles) * 3 + 1).reshape(-1, 3), fmt='f %d %d %d')
        reader.FACE_BLOCK = 17
        for path in (stl, obj):
            loaded = np.concatenate(list(reader.surface_triangles(path)))
            np.testing.assert_array_equal(loaded, triangles)
        a, lo_a, va = reader.voxelize(stl, .1)
        b, lo_b, vb = reader.voxelize(obj, .1)
        np.testing.assert_array_equal(a, b)
        np.testing.assert_array_equal(lo_a, lo_b)
        assert a.sum() == 64 and abs(va - .4**3) < 1e-8 and va == vb
        damaged = root / 'damaged.stl'; damaged.write_bytes(stl.read_bytes()[:-1])
        try:
            list(reader.surface_triangles(damaged))
        except ValueError:
            pass
        else:
            raise AssertionError('Truncated STL accepted')
        invalid = records.copy(); invalid['vertices'][0, 0, 0] = np.nan
        damaged.write_bytes(bytes(80) + struct.pack('<I', len(records)) + invalid.tobytes())
        try:
            list(reader.surface_triangles(damaged))
        except ValueError:
            pass
        else:
            raise AssertionError('Nonfinite STL accepted')
        saved_root, saved_source = run.ROOT, run.SOURCES['modal-8']
        try:
            run.ROOT = root; run.SOURCES['modal-8'] = root / 'modal-8.stl'
            archive = root / 'modal-8.stl.zip'
            with ZipFile(archive, 'w') as z:
                z.write(stl, 'modal-8.stl')
            assert run.ensure_surface('modal-8').read_bytes() == stl.read_bytes()
            assert not archive.exists()
            run.SOURCES['modal-8'].unlink()
            with ZipFile(archive, 'w') as z:
                z.writestr('../modal-8.stl', stl.read_bytes())
            try:
                run.ensure_surface('modal-8')
            except RuntimeError:
                pass
            else:
                raise AssertionError('Unsafe release ZIP accepted')
        finally:
            run.ROOT = saved_root; run.SOURCES['modal-8'] = saved_source
        assert tuple(reader.SCALES) == run.STEMS and tuple(reader.VOXELS) == run.STEMS
        assert len(modal.CASES) == 2 * len(run.STEMS)
        saved_runs = modal.RUNS
        try:
            modal.RUNS = root
            for case in modal.CASES:
                text = modal.render(*case[:3]).read_text()
                assert f'MESHCASE={case[0]}' in text and f'MATCASE={case[1]}' in text
                assert f'../{run.STEMS[case[0] - 1]}_ansys_tet10,cdb' in text
            text = modal.TEMPLATE.read_text()
            if 'fixed' in modal.TEMPLATE.name:
                assert 'NSEL,S,LOC,X,-1.0E9,0' in text and 'MODOPT,LANPCG,6' in text
            else:
                assert 'NSEL,S,LOC,X' not in text and 'MODOPT,LANPCG,12' in text
        finally:
            modal.RUNS = saved_runs
    print('PASS: STL/OBJ, voxel equivalence, invalid input, safe extraction, all decks and branch boundary conditions')


if __name__ == '__main__':
    check()
