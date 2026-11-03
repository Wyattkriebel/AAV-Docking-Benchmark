import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('comparison', ROOT / 'pilot/pilot-script.py')
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


class ComparisonTests(unittest.TestCase):
    def test_geometric_controls(self):
        native = [(0, 0, 0), (1, 2, 3)]
        moved = [(x + 3, y + 4, z) for x, y, z in native]
        self.assertEqual(comparison.rmsd(native, native), 0)
        self.assertEqual(comparison.rmsd(moved, native), 5)
        transform = lambda points: [(-y + 7, x - 2, z + 8) for x, y, z in points]
        self.assertEqual(comparison.rmsd(transform(moved), transform(native)), 5)

    def test_atom_order_does_not_change_pose(self):
        inputs = comparison.ROOT / 'inputs'
        names = comparison.read_json(inputs / 'reference.json')['names']
        text = (inputs / 'holo_complete_104729.pdbqt').read_text()
        shuffled = []
        for model in text.split('ENDMDL')[:-1]:
            lines = model.splitlines()
            atoms = [line for line in lines if line.startswith(('ATOM  ', 'HETATM'))]
            other = [line for line in lines if not line.startswith(('ATOM  ', 'HETATM'))]
            shuffled.extend([*other, *reversed(atoms), 'ENDMDL'])
        self.assertEqual(comparison.pose_models(text, names), comparison.pose_models('\n'.join(shuffled), names))

    def test_invalid_models_fail(self):
        inputs = comparison.ROOT / 'inputs'
        names = comparison.read_json(inputs / 'reference.json')['names']
        text = (inputs / 'holo_complete_104729.pdbqt').read_text().split('ENDMDL')[0] + 'ENDMDL\n'
        atom = next(line for line in text.splitlines() if line.startswith('ATOM  ') and line[12:16].strip() == names[0])
        variants = [text.replace(atom + '\n', ''), text.replace(atom, atom + '\n' + atom),
                    text.replace('ENDMDL', ''), text.replace('MODEL 1', 'MODEL 2'),
                    text.replace(atom, atom[:30] + '     nan' + atom[38:]), '']
        for invalid in variants:
            with self.subTest(invalid=invalid[:30]), self.assertRaises(ValueError):
                comparison.pose_models(invalid, names)

    def test_missing_runs_and_changed_inputs_fail(self):
        inputs = comparison.ROOT / 'inputs'
        with tempfile.TemporaryDirectory() as directory:
            copy = Path(directory)
            for path in inputs.iterdir():
                (copy / path.name).write_bytes(path.read_bytes())
            runs = comparison.read_json(copy / 'runs.json')[:-1]
            (copy / 'runs.json').write_text(json.dumps(runs))
            with self.assertRaisesRegex(ValueError, 'Changed archived input'):
                comparison.analyze(copy)
            sources = comparison.read_json(copy / 'SOURCES.json')
            sources['files']['runs.json']['sha256'] = comparison.sha256(copy / 'runs.json')
            (copy / 'SOURCES.json').write_text(json.dumps(sources))
            with self.assertRaisesRegex(ValueError, 'all six matched'):
                comparison.analyze(copy)

    def test_archive_replays_from_another_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            run = subprocess.run([sys.executable, str(ROOT / 'pilot/pilot-script.py'), '--check'],
                                 cwd=directory, capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)


if __name__ == '__main__':
    unittest.main()
