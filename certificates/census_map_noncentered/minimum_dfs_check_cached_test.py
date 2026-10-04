"""Fresh-process equivalence tests for the exact-cache replay entry point."""
import copy
import json
import os
import subprocess
import sys
import tempfile
import unittest


HERE = os.path.dirname(os.path.abspath(__file__))
PLAIN = os.path.join(HERE, "minimum_dfs_check.py")
CACHED = os.path.join(HERE, "minimum_dfs_check_cached.py")


def replay(script, artifact):
    return subprocess.run(
        [sys.executable, script, artifact], cwd=HERE, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)


class CachedReplayTests(unittest.TestCase):
    def test_real_forest_has_identical_uncached_and_cached_replay(self):
        artifact = os.path.join(
            HERE, "minimum_dfs_f4probe_census_p000.json")
        plain = replay(PLAIN, artifact)
        cached = replay(CACHED, artifact)
        self.assertEqual((cached.returncode, cached.stdout),
                         (plain.returncode, plain.stdout))
        self.assertEqual(cached.returncode, 0, cached.stdout)
        self.assertIn("CERTIFICATE VALID AND COMPLETE", cached.stdout)

    def test_cached_entry_point_rejects_exactly_what_uncached_rejects(self):
        source = os.path.join(
            HERE, "minimum_dfs_f4probe_census_p000.json")
        with open(source) as handle:
            invalid = copy.deepcopy(json.load(handle))
        invalid["lo_idx"] = -1
        with tempfile.TemporaryDirectory() as directory:
            artifact = os.path.join(directory, "invalid.json")
            with open(artifact, "w") as handle:
                json.dump(invalid, handle, separators=(",", ":"))
            plain = replay(PLAIN, artifact)
            cached = replay(CACHED, artifact)
        self.assertEqual((cached.returncode, cached.stdout),
                         (plain.returncode, plain.stdout))
        self.assertNotEqual(cached.returncode, 0)
        self.assertIn("REPLAY FAILED", cached.stdout)


if __name__ == "__main__":
    unittest.main()
