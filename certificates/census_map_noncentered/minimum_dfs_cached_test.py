"""Equivalence and independent-replay tests for the exact atom cache."""
import copy
import json
import os
import subprocess
import sys
import tempfile
import unittest

import minimum_dfs as D
import minimum_dfs_cached as C


HERE = os.path.dirname(os.path.abspath(__file__))


class CachedDfsTests(unittest.TestCase):
    def test_selected_forest_is_bit_exact_and_replays_uncached(self):
        # This current-format census shard is a compact source of a real F4
        # spec.  Removing ``scope`` asks the default selected question; the
        # uncached and cached searches below start from the same exact spec.
        golden_path = os.path.join(
            HERE, "minimum_dfs_f4probe_census_p000.json")
        with open(golden_path) as handle:
            golden = json.load(handle)
        spec = copy.deepcopy(golden["spec"])
        spec.pop("scope", None)

        plain = D.run(spec, 0, 1, 10000)
        C.install_atom_cache()
        cached = D.run(spec, 0, 1, 10000)

        for field in ("spec", "lo_idx", "hi_idx", "n_nodes", "bits_b64"):
            self.assertEqual(cached[field], plain[field])
        for field in ("visited_nodes", "counts", "unresolved_volume"):
            self.assertEqual(cached["stats"][field], plain["stats"][field])
        info = C.cache_info()
        self.assertGreater(info["hits"], 0)
        self.assertLessEqual(info["size"], info["limit"])

        # The trusted replay path deliberately imports minimum_dfs.py in a
        # fresh process, not this cache wrapper.
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "cached.json")
            with open(path, "w") as handle:
                json.dump(cached, handle, separators=(",", ":"))
            result = subprocess.run(
                [sys.executable, os.path.join(HERE, "minimum_dfs_check.py"),
                 path], cwd=HERE, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, timeout=120)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("CERTIFICATE VALID AND COMPLETE", result.stdout)


if __name__ == "__main__":
    unittest.main()
