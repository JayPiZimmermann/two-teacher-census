"""Structural and independent-replay tests for within-cell frontiers."""
import copy
import base64
import json
import os
import subprocess
import sys
import tempfile
import unittest

import minimum_dfs as D
import minimum_dfs_cached as C
import minimum_dfs_frontier as W
import minimum_frontier_shard as R


HERE = os.path.dirname(os.path.abspath(__file__))


def fake_fragment(spec, path, verdict="MAP_PLAIN"):
    counts = {name: 0 for name in D.VERDICTS}
    counts[verdict] = 1
    code = D.CODE[verdict]
    bits = [0, (code >> 2) & 1, (code >> 1) & 1, code & 1]
    return {"format": W.FORMAT, "spec": spec, "cell_idx": 0,
            "path": path, "n_nodes": len(bits),
            "bits_b64": D.encode(bits),
            "stats": {"visited_nodes": 1, "counts": counts,
                      "unresolved_volume": 0.0, "seconds": 0.1}}


class FrontierTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(os.path.join(
                HERE, "minimum_dfs_f4probe_census_p000.json")) as handle:
            golden = json.load(handle)
        cls.spec = copy.deepcopy(golden["spec"])
        cls.spec.pop("scope", None)

    def test_box_path_is_exact_deterministic_child_chain(self):
        expected = D.root_cells(self.spec)[0]
        for bit in "01011":
            expected = D.child(self.spec, expected, int(bit))
        self.assertEqual(W.box_at_path(self.spec, 0, "01011"), expected)

    def test_fixed_frontier_is_encoded_in_preorder(self):
        fragments = [fake_fragment(self.spec, path)
                     for path in ("00", "01", "10", "11")]
        merged = W.merge_fragments(self.spec, 0, 2, fragments)
        bits = W.decode_bits(merged["bits_b64"], merged["n_nodes"])
        leaf = [0, 0, 0, 0]
        self.assertEqual(bits, [1, 1] + leaf + leaf + [1] + leaf + leaf)
        self.assertEqual(merged["stats"]["visited_nodes"], 7)

    def test_one_frontier_leaf_can_be_replaced_only_by_its_children(self):
        fragments = [fake_fragment(self.spec, path)
                     for path in ("000", "001", "01", "10", "11")]
        merged = W.merge_fragments(self.spec, 0, 2, fragments)
        self.assertEqual(merged["stats"]["frontier_base_depth"], 2)
        self.assertEqual(merged["stats"]["frontier_max_depth"], 3)
        self.assertEqual(merged["stats"]["frontier_fragments"], 5)
        self.assertEqual(merged["stats"]["visited_nodes"], 9)

    def test_prefix_overlaps_are_rejected(self):
        fragments = [fake_fragment(self.spec, path)
                     for path in ("00", "000", "001", "01", "10", "11")]
        with self.assertRaisesRegex(ValueError, "overlap"):
            W.merge_fragments(self.spec, 0, 2, fragments)

    def test_missing_duplicate_and_undecided_frontiers_are_rejected(self):
        fragments = [fake_fragment(self.spec, path)
                     for path in ("00", "01", "10", "11")]
        with self.assertRaisesRegex(ValueError, "incomplete"):
            W.merge_fragments(self.spec, 0, 2, fragments[:-1])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            W.merge_fragments(self.spec, 0, 2,
                              fragments + [fragments[0]])
        bad = copy.deepcopy(fragments)
        bad[0] = fake_fragment(self.spec, "00", verdict="UNDECIDED")
        with self.assertRaisesRegex(ValueError, "undecided"):
            W.merge_fragments(self.spec, 0, 2, bad)

    def test_wrong_metadata_path_and_statistics_are_rejected(self):
        good = [fake_fragment(self.spec, path)
                for path in ("00", "01", "10", "11")]
        for mutate, message in (
                (lambda item: item.update(cell_idx=1), "metadata"),
                (lambda item: item.update(path="001"), "incomplete"),
                (lambda item: item["spec"].update(y1=-0.34), "metadata"),
                (lambda item: item["stats"].update(visited_nodes=2),
                 "node statistics"),
                (lambda item: item["stats"]["counts"].update(MAP_PLAIN=2),
                 "verdict statistics")):
            bad = copy.deepcopy(good)
            mutate(bad[0])
            with self.assertRaisesRegex(ValueError, message):
                W.merge_fragments(self.spec, 0, 2, bad)

    def test_truncated_extra_and_nonzero_padding_bits_are_rejected(self):
        good = fake_fragment(self.spec, "")
        truncated = copy.deepcopy(good)
        truncated["n_nodes"] = 9
        with self.assertRaisesRegex(ValueError, "invalid bit count"):
            W.merge_fragments(self.spec, 0, 0, [truncated])
        extra = copy.deepcopy(good)
        extra["bits_b64"] = base64.b64encode(b"\x00\x00").decode()
        with self.assertRaisesRegex(ValueError, "extra encoded bytes"):
            W.merge_fragments(self.spec, 0, 0, [extra])
        padding = copy.deepcopy(good)
        padding["bits_b64"] = base64.b64encode(b"\x01").decode()
        with self.assertRaisesRegex(ValueError, "nonzero padding"):
            W.merge_fragments(self.spec, 0, 0, [padding])

    def test_resume_plan_reuses_deep_siblings_and_splits_only_failed_path(self):
        available = {"000": "left", "010": "mid-left",
                     "011": "mid-right", "10": "right-left",
                     "11": "right-right"}
        selected, pending = R.plan_frontier(
            2, 5, available, {"00"})
        self.assertEqual(selected,
                         {"000": "left", "010": "mid-left",
                          "011": "mid-right", "10": "right-left",
                          "11": "right-right"})
        self.assertEqual(pending, ["001"])

    def test_structurally_truncated_and_extra_trees_are_rejected(self):
        truncated = fake_fragment(self.spec, "")
        truncated["n_nodes"] = 1
        truncated["bits_b64"] = D.encode([0])
        with self.assertRaisesRegex(ValueError, "truncated"):
            W.merge_fragments(self.spec, 0, 0, [truncated])
        extra = fake_fragment(self.spec, "")
        bits = W.decode_bits(extra["bits_b64"], extra["n_nodes"]) + [0]
        extra["n_nodes"] = len(bits)
        extra["bits_b64"] = D.encode(bits)
        with self.assertRaisesRegex(ValueError, "extra tree bits"):
            W.merge_fragments(self.spec, 0, 0, [extra])

    def test_real_forced_frontier_replays_with_unchanged_uncached_checker(self):
        # Four independently searched subtrees are merged below two forced
        # deterministic branch levels.  The checker process imports only
        # minimum_dfs.py and therefore supplies an independent replay of both
        # the prefix geometry and every reordered map-first leaf claim.
        C.install_atom_cache()
        fragments = [W.make_fragment(self.spec, 0, path, 10000)
                     for path in ("00", "01", "10", "11")]
        merged = W.merge_fragments(self.spec, 0, 2, fragments)
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "frontier.json")
            with open(path, "w") as handle:
                json.dump(merged, handle, separators=(",", ":"))
            result = subprocess.run(
                [sys.executable, os.path.join(HERE, "minimum_dfs_check.py"),
                 path], cwd=HERE, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, timeout=120)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("CERTIFICATE VALID AND COMPLETE", result.stdout)

    def test_real_nonuniform_frontier_replays_uncached(self):
        C.install_atom_cache()
        paths = ("000", "001", "01", "10", "11")
        fragments = [W.make_fragment(self.spec, 0, path, 10000)
                     for path in paths]
        merged = W.merge_fragments(self.spec, 0, 2, fragments)
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "refined-frontier.json")
            with open(path, "w") as handle:
                json.dump(merged, handle, separators=(",", ":"))
            result = subprocess.run(
                [sys.executable, os.path.join(HERE, "minimum_dfs_check.py"),
                 path], cwd=HERE, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, timeout=120)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("CERTIFICATE VALID AND COMPLETE", result.stdout)


if __name__ == "__main__":
    unittest.main()
