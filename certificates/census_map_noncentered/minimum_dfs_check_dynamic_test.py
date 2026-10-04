"""Adversarial tests for dynamically queued exact minimum-DFS replay."""
import copy
import os
import tempfile
import unittest
from unittest import mock

import minimum_atlas_check as A
import minimum_atlas_sweep as S
import minimum_dfs as D
import minimum_dfs_check_dynamic as Q
import minimum_dfs_check_parallel as P
import minimum_frontier_shard as F
from minimum_dfs_check_parallel_test import (
    GOLDEN, SERIAL, bit_list, first_leaf, output_line, run, write_document)


HERE = os.path.dirname(os.path.abspath(__file__))
DYNAMIC = os.path.join(HERE, "minimum_dfs_check_dynamic.py")


class DynamicReplayTests(unittest.TestCase):
    def tearDown(self):
        Q._WORKER_SCANNED = None
        Q._WORKER_DIGEST = None

    def test_cached_uncached_and_serial_fresh_replays_agree(self):
        with tempfile.TemporaryDirectory() as directory:
            path = write_document(directory, GOLDEN)
            serial = run(SERIAL, path)
            cached = run(DYNAMIC, path, "--depth", "2", "--workers", "2")
            uncached = run(DYNAMIC, path, "--depth", "2", "--workers", "2",
                           "--no-cache")
        for result in (serial, cached, uncached):
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertIn("CERTIFICATE VALID AND COMPLETE", result.stdout)
        self.assertEqual(output_line(serial.stdout, "leaves"),
                         output_line(cached.stdout, "leaves"))
        self.assertEqual(output_line(cached.stdout, "leaves"),
                         output_line(uncached.stdout, "leaves"))
        self.assertIn("4 at depth <= 2 in 4 queued groups", cached.stdout)

    def test_false_leaf_claim_is_rejected_by_dynamic_replay(self):
        document = copy.deepcopy(GOLDEN)
        bits = bit_list(document)
        code_start, old_code = first_leaf(bits)
        old_name = D.NAME[old_code]
        replacement = "NEG_T00"
        code = D.CODE[replacement]
        bits[code_start:code_start + 3] = [
            (code >> 2) & 1, (code >> 1) & 1, code & 1]
        document["bits_b64"] = D.encode(bits)
        document["stats"]["counts"][old_name] -= 1
        document["stats"]["counts"][replacement] += 1
        with tempfile.TemporaryDirectory() as directory:
            path = write_document(directory, document)
            result = run(
                DYNAMIC, path, "--depth", "2", "--workers", "2")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("DYNAMIC REPLAY FAILED", result.stdout)
        self.assertIn("leaf method mismatch", result.stdout)

    def test_queue_is_longest_first_and_prefix_exact(self):
        units = {"small": {"nodes": 1}, "largest": {"nodes": 11},
                 "middle": {"nodes": 5}, "tie": {"nodes": 5}}
        groups = Q.queued_groups(units, batch_size=2)
        self.assertEqual(groups,
                         [["largest", "middle"], ["tie", "small"]])
        P.validate_assignment(groups, units)
        for bad in (0, -1, 1.5):
            with self.subTest(batch=bad), self.assertRaises(ValueError):
                Q.queued_groups(units, bad)

    @mock.patch.object(Q.C, "install_atom_cache")
    @mock.patch.object(Q.P, "replay_unit")
    @mock.patch.object(Q.P, "scan_document")
    @mock.patch.object(Q.P, "load_hashed")
    def test_worker_loads_and_scans_once_then_handles_multiple_jobs(
            self, load_hashed, scan_document, replay_unit, install_cache):
        counts = P.empty_counts()
        counts["MAP_PLAIN"] = 1
        scanned = {"units": {"a": object(), "b": object()}}
        load_hashed.return_value = ({"doc": True}, "abc")
        scan_document.return_value = scanned
        replay_unit.return_value = (1, counts, 0.25)

        Q.initialize_worker("forest.json", "abc", 10, True)
        first = Q.replay_queued_group(["a"])
        second = Q.replay_queued_group(["b"])

        load_hashed.assert_called_once_with("forest.json", "abc")
        scan_document.assert_called_once_with({"doc": True}, 10)
        install_cache.assert_called_once_with()
        self.assertEqual(replay_unit.call_count, 2)
        self.assertEqual(first["digest"], "abc")
        self.assertEqual(second["digest"], "abc")

    def test_recomposition_rejects_wrong_hash_duplicate_and_bad_totals(self):
        scanned = P.scan_document(copy.deepcopy(GOLDEN), 2)
        digest = "sha256-placeholder"
        results = []
        for key, unit in scanned["units"].items():
            results.append({
                "keys": [key], "nodes": unit["nodes"],
                "counts": dict(unit["counts"]),
                "leaf_volume": D.volume(unit["box"]),
                "pid": 1, "digest": digest})
        totals = Q.recompose(scanned, results, digest)
        self.assertEqual(totals["nodes"],
                         GOLDEN["stats"]["visited_nodes"])
        self.assertEqual(totals["counts"], GOLDEN["stats"]["counts"])

        wrong_hash = copy.deepcopy(results)
        wrong_hash[0]["digest"] = "wrong"
        with self.assertRaisesRegex(ValueError, "wrong artifact hash"):
            Q.recompose(scanned, wrong_hash, digest)
        duplicate = copy.deepcopy(results)
        duplicate[-1]["keys"] = duplicate[0]["keys"]
        with self.assertRaisesRegex(ValueError, "duplicates"):
            Q.recompose(scanned, duplicate, digest)
        bad_nodes = copy.deepcopy(results)
        bad_nodes[0]["nodes"] += 1
        with self.assertRaisesRegex(ValueError, "do not recompose"):
            Q.recompose(scanned, bad_nodes, digest)

    def test_dynamic_checker_is_explicitly_allowlisted(self):
        name = "minimum_dfs_check_dynamic.py"
        self.assertIn(name, S.DFS_CHECKERS)
        self.assertIn(name, A.DFS_CHECKERS)
        self.assertIn(name, F.CHECKERS)
        self.assertEqual(S.validate_dfs_checker(name), name)


if __name__ == "__main__":
    unittest.main()
