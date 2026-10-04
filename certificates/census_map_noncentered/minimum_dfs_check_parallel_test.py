"""Adversarial tests for exact prefix-parallel minimum-DFS replay."""
import base64
import copy
import json
import os
import subprocess
import sys
import tempfile
import unittest

import minimum_dfs as D
import minimum_dfs_check_parallel as P


HERE = os.path.dirname(os.path.abspath(__file__))
SERIAL = os.path.join(HERE, "minimum_dfs_check.py")
PARALLEL = os.path.join(HERE, "minimum_dfs_check_parallel.py")

# A small ordinary census forest.  Keeping the independently replayable
# artifact inline makes this test independent of the ignored live census
# outputs in the working directory.
GOLDEN = {
    "spec": {
        "b0": 0.7, "b1": 0.8, "y0": -0.45, "y1": -0.35,
        "seam": 0.137, "delta": 0.001, "dmax": 3.13,
        "n_parts": 64, "minw": 1e-6, "mv_maxw": 0.05,
        "row_mv_maxw": 0.1, "kraw_maxw": 0.1,
        "coord": "relative",
        "split": "relative-u-pi-endpoints-then-widest-midpoint",
        "scope": "census",
    },
    "lo_idx": 0,
    "hi_idx": 1,
    "n_nodes": 154,
    "bits_b64": "55nmeIxH4lEsikXIpFkUi4AAQAA=",
    "stats": {
        "visited_nodes": 61,
        "counts": {
            "MAP_PLAIN": 6, "MAP_CENTERED": 6, "MAP_AUX": 14,
            "MAP_COLLAR": 5, "NEG_T00": 0, "NEG_DET": 0,
            "POS_DET": 0, "UNDECIDED": 0,
        },
        "unresolved_volume": 0.0,
        "seconds": 1.0,
    },
}


def bit_list(document):
    bits = P.decode_bits(document["bits_b64"], document["n_nodes"])
    return [bits[index] for index in range(len(bits))]


def first_leaf(bits):
    reader = P.Reader(bits)

    def walk():
        marker = reader.read()
        if marker:
            return walk()
        start = reader.pos
        return start, reader.read(3)

    return walk()


def write_document(directory, document, name="forest.json"):
    path = os.path.join(directory, name)
    with open(path, "w") as handle:
        json.dump(document, handle, separators=(",", ":"))
    return path


def run(script, artifact, *arguments):
    return subprocess.run(
        [sys.executable, script, artifact] + list(arguments), cwd=HERE,
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        timeout=120)


def output_line(output, prefix):
    return next(line for line in output.splitlines()
                if line.startswith(prefix))


class ParallelReplayTests(unittest.TestCase):
    def test_cached_uncached_and_serial_fresh_replays_agree(self):
        with tempfile.TemporaryDirectory() as directory:
            path = write_document(directory, GOLDEN)
            serial = run(SERIAL, path)
            cached = run(PARALLEL, path, "--depth", "2", "--workers", "2")
            uncached = run(PARALLEL, path, "--depth", "2", "--workers", "2",
                           "--no-cache")
        for result in (serial, cached, uncached):
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertIn("CERTIFICATE VALID AND COMPLETE", result.stdout)
        self.assertEqual(output_line(serial.stdout, "leaves"),
                         output_line(cached.stdout, "leaves"))
        self.assertEqual(output_line(cached.stdout, "leaves"),
                         output_line(uncached.stdout, "leaves"))
        self.assertIn("4 at depth <= 2 over 2 fresh workers", cached.stdout)

    def test_cut_is_deterministic_prefix_free_complete_cover(self):
        scanned = P.scan_document(copy.deepcopy(GOLDEN), 2)
        self.assertEqual(set(scanned["units"]),
                         {"0:00", "0:01", "0:10", "0:11"})
        self.assertEqual(scanned["prefix_nodes"], 3)
        self.assertEqual(scanned["structural_nodes"],
                         GOLDEN["stats"]["visited_nodes"])
        groups = P.balanced_groups(scanned["units"], 3)
        P.validate_assignment(groups, scanned["units"])
        with self.assertRaisesRegex(ValueError, "duplicates"):
            P.validate_assignment(groups + [[groups[0][0]]],
                                  scanned["units"])
        with self.assertRaisesRegex(ValueError, "misses"):
            P.validate_assignment(groups[1:], scanned["units"])

    def test_tree_stats_extra_truncated_and_padding_are_rejected(self):
        bad_stats = copy.deepcopy(GOLDEN)
        bad_stats["stats"]["visited_nodes"] += 1
        with self.assertRaisesRegex(ValueError, "statistics"):
            P.scan_document(bad_stats, 2)

        truncated = copy.deepcopy(GOLDEN)
        truncated["n_nodes"] = 1
        truncated["bits_b64"] = D.encode([1])
        with self.assertRaisesRegex(ValueError, "truncated"):
            P.scan_document(truncated, 2)

        extra = copy.deepcopy(GOLDEN)
        extra_bits = bit_list(extra) + [0, 0, 0, 0]
        extra["n_nodes"] = len(extra_bits)
        extra["bits_b64"] = D.encode(extra_bits)
        with self.assertRaisesRegex(ValueError, "extra tree bits"):
            P.scan_document(extra, 2)

        padding = copy.deepcopy(GOLDEN)
        raw = bytearray(base64.b64decode(padding["bits_b64"]))
        raw[-1] |= 1
        padding["bits_b64"] = base64.b64encode(raw).decode()
        with self.assertRaisesRegex(ValueError, "nonzero padding"):
            P.scan_document(padding, 2)

    def test_undecided_and_false_leaf_claims_are_rejected(self):
        undecided = copy.deepcopy(GOLDEN)
        bits = bit_list(undecided)
        code_start, old_code = first_leaf(bits)
        old_name = D.NAME[old_code]
        bits[code_start:code_start + 3] = [1, 1, 1]
        undecided["bits_b64"] = D.encode(bits)
        undecided["stats"]["counts"][old_name] -= 1
        undecided["stats"]["counts"]["UNDECIDED"] += 1
        with self.assertRaisesRegex(ValueError, "UNDECIDED"):
            P.scan_document(undecided, 2)

        false_claim = copy.deepcopy(GOLDEN)
        bits = bit_list(false_claim)
        code_start, old_code = first_leaf(bits)
        old_name = D.NAME[old_code]
        replacement = "NEG_T00"
        replacement_code = D.CODE[replacement]
        bits[code_start:code_start + 3] = [
            (replacement_code >> 2) & 1,
            (replacement_code >> 1) & 1,
            replacement_code & 1,
        ]
        false_claim["bits_b64"] = D.encode(bits)
        false_claim["stats"]["counts"][old_name] -= 1
        false_claim["stats"]["counts"][replacement] += 1
        with tempfile.TemporaryDirectory() as directory:
            path = write_document(directory, false_claim)
            serial = run(SERIAL, path)
            result = run(PARALLEL, path, "--depth", "2", "--workers", "2")
        self.assertNotEqual(serial.returncode, 0)
        self.assertIn("REPLAY FAILED", serial.stdout)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("PARALLEL REPLAY FAILED", result.stdout)
        self.assertIn("leaf method mismatch", result.stdout)

    def test_worker_reloads_only_the_named_artifact_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            path = write_document(directory, GOLDEN)
            _, digest = P.load_hashed(path)
            changed = copy.deepcopy(GOLDEN)
            changed["stats"]["seconds"] = 2.0
            write_document(directory, changed)
            with self.assertRaisesRegex(ValueError, "changed"):
                P.load_hashed(path, digest)


if __name__ == "__main__":
    unittest.main()
