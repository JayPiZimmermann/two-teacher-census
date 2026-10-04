"""Structural tests for ``minimum_dfs_merge.py`` (no interval arithmetic)."""
import copy
import unittest

import minimum_dfs as D
import minimum_dfs_merge as M


def spec():
    return {"b0": 0.7, "b1": 0.8, "y0": -0.45, "y1": -0.35,
            "seam": 0.137, "delta": 0.001, "dmax": 3.13,
            "n_parts": 4, "minw": 1e-6, "mv_maxw": 0.05,
            "row_mv_maxw": 0.1, "kraw_maxw": 0.1,
            "coord": "relative",
            "split": "relative-u-pi-endpoints-then-widest-midpoint"}


def shard(lo_idx, hi_idx, bits):
    counts = {name: 0 for name in D.VERDICTS}
    counts["MAP_PLAIN"] = hi_idx - lo_idx
    return {"spec": spec(), "lo_idx": lo_idx, "hi_idx": hi_idx,
            "n_nodes": len(bits), "bits_b64": M.encode_bits(bits),
            "stats": {"visited_nodes": len(bits), "counts": counts,
                      "unresolved_volume": 0.0, "seconds": 1.0}}


class MergeTests(unittest.TestCase):
    def test_round_trip_non_byte_aligned(self):
        bits = [1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1]
        self.assertEqual(M.decode_bits(M.encode_bits(bits), len(bits)), bits)

    def test_merge_sorts_and_concatenates(self):
        left = shard(0, 2, [1, 0, 0])
        right = shard(2, 4, [0, 1, 1, 0, 1])
        merged = M.merge_docs([right, left])
        self.assertEqual((merged["lo_idx"], merged["hi_idx"]), (0, 4))
        self.assertEqual(M.decode_bits(merged["bits_b64"],
                                       merged["n_nodes"]),
                         [1, 0, 0, 0, 1, 1, 0, 1])
        self.assertEqual(merged["stats"]["merged_ranges"], [[0, 2], [2, 4]])
        self.assertEqual(merged["stats"]["counts"]["MAP_PLAIN"], 4)

    def test_gap_and_overlap_rejected(self):
        for bad in ([shard(0, 1, [0]), shard(2, 4, [0])],
                    [shard(0, 3, [0]), shard(2, 4, [0])]):
            with self.assertRaisesRegex(ValueError, "gap/overlap"):
                M.merge_docs(bad)

    def test_spec_mismatch_rejected(self):
        left = shard(0, 2, [0])
        right = shard(2, 4, [0])
        right["spec"] = copy.deepcopy(right["spec"])
        right["spec"]["y1"] = -0.34
        with self.assertRaisesRegex(ValueError, "exact spec"):
            M.merge_docs([left, right])

    def test_undecided_rejected(self):
        left = shard(0, 2, [0])
        right = shard(2, 4, [0])
        right["stats"]["counts"]["UNDECIDED"] = 1
        with self.assertRaisesRegex(ValueError, "undecided"):
            M.merge_docs([left, right])


if __name__ == "__main__":
    unittest.main()
