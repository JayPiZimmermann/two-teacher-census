"""Structural tests for the read-only recovered-root embedding audit."""
import copy
from decimal import Decimal
import json
import os
import tempfile
import unittest
from unittest import mock

import minimum_atlas_embedding_audit as A
import minimum_atlas_sweep as S
import minimum_dfs as D
import minimum_dfs_merge as M


def leaf_bits(name):
    code = D.CODE[name]
    return [0, (code >> 2) & 1, (code >> 1) & 1, code & 1]


class EmbeddingAuditTests(unittest.TestCase):
    tile = {"b0": Decimal("0.7"), "b1": Decimal("0.8"),
            "y0": Decimal("-0.15"), "y1": Decimal("-0.1")}
    scope = "census"
    n_parts = 2

    def spec(self):
        return {"b0": 0.7, "b1": 0.8, "y0": -0.15, "y1": -0.1,
                "seam": 0.137, "delta": 0.001, "dmax": 3.13,
                "n_parts": self.n_parts, "minw": 1e-6,
                "mv_maxw": 0.05, "row_mv_maxw": 0.1,
                "kraw_maxw": 0.1, "coord": "relative",
                "split": "relative-u-pi-endpoints-then-widest-midpoint",
                "scope": self.scope}

    def shard(self, lo_idx, hi_idx, verdict="MAP_PLAIN"):
        bits = []
        for _ in range(lo_idx, hi_idx):
            bits.extend(leaf_bits(verdict))
        counts = {name: 0 for name in D.VERDICTS}
        counts[verdict] = hi_idx - lo_idx
        return {"spec": self.spec(), "lo_idx": lo_idx, "hi_idx": hi_idx,
                "n_nodes": len(bits), "bits_b64": M.encode_bits(bits),
                "stats": {"visited_nodes": hi_idx - lo_idx,
                          "counts": counts, "unresolved_volume": 0.0,
                          "seconds": 1.0}}

    def write_shard(self, directory, document, name=None):
        tag = S.tile_id(self.tile, self.scope)
        if name is None:
            name = ("minimum_dfs_%s_p%03d.json" %
                    (tag, document["lo_idx"]))
        path = os.path.join(directory, name)
        with open(path, "w") as handle:
            json.dump(document, handle)
        return path

    def complete(self, directory):
        self.write_shard(directory, self.shard(0, 1))
        self.write_shard(directory, self.shard(1, 2))

    def audit(self, directory):
        with mock.patch.object(S, "HERE", directory):
            return A.audit_recovered_forest(
                self.tile, self.scope, directory, self.n_parts)

    def test_exact_singletons_dry_merge_and_match_sweep_cover(self):
        with tempfile.TemporaryDirectory() as directory:
            self.complete(directory)
            result = self.audit(directory)
        self.assertEqual(result["ranges"], [(0, 1), (1, 2)])
        self.assertEqual(result["merged_nodes"], 2)
        self.assertIsNone(result["discovered"])

    def test_hardlinked_adopted_stage_is_an_ordinary_canonical_shard(self):
        with tempfile.TemporaryDirectory() as directory:
            tag = S.tile_id(self.tile, self.scope)
            stage = os.path.join(directory,
                                 "minimum_dfs_auxstage_test.json")
            with open(stage, "w") as handle:
                json.dump(self.shard(0, 1), handle)
            os.link(stage, os.path.join(
                directory, "minimum_dfs_%s_p000.json" % tag))
            self.write_shard(directory, self.shard(1, 2))
            result = self.audit(directory)
        self.assertEqual(len(result["paths"]), 2)

    def test_missing_root_is_reported_before_merge(self):
        with tempfile.TemporaryDirectory() as directory:
            self.write_shard(directory, self.shard(0, 1))
            with self.assertRaisesRegex(A.AuditError,
                                        r"missing \[\(1, 2\)\]"):
                self.audit(directory)

    def test_strict_scan_rejects_stats_that_loose_discovery_accepts(self):
        with tempfile.TemporaryDirectory() as directory:
            bad = self.shard(0, 1)
            bad["stats"]["visited_nodes"] = 99
            self.write_shard(directory, bad)
            self.write_shard(directory, self.shard(1, 2))
            with self.assertRaisesRegex(A.AuditError,
                                        "stored statistics"):
                self.audit(directory)

    def test_extra_noncanonical_matching_file_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            self.complete(directory)
            tag = S.tile_id(self.tile, self.scope)
            self.write_shard(
                directory, self.shard(0, 1),
                "minimum_dfs_%s_p000_copy.json" % tag)
            with self.assertRaisesRegex(A.AuditError,
                                        "noncanonical shard name"):
                self.audit(directory)

    def test_stale_merged_forest_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            self.complete(directory)
            tag = S.tile_id(self.tile, self.scope)
            stale = M.merge_docs([self.shard(0, 1, "MAP_AUX"),
                                  self.shard(1, 2, "MAP_AUX")])
            with open(os.path.join(
                    directory, "minimum_dfs_%s_merged.json" % tag),
                    "w") as handle:
                json.dump(stale, handle)
            with self.assertRaisesRegex(A.AuditError, "stale or differs"):
                self.audit(directory)

    def test_matching_merged_forest_is_discovered_and_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            left, right = self.shard(0, 1), self.shard(1, 2)
            self.write_shard(directory, left)
            self.write_shard(directory, right)
            tag = S.tile_id(self.tile, self.scope)
            with open(os.path.join(
                    directory, "minimum_dfs_%s_merged.json" % tag),
                    "w") as handle:
                json.dump(M.merge_docs([left, right]), handle)
            result = self.audit(directory)
        self.assertTrue(result["merged_present"])
        self.assertTrue(result["discovered"].endswith("_merged.json"))

    def test_exact_spec_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            bad = copy.deepcopy(self.shard(0, 1))
            bad["spec"]["scope"] = "selected"
            self.write_shard(directory, bad)
            self.write_shard(directory, self.shard(1, 2))
            with self.assertRaisesRegex(A.AuditError,
                                        "tile/spec mismatch"):
                self.audit(directory)


if __name__ == "__main__":
    unittest.main()
