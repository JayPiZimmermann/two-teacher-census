"""Structural tests for singleton-shard escalation scheduling."""
import collections
import unittest
from unittest import mock

import minimum_shard_escalate as E


def tile(tag):
    return {"b0": tag, "b1": tag + 0.1, "y0": -0.9, "y1": -0.85}


class EscalateTests(unittest.TestCase):
    def test_cached_generator_is_default_only_for_refined_partitions(self):
        self.assertEqual(E.generator_script(64), "minimum_dfs.py")
        self.assertEqual(E.generator_script(512), "minimum_dfs_cached.py")
        self.assertEqual(
            E.generator_script(64, "minimum_dfs_cached.py"),
            "minimum_dfs_cached.py")
        self.assertEqual(
            E.generator_script(512, "minimum_dfs.py"), "minimum_dfs.py")

    def test_unknown_generator_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "unknown DFS script"):
            E.generator_script(512, "other.py")

    def test_refined_partition_has_disjoint_artifact_names(self):
        item = tile(0.1)
        self.assertNotEqual(E.tile_tag(item, 64), E.tile_tag(item, 512))
        self.assertNotEqual(E.state_path(64), E.state_path(512))
        self.assertNotEqual(E.log_dir(64), E.log_dir(512))

    @mock.patch.object(E.F, "usable_shards",
                       return_value=[(0, 128), (128, 256)])
    def test_refined_missing_range_uses_whole_new_partition(self, usable):
        item = tile(0.1)
        self.assertEqual(
            E.F.missing_ranges(item, n_parts=512,
                               tag=E.tile_tag(item, 512)),
            [(256, 512)])
        usable.assert_called_once_with(
            item, "selected", 512, E.tile_tag(item, 512))

    @mock.patch.object(E, "missing_cells", return_value=[1, 3, 7])
    def test_refined_cells_run_high_angle_first(self, missing):
        item = tile(0.1)
        self.assertEqual(E.scheduled_cells(item, 512), [7, 3, 1])
        missing.assert_called_once_with(item, 512)

    @mock.patch.object(E, "missing_cells", return_value=[1, 3, 7])
    def test_default_partition_keeps_recovery_order(self, missing):
        item = tile(0.1)
        self.assertEqual(E.scheduled_cells(item, 64), [1, 3, 7])
        missing.assert_called_once_with(item, 64)

    @mock.patch.object(E.F, "tile_id", side_effect=lambda item: str(item["b0"]))
    @mock.patch.object(E, "missing_cells",
                       side_effect=lambda item, n_parts=64: []
                       if item["b0"] == 0.1
                       else [42])
    def test_complete_unjournaled_tile_is_finished_on_recovery(
            self, missing, tile_id):
        complete = tile(0.1)
        incomplete = tile(0.2)
        outstanding = collections.Counter({"0.2": 1})
        self.assertEqual(
            E.ready_to_finish([complete, incomplete], outstanding),
            [complete])

    @mock.patch.object(E.F, "replay_boundary")
    @mock.patch.object(E.F, "generate_boundary", return_value="face.json")
    @mock.patch.object(E.F, "discover", return_value=({}, {}))
    @mock.patch.object(E.F, "replay_regularity")
    @mock.patch.object(E, "merge_tile", return_value="merged.json")
    def test_finish_tile_forwards_allowlisted_checker(
            self, merge, replay, discover, generate, replay_boundary):
        outcome = E.finish_tile(
            tile(0.1), 64, "minimum_dfs_check_parallel.py")
        self.assertEqual(outcome["kind"], "passed")
        replay.assert_called_once_with(
            "merged.json", "minimum_dfs_check_parallel.py")


if __name__ == "__main__":
    unittest.main()
