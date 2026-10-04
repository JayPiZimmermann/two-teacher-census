"""Structural tests for exact atlas scheduling cutoffs."""
import json
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from decimal import Decimal
from unittest import mock

import minimum_atlas_check as C
import minimum_atlas_sweep as S


class DfsScriptTests(unittest.TestCase):
    tile = {"b0": Decimal("0.7"), "b1": Decimal("0.8"),
            "y0": Decimal("-0.45"), "y1": Decimal("-0.35")}

    def test_only_allowed_in_directory_basenames_are_accepted(self):
        for script in S.DFS_SCRIPTS:
            self.assertEqual(S.validate_dfs_script(script), script)
        for script in ("../minimum_dfs_cached.py",
                       "/tmp/minimum_dfs_cached.py",
                       "minimum_boundary_dfs.py", "other.py", ""):
            with self.assertRaisesRegex(ValueError, "DFS script"):
                S.validate_dfs_script(script)

    def test_only_allowed_in_directory_checker_basenames_are_accepted(self):
        for script in S.DFS_CHECKERS:
            self.assertEqual(S.validate_dfs_checker(script), script)
        for script in ("../minimum_dfs_check_cached.py",
                       "/tmp/minimum_dfs_check_cached.py",
                       "minimum_boundary_check.py", "other.py", ""):
            with self.assertRaisesRegex(ValueError, "DFS checker"):
                S.validate_dfs_checker(script)

    def test_regularity_command_uses_requested_generator_and_scope(self):
        argv = S.regularity_command(
            self.tile, "probe", 64, 8, 16, 12345, "census",
            "minimum_dfs_cached.py")
        self.assertEqual(os.path.basename(argv[1]),
                         "minimum_dfs_cached.py")
        self.assertEqual(argv[-2:], ["relative", "census"])
        self.assertEqual(argv[9:12], ["64", "8", "16"])
        ordinary = S.regularity_command(
            self.tile, "probe", 64, 0, 1, 100)
        self.assertEqual(os.path.basename(ordinary[1]), "minimum_dfs.py")

    def test_regularity_command_rejects_path_even_programmatically(self):
        with self.assertRaisesRegex(ValueError, "DFS script"):
            S.regularity_command(
                self.tile, "probe", 64, 0, 1, 100, "selected",
                "../minimum_dfs_cached.py")

    def generator_cmdline(self, script):
        argv = [sys.executable, script, "0.7", "0.8", "-0.45", "-0.35",
                "0.001", "3.13", "probe", "64", "0", "1", "100",
                "1e-6", "relative", "selected"]
        return b"\0".join(arg.encode() for arg in argv) + b"\0"

    @mock.patch.object(S, "process_is_alive", return_value=True)
    @mock.patch.object(S.os, "readlink", return_value=S.HERE)
    def test_adoption_recognizes_only_local_cached_generator(
            self, readlink, alive):
        by_key = {S.tile_key(self.tile): self.tile}
        local = os.path.join(S.HERE, "minimum_dfs_cached.py")
        with mock.patch("builtins.open", mock.mock_open(
                read_data=self.generator_cmdline(local))):
            self.assertIs(S.tile_from_generator_process(123, by_key),
                          self.tile)
        outside = os.path.join("/tmp", "minimum_dfs_cached.py")
        with mock.patch("builtins.open", mock.mock_open(
                read_data=self.generator_cmdline(outside))):
            self.assertIsNone(S.tile_from_generator_process(123, by_key))

    @mock.patch.object(S, "shard_matches", return_value=True)
    @mock.patch.object(S, "artifact_from_output",
                       side_effect=["shard0.json", "shard1.json"])
    @mock.patch.object(
        S, "run_checked",
        side_effect=["first", "second", "2 shards -> merged.json"])
    @mock.patch.object(S, "missing_ranges", return_value=[(0, 2)])
    def test_every_shard_uses_cached_generator_but_merge_does_not(
            self, missing, run_checked, artifact, matches):
        merged = S.generate_regularity_sharded(
            self.tile, shard_size=1, shard_budget=100,
            scope="census", dfs_script="minimum_dfs_cached.py")
        calls = [item.args[0] for item in run_checked.call_args_list]
        self.assertEqual(
            [os.path.basename(call[1]) for call in calls],
            ["minimum_dfs_cached.py", "minimum_dfs_cached.py",
             "minimum_dfs_merge.py"])
        self.assertTrue(all(call[-2:] == ["relative", "census"]
                            for call in calls[:2]))
        self.assertEqual(os.path.basename(merged), "merged.json")

    @mock.patch.object(S, "artifact_from_output", return_value="face.json")
    @mock.patch.object(S, "run_checked", return_value="wrote face.json")
    def test_boundary_generator_is_not_replaced(self, run_checked, artifact):
        S.generate_boundary(self.tile, "seam_lo", "census")
        argv = run_checked.call_args.args[0]
        self.assertEqual(os.path.basename(argv[1]),
                         "minimum_boundary_dfs.py")

    @mock.patch.object(S, "artifact_from_output", return_value="forest.json")
    @mock.patch.object(S, "run_checked", return_value="wrote forest.json")
    def test_unsharded_fixed_census_uses_cached_generator(
            self, run_checked, artifact):
        S.generate_regularity(
            self.tile, "census", "minimum_dfs_cached.py")
        argv = run_checked.call_args.args[0]
        self.assertEqual(os.path.basename(argv[1]),
                         "minimum_dfs_cached.py")
        self.assertEqual(argv[-2:], ["relative", "census"])

    @mock.patch.object(S, "run_checked")
    def test_replay_checker_defaults_to_ordinary_and_boundary_is_fixed(
            self, run_checked):
        S.replay_regularity("regularity.json")
        S.replay_boundary("boundary.json")
        scripts = [os.path.basename(item.args[0][1])
                   for item in run_checked.call_args_list]
        self.assertEqual(scripts,
                         ["minimum_dfs_check.py",
                          "minimum_boundary_check.py"])

    @mock.patch.object(S, "run_checked")
    def test_replay_checker_uses_requested_cached_entry_point(
            self, run_checked):
        S.replay_regularity(
            "regularity.json", "minimum_dfs_check_cached.py")
        argv = run_checked.call_args.args[0]
        self.assertEqual(os.path.basename(argv[1]),
                         "minimum_dfs_check_cached.py")
        self.assertEqual(argv[2], "regularity.json")

    @mock.patch.object(S, "run_checked")
    def test_replay_checker_uses_requested_parallel_entry_point(
            self, run_checked):
        S.replay_regularity(
            "regularity.json", "minimum_dfs_check_parallel.py")
        argv = run_checked.call_args.args[0]
        self.assertEqual(os.path.basename(argv[1]),
                         "minimum_dfs_check_parallel.py")
        self.assertEqual(argv[2], "regularity.json")

    def test_replay_checker_rejects_path_even_programmatically(self):
        with self.assertRaisesRegex(ValueError, "DFS checker"):
            S.replay_regularity(
                "regularity.json", "../minimum_dfs_check_cached.py")


class TileSelectionTests(unittest.TestCase):
    def test_census_candidate_rectangle_is_exactly_68_tiles(self):
        chosen = S.selected_tiles(y_min="-0.45", beta_min="0.5")
        self.assertEqual(len(chosen), 68)
        self.assertEqual(min(tile["b0"] for tile in chosen), Decimal("0.5"))
        self.assertEqual(max(tile["b1"] for tile in chosen), Decimal("2.2"))
        self.assertEqual(min(tile["y0"] for tile in chosen), Decimal("-0.45"))
        self.assertEqual(max(tile["y1"] for tile in chosen), Decimal("-0.1"))
        self.assertEqual(
            len({(tile["b0"], tile["b1"]) for tile in chosen}), 17)
        self.assertEqual(
            len({(tile["y0"], tile["y1"]) for tile in chosen}), 4)

    def test_cutoffs_do_not_change_the_default_189_tile_target(self):
        self.assertEqual(S.selected_tiles(), S.tiles())

    def test_census_manifest_records_exact_scope_count_and_target(self):
        chosen = S.selected_tiles(y_min="-0.45", beta_min="0.5")
        state = {"passed": {}}
        for tile in chosen:
            state["passed"][S.state_key(tile)] = {
                "regularity": "regularity.json",
                "boundaries": {face: face + ".json" for face in S.FACES},
            }
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "atlas.json")
            with open(os.devnull, "w") as sink, redirect_stdout(sink):
                S.write_atlas(state, path, scope="census", y_min="-0.45",
                              beta_min="0.5")
            with open(path) as handle:
                manifest = json.load(handle)
        self.assertEqual(manifest["scope"], "census")
        self.assertEqual(manifest["zero_count"], 2)
        self.assertNotIn("selected_count", manifest)
        self.assertEqual(manifest["witness"],
                         "minimum_witness_f4probe_census.json")
        self.assertEqual(manifest["base_box"], "f4probe_census")
        self.assertEqual(manifest["target"], {
            "b0": 0.5, "b1": 2.2, "y0": -0.45, "y1": -0.1})
        self.assertEqual(len(manifest["boxes"]), 68)


class ScopeValidationTests(unittest.TestCase):
    strip = {"seam": 0.137, "delta": 0.001, "dmax": 3.13}
    box = {
        "id": "box", "b0": 0.5, "b1": 0.6, "y0": -0.45, "y1": -0.35,
        "regularity": "regularity.json",
        "boundaries": {face: face + ".json" for face in C.FACES},
    }

    def artifact(self, coord, scope=None, face=None):
        spec = dict(self.strip, b0=0.5, b1=0.6, y0=-0.45, y1=-0.35,
                    coord=coord, n_parts=64)
        if scope is not None:
            spec["scope"] = scope
        if face is not None:
            spec["face"] = face
        return {"spec": spec, "lo_idx": 0, "hi_idx": 64}

    def test_selected_forest_cannot_satisfy_census_manifest(self):
        selected = self.artifact("relative")
        old_load = C.load
        try:
            C.load = lambda _path: selected
            with self.assertRaisesRegex(ValueError, "regularity spec mismatch"):
                C.validate_box(self.box, self.strip, "census")
        finally:
            C.load = old_load

    def test_census_forest_and_four_census_boundaries_are_accepted(self):
        docs = {"regularity.json": self.artifact("relative", "census")}
        docs.update({face + ".json": self.artifact(
            "absolute", "census", face) for face in C.FACES})
        old_load = C.load
        try:
            C.load = lambda path: docs[path]
            C.validate_box(self.box, self.strip, "census")
        finally:
            C.load = old_load


if __name__ == "__main__":
    unittest.main()
