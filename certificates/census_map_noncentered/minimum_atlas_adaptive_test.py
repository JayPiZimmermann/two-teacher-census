"""Focused command-plumbing tests for adaptive regularity generation."""
import os
import unittest
from unittest import mock

import minimum_atlas_adaptive as A


class AdaptiveDfsScriptTests(unittest.TestCase):
    rect = (0, 1, 0, 1)

    @mock.patch.object(A, "save_log")
    @mock.patch.object(
        A, "run_process",
        return_value=(2, "budget exhausted at cell 0", 0.1, False))
    def test_unsharded_attempt_uses_requested_generator(
            self, run_process, save_log):
        outcome = A.certify(
            self.rect, 100, None, {}, {},
            dfs_script="minimum_dfs_cached.py")
        argv = run_process.call_args.args[0]
        self.assertEqual(os.path.basename(argv[1]),
                         "minimum_dfs_cached.py")
        self.assertEqual(argv[-2:], ["relative", "selected"])
        self.assertEqual(outcome["kind"], "budget")

    @mock.patch.object(A, "save_log")
    @mock.patch.object(
        A, "run_process",
        return_value=(2, "budget exhausted at cell 0", 0.1, False))
    def test_unsharded_default_remains_ordinary_generator(
            self, run_process, save_log):
        A.certify(self.rect, 100, None, {}, {})
        argv = run_process.call_args.args[0]
        self.assertEqual(os.path.basename(argv[1]), "minimum_dfs.py")

    @mock.patch.object(
        A.fixed, "generate_regularity_sharded",
        side_effect=RuntimeError("stop after command plumbing"))
    def test_sharded_attempt_forwards_requested_generator(self, generate):
        outcome = A.certify(
            self.rect, 100, None, {}, {}, sharded=True,
            shard_size=1, shard_budget=200,
            dfs_script="minimum_dfs_cached.py")
        self.assertEqual(outcome["kind"], "error")
        self.assertEqual(
            generate.call_args.kwargs["dfs_script"],
            "minimum_dfs_cached.py")


class AdaptiveDfsCheckerTests(unittest.TestCase):
    rect = (0, 1, 0, 1)

    def boundary_artifacts(self):
        key = A.fixed.tile_key(A.tile_of(self.rect))
        return {(key, face): "/tmp/%s.json" % face
                for face in A.fixed.FACES}

    @mock.patch.object(A.fixed, "replay_boundary")
    @mock.patch.object(A.fixed, "replay_regularity")
    def test_reused_artifact_uses_requested_checker(
            self, replay_regularity, replay_boundary):
        key = A.fixed.tile_key(A.tile_of(self.rect))
        outcome = A.certify(
            self.rect, 100, None, {key: "/tmp/existing.json"},
            self.boundary_artifacts(),
            dfs_checker="minimum_dfs_check_cached.py")
        self.assertEqual(outcome["kind"], "passed")
        self.assertTrue(outcome["reused"])
        replay_regularity.assert_called_once_with(
            "/tmp/existing.json", "minimum_dfs_check_cached.py")
        self.assertEqual(replay_boundary.call_count, len(A.fixed.FACES))

    @mock.patch.object(A.fixed, "replay_boundary")
    @mock.patch.object(A.fixed, "replay_regularity")
    def test_reused_artifact_default_remains_ordinary_checker(
            self, replay_regularity, replay_boundary):
        key = A.fixed.tile_key(A.tile_of(self.rect))
        outcome = A.certify(
            self.rect, 100, None, {key: "/tmp/existing.json"},
            self.boundary_artifacts())
        self.assertEqual(outcome["kind"], "passed")
        replay_regularity.assert_called_once_with(
            "/tmp/existing.json", "minimum_dfs_check.py")

    @mock.patch.object(A, "load_json", return_value={"stats": {
        "visited_nodes": 123}})
    @mock.patch.object(A.fixed, "replay_boundary")
    @mock.patch.object(A.fixed, "replay_regularity")
    @mock.patch.object(A.fixed, "generate_regularity_sharded",
                       return_value="/tmp/sharded.json")
    def test_sharded_artifact_uses_requested_checker(
            self, generate, replay_regularity, replay_boundary, load_json):
        outcome = A.certify(
            self.rect, 100, None, {}, self.boundary_artifacts(),
            sharded=True, shard_size=1, shard_budget=200,
            dfs_checker="minimum_dfs_check_cached.py")
        self.assertEqual(outcome["kind"], "passed")
        self.assertTrue(outcome["sharded"])
        replay_regularity.assert_called_once_with(
            "/tmp/sharded.json", "minimum_dfs_check_cached.py")

    @mock.patch.object(A, "save_log")
    @mock.patch.object(
        A, "run_process",
        return_value=(0, "wrote /tmp/unsharded.json", 0.1, False))
    @mock.patch.object(A.fixed, "artifact_from_output",
                       return_value="/tmp/unsharded.json")
    @mock.patch.object(A.fixed, "replay_boundary")
    @mock.patch.object(A.fixed, "replay_regularity")
    def test_unsharded_artifact_uses_requested_checker(
            self, replay_regularity, replay_boundary, artifact,
            run_process, save_log):
        outcome = A.certify(
            self.rect, 100, None, {}, self.boundary_artifacts(),
            dfs_checker="minimum_dfs_check_cached.py")
        self.assertEqual(outcome["kind"], "passed")
        self.assertFalse(outcome["reused"])
        replay_regularity.assert_called_once_with(
            "/tmp/unsharded.json", "minimum_dfs_check_cached.py")


if __name__ == "__main__":
    unittest.main()
