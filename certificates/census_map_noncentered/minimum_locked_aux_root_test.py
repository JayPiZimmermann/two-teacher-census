import fcntl
import os
import tempfile
import unittest
from unittest import mock

import minimum_atlas_aux_root as A
import minimum_atlas_sweep as S
import minimum_frontier_shard as F
import minimum_locked_aux_root as L


class LockedAuxRootTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.here = self.temp.name
        self.job = A.make_job(
            ["0.6", "0.7", "-0.25", "-0.15", "55"], "census")
        self.patch_here = mock.patch.object(S, "HERE", self.here)
        self.patch_here.start()
        self.addCleanup(self.patch_here.stop)
        self.verify = mock.patch.object(A, "verify_scheduler")
        self.verify.start()
        self.addCleanup(self.verify.stop)
        children = [{"pid": 1000 + index, "script": "minimum_dfs.py",
                     "tile_key": (S.tiles()[index + 20]["b0"],
                                  S.tiles()[index + 20]["b1"],
                                  S.tiles()[index + 20]["y0"],
                                  S.tiles()[index + 20]["y1"])}
                    for index in range(8)]
        self.children = mock.patch.object(
            A, "direct_live_children", return_value=children)
        self.children.start()
        self.addCleanup(self.children.stop)
        self.writers = mock.patch.object(
            F, "live_ordinary_writers", return_value=[])
        self.writers.start()
        self.addCleanup(self.writers.stop)

    def run_success(self):
        events = []
        tag, canonical, _, lock, _ = L.target_paths(self.job)

        def assert_locked():
            handle = open(lock, "a+")
            try:
                with self.assertRaises(BlockingIOError):
                    fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            finally:
                handle.close()

        def stage(job, script, budget, nonce):
            events.append("stage")
            assert_locked()
            staged = os.path.join(self.here, "staged.json")
            with open(staged, "wb") as handle:
                handle.write(b"exact-stage")
            return A.Publication(job, staged, canonical)

        def replay(publication, checker):
            events.append("replay")
            assert_locked()
            return L.file_hash(publication.staged_path)

        def publish(publications, pid, workers, scope):
            events.append("publish")
            assert_locked()
            os.link(publications[0].staged_path, canonical)

        publication, digest = L.run_locked(
            self.job, 42, 8, "minimum_dfs_cached.py", 600000,
            "minimum_dfs_check_parallel.py", 24.0, load_one=lambda: 5.0,
            stage_fn=stage, replay_fn=replay, publish_fn=publish)
        self.assertEqual(events, ["stage", "replay", "publish"])
        self.assertEqual(publication.canonical_path, canonical)
        self.assertEqual(digest, L.file_hash(canonical))
        handle = open(lock, "a+")
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        finally:
            handle.close()

    def test_lock_is_held_through_stage_replay_and_publication(self):
        self.run_success()

    def test_existing_canonical_is_rejected_before_stage(self):
        _, canonical, _, _, _ = L.target_paths(self.job)
        with open(canonical, "w") as handle:
            handle.write("occupied")
        with self.assertRaises(A.NoClobberError):
            L.run_locked(
                self.job, 42, 8, "minimum_dfs_cached.py", 600000,
                "minimum_dfs_check_parallel.py", 24.0,
                load_one=lambda: 5.0,
                stage_fn=mock.Mock(side_effect=AssertionError("ran")))

    def test_existing_frontier_state_is_rejected(self):
        tag, _, _, _, _ = L.target_paths(self.job)
        os.mkdir(os.path.join(
            self.here, "minimum_frontier_%s_p055_d10" % tag))
        with self.assertRaises(A.PublicationError):
            L.run_locked(
                self.job, 42, 8, "minimum_dfs_cached.py", 600000,
                "minimum_dfs_check_parallel.py", 24.0,
                load_one=lambda: 5.0)

    def test_live_ordinary_writer_is_rejected(self):
        with mock.patch.object(F, "live_ordinary_writers",
                               return_value=[1234]):
            with self.assertRaises(A.ActiveTargetError):
                L.run_locked(
                    self.job, 42, 8, "minimum_dfs_cached.py", 600000,
                    "minimum_dfs_check_parallel.py", 24.0,
                    load_one=lambda: 5.0)

    def test_held_frontier_lock_is_rejected(self):
        _, _, _, lock, _ = L.target_paths(self.job)
        handle = open(lock, "a+")
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            with self.assertRaises(A.ActiveTargetError):
                L.run_locked(
                    self.job, 42, 8, "minimum_dfs_cached.py", 600000,
                    "minimum_dfs_check_parallel.py", 24.0,
                    load_one=lambda: 5.0)
        finally:
            handle.close()

    def test_load_ceiling_is_strict(self):
        with self.assertRaises(A.PublicationError):
            L.run_locked(
                self.job, 42, 8, "minimum_dfs_cached.py", 600000,
                "minimum_dfs_check_parallel.py", 24.0,
                load_one=lambda: 24.0)

    def test_active_scheduler_tile_is_rejected(self):
        children = [{"pid": 1, "script": "minimum_dfs.py",
                     "tile_key": self.job.tile_key}] * 8
        with mock.patch.object(A, "direct_live_children",
                               return_value=children):
            with self.assertRaises(A.ActiveTargetError):
                L.run_locked(
                    self.job, 42, 8, "minimum_dfs_cached.py", 600000,
                    "minimum_dfs_check_parallel.py", 24.0,
                    load_one=lambda: 5.0)


if __name__ == "__main__":
    unittest.main()
