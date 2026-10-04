import fcntl
import json
import os
import tempfile
import unittest
from unittest import mock

import minimum_atlas_aux_root as A
import minimum_atlas_sweep as S
import minimum_frontier_shard as F
import minimum_locked_aux_adopt as P
import minimum_locked_aux_root as L


class PostSweepAuxAdoptTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.here = self.temp.name
        self.job = A.make_job(
            ["0.6", "0.7", "-0.25", "-0.15", "55"], "census")
        self.patch_here = mock.patch.object(S, "HERE", self.here)
        self.patch_here.start()
        self.addCleanup(self.patch_here.stop)
        self.writers = mock.patch.object(
            F, "live_ordinary_writers", return_value=[])
        self.writers.start()
        self.addCleanup(self.writers.stop)
        tag = S.tile_id(self.job.tile, self.job.scope)
        self.stage = os.path.join(
            self.here, "minimum_dfs_auxstage_test_%s_r055_p055.json" % tag)
        self.write_stage(self.stage)

    def write_stage(self, path, undecided=0):
        spec = {
            "b0": 0.6, "b1": 0.7, "y0": -0.25, "y1": -0.15,
            "seam": 0.137, "delta": 0.001, "dmax": 3.13,
            "n_parts": 64, "minw": 1e-6, "mv_maxw": 0.05,
            "row_mv_maxw": 0.1, "kraw_maxw": 0.1,
            "coord": "relative",
            "split": "relative-u-pi-endpoints-then-widest-midpoint",
            "scope": "census",
        }
        doc = {
            "spec": spec, "lo_idx": 55, "hi_idx": 56,
            "n_nodes": 1, "bits_b64": "AA==",
            "stats": {"visited_nodes": 1,
                      "counts": {"UNDECIDED": undecided}},
        }
        with open(path, "w") as handle:
            json.dump(doc, handle)

    def replay(self, publication, checker):
        self.assertEqual(publication.staged_path, self.stage)
        self.assertEqual(checker, "minimum_dfs_check_parallel.py")
        return L.file_hash(self.stage)

    def run_success(self, sweep_fn=lambda scope: []):
        return P.adopt_stage(
            self.job, self.stage, 424242,
            "minimum_dfs_check_parallel.py", sweep_fn=sweep_fn,
            replay_fn=self.replay, pid_exists=lambda path: False)

    def test_success_replays_and_hard_links_exact_stage(self):
        canonical, digest = self.run_success()
        self.assertEqual(digest, L.file_hash(self.stage))
        self.assertEqual(os.stat(canonical).st_ino,
                         os.stat(self.stage).st_ino)

    def test_named_scheduler_must_be_absent(self):
        with self.assertRaisesRegex(A.AccountingError, "has not exited"):
            P.adopt_stage(
                self.job, self.stage, os.getpid(),
                "minimum_dfs_check_parallel.py", sweep_fn=lambda scope: [],
                replay_fn=self.replay)

    def test_any_new_local_sweep_is_rejected(self):
        with self.assertRaisesRegex(A.AccountingError, "sweep still alive"):
            P.adopt_stage(
                self.job, self.stage, 424242,
                "minimum_dfs_check_parallel.py", sweep_fn=lambda scope: [91],
                replay_fn=self.replay, pid_exists=lambda path: False)

    def test_sweep_reappearing_after_replay_is_rejected(self):
        sweeps = mock.Mock(side_effect=[[], [92]])
        with self.assertRaisesRegex(A.AccountingError, "sweep still alive"):
            P.adopt_stage(
                self.job, self.stage, 424242,
                "minimum_dfs_check_parallel.py", sweep_fn=sweeps,
                replay_fn=self.replay, pid_exists=lambda path: False)
        self.assertFalse(os.path.exists(L.target_paths(self.job)[1]))

    def test_named_pid_reappearing_after_replay_is_rejected(self):
        exists = mock.Mock(side_effect=[False, True])
        with self.assertRaisesRegex(A.AccountingError, "has not exited"):
            P.adopt_stage(
                self.job, self.stage, 424242,
                "minimum_dfs_check_parallel.py",
                sweep_fn=lambda scope: [], replay_fn=self.replay,
                pid_exists=exists)
        self.assertFalse(os.path.exists(L.target_paths(self.job)[1]))

    def test_existing_canonical_candidate_and_frontier_are_rejected(self):
        _, canonical, candidate, _, frontier_glob = L.target_paths(self.job)
        for path in (canonical, candidate, frontier_glob.replace("*", "10")):
            with self.subTest(path=path):
                if os.path.isdir(path):
                    os.rmdir(path)
                elif os.path.exists(path):
                    os.unlink(path)
                if path.endswith("10"):
                    os.mkdir(path)
                else:
                    with open(path, "w") as handle:
                        handle.write("occupied")
                with self.assertRaises((A.NoClobberError,
                                        A.PublicationError)):
                    P.adopt_stage(
                        self.job, self.stage, 424242,
                        "minimum_dfs_check_parallel.py",
                        sweep_fn=lambda scope: [], replay_fn=self.replay,
                        pid_exists=lambda path: False)
                if os.path.isdir(path):
                    os.rmdir(path)
                else:
                    os.unlink(path)

    def test_live_writer_and_held_lock_are_rejected(self):
        with mock.patch.object(F, "live_ordinary_writers",
                               return_value=[93]):
            with self.assertRaises(A.ActiveTargetError):
                P.adopt_stage(
                    self.job, self.stage, 424242,
                    "minimum_dfs_check_parallel.py",
                    sweep_fn=lambda scope: [], replay_fn=self.replay,
                    pid_exists=lambda path: False)
        lock = L.target_paths(self.job)[3]
        handle = open(lock, "a+")
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            with self.assertRaisesRegex(A.ActiveTargetError, "lock"):
                P.adopt_stage(
                    self.job, self.stage, 424242,
                    "minimum_dfs_check_parallel.py",
                    sweep_fn=lambda scope: [], replay_fn=self.replay,
                    pid_exists=lambda path: False)
        finally:
            handle.close()

    def test_wrong_name_incomplete_and_changed_stage_are_rejected(self):
        wrong = os.path.join(self.here, "wrong.json")
        self.write_stage(wrong)
        with self.assertRaisesRegex(A.PublicationError, "name"):
            P.adopt_stage(
                self.job, wrong, 424242,
                "minimum_dfs_check_parallel.py", sweep_fn=lambda scope: [],
                replay_fn=self.replay, pid_exists=lambda path: False)
        self.write_stage(self.stage, undecided=1)
        with self.assertRaisesRegex(A.PublicationError, "not complete"):
            P.adopt_stage(
                self.job, self.stage, 424242,
                "minimum_dfs_check_parallel.py", sweep_fn=lambda scope: [],
                replay_fn=self.replay, pid_exists=lambda path: False)
        self.write_stage(self.stage)

        def changing_replay(publication, checker):
            before = L.file_hash(self.stage)
            with open(self.stage, "a") as handle:
                handle.write("changed")
            return before

        with self.assertRaisesRegex(A.PublicationError, "changed"):
            P.adopt_stage(
                self.job, self.stage, 424242,
                "minimum_dfs_check_parallel.py", sweep_fn=lambda scope: [],
                replay_fn=changing_replay, pid_exists=lambda path: False)

    def test_link_race_never_clobbers_existing_canonical(self):
        canonical = L.target_paths(self.job)[1]

        def racing_link(source, target):
            self.assertEqual(target, canonical)
            with open(target, "w") as handle:
                handle.write("racer")
            raise FileExistsError(target)

        with self.assertRaises(A.NoClobberError):
            P.adopt_stage(
                self.job, self.stage, 424242,
                "minimum_dfs_check_parallel.py",
                sweep_fn=lambda scope: [], replay_fn=self.replay,
                link_fn=racing_link, pid_exists=lambda path: False)
        with open(canonical) as handle:
            self.assertEqual(handle.read(), "racer")

    def test_post_link_verification_failure_removes_only_our_link(self):
        canonical = L.target_paths(self.job)[1]

        def corrupting_link(source, target):
            os.link(source, target)
            with open(source, "a") as handle:
                handle.write("corrupt")

        with self.assertRaises(Exception):
            P.adopt_stage(
                self.job, self.stage, 424242,
                "minimum_dfs_check_parallel.py",
                sweep_fn=lambda scope: [], replay_fn=self.replay,
                link_fn=corrupting_link, pid_exists=lambda path: False)
        self.assertFalse(os.path.lexists(canonical))

    def test_post_link_failure_does_not_remove_replaced_path(self):
        canonical = L.target_paths(self.job)[1]

        def replacing_link(source, target):
            os.link(source, target)
            os.unlink(target)
            with open(target, "w") as handle:
                handle.write("replacement")

        with self.assertRaises(A.PublicationError):
            P.adopt_stage(
                self.job, self.stage, 424242,
                "minimum_dfs_check_parallel.py",
                sweep_fn=lambda scope: [], replay_fn=self.replay,
                link_fn=replacing_link, pid_exists=lambda path: False)
        with open(canonical) as handle:
            self.assertEqual(handle.read(), "replacement")


if __name__ == "__main__":
    unittest.main()
