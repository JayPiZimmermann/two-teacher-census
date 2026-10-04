"""Race and artifact guards for staged auxiliary atlas roots."""
import json
import os
import signal
import tempfile
import unittest
from decimal import Decimal
from unittest import mock

import minimum_atlas_aux_root as A


class AuxiliaryRootTests(unittest.TestCase):
    tile_key = (Decimal("2.1"), Decimal("2.2"),
                Decimal("-0.25"), Decimal("-0.15"))

    def job(self):
        return A.RootJob(self.tile_key, 48, "census")

    def recognized_children(self, active=False):
        keys = []
        for index in range(8):
            keys.append((Decimal("0.5") + Decimal(index) / 10,
                         Decimal("0.6") + Decimal(index) / 10,
                         Decimal("-0.45"), Decimal("-0.35")))
        if active:
            keys[3] = self.tile_key
        return [{"pid": 1000 + index, "script": "minimum_dfs_cached.py",
                 "tile_key": key}
                for index, key in enumerate(keys)]

    def publication(self, directory):
        staged = os.path.join(directory, "staged.json")
        canonical = os.path.join(directory, "canonical.json")
        with open(staged, "w") as handle:
            handle.write("staged")
        return A.Publication(self.job(), staged, canonical)

    @mock.patch.object(A, "validate_staged")
    @mock.patch.object(A, "direct_live_children")
    @mock.patch.object(A, "wait_stopped")
    @mock.patch.object(A, "verify_scheduler")
    @mock.patch.object(A.S, "process_is_alive", return_value=True)
    @mock.patch.object(A.os, "kill")
    def test_active_target_refuses_publication_and_always_resumes(
            self, kill, alive, verify, wait, children, validate):
        children.return_value = self.recognized_children(active=True)
        with tempfile.TemporaryDirectory() as directory:
            publication = self.publication(directory)
            with self.assertRaises(A.ActiveTargetError):
                A.guarded_publish([publication], 77, 8, "census")
            self.assertFalse(os.path.exists(publication.canonical_path))
        self.assertEqual(kill.call_args_list, [
            mock.call(77, signal.SIGSTOP),
            mock.call(77, signal.SIGCONT),
        ])

    @mock.patch.object(A, "validate_staged")
    @mock.patch.object(A, "direct_live_children")
    @mock.patch.object(A, "wait_stopped")
    @mock.patch.object(A, "verify_scheduler")
    @mock.patch.object(A.S, "process_is_alive", return_value=True)
    @mock.patch.object(A.os, "kill")
    def test_existing_canonical_is_never_clobbered_and_parent_resumes(
            self, kill, alive, verify, wait, children, validate):
        children.return_value = self.recognized_children(active=False)
        with tempfile.TemporaryDirectory() as directory:
            publication = self.publication(directory)
            with open(publication.canonical_path, "w") as handle:
                handle.write("original")
            with self.assertRaises(A.NoClobberError):
                A.guarded_publish([publication], 77, 8, "census")
            with open(publication.canonical_path) as handle:
                self.assertEqual(handle.read(), "original")
        self.assertEqual(kill.call_args_list, [
            mock.call(77, signal.SIGSTOP),
            mock.call(77, signal.SIGCONT),
        ])

    @mock.patch.object(A, "validate_staged")
    @mock.patch.object(A, "direct_live_children")
    @mock.patch.object(A, "wait_stopped")
    @mock.patch.object(A, "verify_scheduler")
    @mock.patch.object(A.S, "process_is_alive", return_value=True)
    @mock.patch.object(A.os, "kill")
    def test_success_uses_hard_link_and_resumes_parent(
            self, kill, alive, verify, wait, children, validate):
        children.return_value = self.recognized_children(active=False)
        with tempfile.TemporaryDirectory() as directory:
            publication = self.publication(directory)
            A.guarded_publish([publication], 77, 8, "census")
            self.assertEqual(
                os.stat(publication.staged_path).st_ino,
                os.stat(publication.canonical_path).st_ino)
        self.assertEqual(kill.call_args_list, [
            mock.call(77, signal.SIGSTOP),
            mock.call(77, signal.SIGCONT),
        ])

    def write_artifact(self, path, scope="census", undecided=0):
        spec = {
            "b0": 2.1, "b1": 2.2, "y0": -0.25, "y1": -0.15,
            "seam": 0.137, "delta": 0.001, "dmax": 3.13,
            "n_parts": 64, "minw": 1e-6, "mv_maxw": 0.05,
            "row_mv_maxw": 0.1, "kraw_maxw": 0.1,
            "coord": "relative",
            "split": "relative-u-pi-endpoints-then-widest-midpoint",
            "scope": scope,
        }
        doc = {
            "spec": spec, "lo_idx": 48, "hi_idx": 49,
            "n_nodes": 4, "bits_b64": "AA==",
            "stats": {"visited_nodes": 1,
                      "counts": {"UNDECIDED": undecided}},
        }
        with open(path, "w") as handle:
            json.dump(doc, handle)

    def test_staged_validation_checks_scope_range_and_completeness(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "shard.json")
            self.write_artifact(path)
            A.validate_staged(self.job(), path)
            self.write_artifact(path, scope="selected")
            with self.assertRaisesRegex(A.PublicationError, "mismatch"):
                A.validate_staged(self.job(), path)
            self.write_artifact(path, undecided=1)
            with self.assertRaisesRegex(A.PublicationError, "not complete"):
                A.validate_staged(self.job(), path)


if __name__ == "__main__":
    unittest.main()
