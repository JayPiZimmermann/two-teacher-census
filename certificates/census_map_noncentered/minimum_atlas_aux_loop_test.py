import os
import tempfile
import unittest
from decimal import Decimal
from unittest import mock

import minimum_atlas_aux_loop as L
import minimum_atlas_aux_root as A
import minimum_atlas_sweep as S


def tile(b0):
    b0 = Decimal(b0)
    return {"b0": b0, "b1": b0 + Decimal("0.1"),
            "y0": Decimal("-0.25"), "y1": Decimal("-0.15")}


class AuxLoopTest(unittest.TestCase):
    def test_stop_sentinel_is_persistent_and_explicit(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "stop")
            self.assertFalse(L.stop_requested(None))
            self.assertFalse(L.stop_requested(path))
            with open(path, "w"):
                pass
            self.assertTrue(L.stop_requested(path))

    def test_active_frontier_requires_exact_worker_accounting(self):
        queue = [tile("0.5"), tile("0.6"), tile("0.7")]
        children = [{"tile_key": S.tile_key(queue[1])}]
        with mock.patch.object(A, "verify_scheduler"), \
                mock.patch.object(A, "direct_live_children",
                                  return_value=children):
            active, frontier = L.active_frontier(12, 1, "census", queue)
            self.assertEqual(active, {S.tile_key(queue[1])})
            self.assertEqual(frontier, 1)
            with self.assertRaises(A.AccountingError):
                L.active_frontier(12, 2, "census", queue)

    def test_choose_jobs_uses_two_distinct_farthest_first_missing_roots(self):
        queue = [tile(str(Decimal("0.5") + Decimal(i) / 10))
                 for i in range(5)]
        gaps = {
            S.tile_key(queue[4]): [(7, 64)],
            S.tile_key(queue[3]): [],
            S.tile_key(queue[2]): [(4, 64)],
        }
        with mock.patch.object(S, "missing_ranges",
                               side_effect=lambda item, scope:
                               gaps.get(S.tile_key(item), [])), \
                mock.patch.object(S, "expected_shard_path",
                                  return_value="/definitely/absent"):
            jobs = L.choose_jobs(queue, set(), 0, 1, "census")
        self.assertEqual([(index, job.root) for index, job in jobs],
                         [(4, 7), (2, 4)])
        self.assertNotEqual(jobs[0][1].tile_key, jobs[1][1].tile_key)

    def test_choose_jobs_refuses_occupied_first_missing_path(self):
        queue = [tile("0.5"), tile("0.6"), tile("0.7")]
        with tempfile.NamedTemporaryFile() as handle, \
                mock.patch.object(S, "missing_ranges", return_value=[(3, 64)]), \
                mock.patch.object(S, "expected_shard_path",
                                  return_value=handle.name):
            with self.assertRaises(A.NoClobberError):
                L.choose_jobs(queue, set(), 0, 1, "census")

    def test_queue_gap_stops_before_approaching_targets(self):
        queue = [tile(str(Decimal("0.5") + Decimal(i) / 10))
                 for i in range(5)]
        with mock.patch.object(S, "missing_ranges", return_value=[(0, 64)]), \
                mock.patch.object(S, "expected_shard_path",
                                  return_value="/definitely/absent"):
            self.assertEqual(L.choose_jobs(queue, set(), 2, 2, "census"), [])


if __name__ == "__main__":
    unittest.main()
