import threading
import unittest
from decimal import Decimal
from unittest import mock

import minimum_atlas_aux_lanes as N
import minimum_atlas_aux_root as A
import minimum_atlas_sweep as S


def tile(b0):
    b0 = Decimal(b0)
    return {"b0": b0, "b1": b0 + Decimal("0.1"),
            "y0": Decimal("-0.25"), "y1": Decimal("-0.15")}


def publication(job):
    return A.Publication(job, "staged", "canonical")


class AuxiliaryLanesTest(unittest.TestCase):
    def setUp(self):
        self.queue = [tile("0.5"), tile("0.6"), tile("0.7")]
        self.jobs = [A.RootJob(S.tile_key(self.queue[2]), root, "census")
                     for root in (0, 1)]
        self.other = A.RootJob(S.tile_key(self.queue[1]), 0, "census")

    def test_asymmetric_lane_refills_before_slow_lane_finishes(self):
        fast_refill_started = threading.Event()
        slow_finished = threading.Event()
        selections = [
            (2, self.jobs[0]), (1, self.other), (2, self.jobs[1]), None,
        ]
        lock = threading.Lock()
        active = 0
        max_active = 0

        def choose(*args, **kwargs):
            return selections.pop(0) if selections else None

        def stage(indexed, script, budget):
            nonlocal active, max_active
            _, job = indexed
            with lock:
                active += 1
                max_active = max(max_active, active)
            try:
                if job == self.other:
                    self.assertTrue(fast_refill_started.wait(1.0))
                    slow_finished.set()
                elif job == self.jobs[1]:
                    fast_refill_started.set()
                return indexed[0], publication(job)
            finally:
                with lock:
                    active -= 1

        published = []
        with mock.patch.object(N.L, "stop_requested",
                               side_effect=lambda path:
                               slow_finished.is_set()), \
                mock.patch.object(N.L, "active_frontier",
                                  return_value=(set(), 0)), \
                mock.patch.object(N, "choose_one", side_effect=choose), \
                mock.patch.object(N, "stage_one", side_effect=stage), \
                mock.patch.object(A, "guarded_publish",
                                  side_effect=lambda items, *args:
                                  published.extend(items)):
            count, started, _ = N.run_lanes(
                self.queue, 10, 8, "census", "minimum_dfs_cached.py",
                300000, 1, 24, "stop", load_one=lambda: 0)
        self.assertTrue(fast_refill_started.is_set())
        self.assertTrue(slow_finished.is_set())
        self.assertEqual((count, started, len(published)), (3, 3, 3))
        self.assertLessEqual(max_active, 2)

    def test_high_load_after_completion_refuses_refill(self):
        choices = [(2, self.jobs[0]), (1, self.other)]
        loads = iter((10, 10, 24, 24))
        with mock.patch.object(N.L, "stop_requested", return_value=False), \
                mock.patch.object(N.L, "active_frontier",
                                  return_value=(set(), 0)), \
                mock.patch.object(N, "choose_one",
                                  side_effect=lambda *args, **kwargs:
                                  choices.pop(0)), \
                mock.patch.object(N, "stage_one",
                                  side_effect=lambda indexed, *args:
                                  (indexed[0], publication(indexed[1]))), \
                mock.patch.object(A, "guarded_publish"):
            count, started, reason = N.run_lanes(
                self.queue, 10, 8, "census", "minimum_dfs_cached.py",
                300000, 1, 24, None, load_one=lambda: next(loads))
        self.assertEqual((count, started), (2, 2))
        self.assertIn("load 24.00", reason)

    def test_stop_sentinel_drains_without_refill(self):
        choices = [(2, self.jobs[0]), (1, self.other)]
        stop_checks = iter((False, False, True, True))
        with mock.patch.object(N.L, "stop_requested",
                               side_effect=lambda path: next(stop_checks)), \
                mock.patch.object(N.L, "active_frontier",
                                  return_value=(set(), 0)), \
                mock.patch.object(N, "choose_one",
                                  side_effect=lambda *args, **kwargs:
                                  choices.pop(0)), \
                mock.patch.object(N, "stage_one",
                                  side_effect=lambda indexed, *args:
                                  (indexed[0], publication(indexed[1]))), \
                mock.patch.object(A, "guarded_publish"):
            count, started, reason = N.run_lanes(
                self.queue, 10, 8, "census", "minimum_dfs_cached.py",
                300000, 1, 24, "stop", load_one=lambda: 0)
        self.assertEqual((count, started), (2, 2))
        self.assertEqual(reason, "batch-boundary sentinel")

    def test_guard_refusal_stops_without_refill(self):
        choices = [(2, self.jobs[0]), (1, self.other)]
        with mock.patch.object(N.L, "stop_requested", return_value=False), \
                mock.patch.object(N.L, "active_frontier",
                                  return_value=(set(), 0)), \
                mock.patch.object(N, "choose_one",
                                  side_effect=lambda *args, **kwargs:
                                  choices.pop(0)), \
                mock.patch.object(N, "stage_one",
                                  side_effect=lambda indexed, *args:
                                  (indexed[0], publication(indexed[1]))), \
                mock.patch.object(A, "guarded_publish",
                                  side_effect=A.ActiveTargetError("active")):
            with self.assertRaises(A.ActiveTargetError):
                N.run_lanes(
                    self.queue, 10, 8, "census", "minimum_dfs_cached.py",
                    300000, 1, 24, None, load_one=lambda: 0)


if __name__ == "__main__":
    unittest.main()
