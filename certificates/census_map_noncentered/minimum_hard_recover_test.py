"""Structural tests for hard-cell recovery work selection."""
from types import SimpleNamespace
import unittest
from unittest import mock

import minimum_hard_recover as H
import minimum_shard_escalate as E


def tile(value):
    return {"b0": value, "b1": value + 0.1,
            "y0": -0.9, "y1": -0.85}


class HardRecoverTests(unittest.TestCase):
    def test_recorded_budget_cells_are_not_resubmitted_as_normal(self):
        first, second = tile(0.1), tile(0.2)
        missing = {E.tile_tag(first, 64): [3, 4],
                   E.tile_tag(second, 64): [9]}

        def missing_fn(item, n_parts):
            self.assertEqual(n_parts, 64)
            return missing[E.tile_tag(item, 64)]

        failed = {E.cell_key(first, 4, 64), E.cell_key(second, 9, 64)}
        normal, frontier = H.partition_work(
            [first, second], failed, missing_fn=missing_fn)
        self.assertEqual(normal, [(first, 3)])
        self.assertEqual(frontier, [(first, 4), (second, 9)])

    def test_fresh_state_keeps_recovery_names_separate(self):
        state = H.fresh_state()
        self.assertEqual(state["format"], H.FORMAT)
        self.assertEqual(state["n_parts"], 64)
        self.assertEqual(state["deferred_frontier"], {})
        self.assertNotEqual(H.STATE_PATH, E.state_path(64))
        self.assertNotEqual(H.LOG_DIR, E.log_dir(64))

    def test_legacy_unpadded_budget_keys_are_normalized(self):
        state = H.fresh_state()
        state["over_budget"] = {"example#7": 10}
        original = H.load_json
        try:
            H.load_json = lambda path: {"over_budget": {"other#61": 20}}
            self.assertEqual(H.recorded_over_budget(state),
                             {"example#007", "other#061"})
        finally:
            H.load_json = original

    def test_seeded_roots_skip_normal_and_enter_frontier(self):
        first, second = tile(0.1), tile(0.2)
        missing = {E.tile_tag(first, 64): [3, 4],
                   E.tile_tag(second, 64): [9]}

        def missing_fn(item, n_parts):
            self.assertEqual(n_parts, 64)
            return missing[E.tile_tag(item, 64)]

        seeds = H.validate_over_budget_seeds(
            [E.cell_key(first, 4, 64), E.cell_key(second, 9, 64)],
            [first, second], missing_fn=missing_fn)
        normal, frontier = H.partition_work(
            [first, second], seeds, missing_fn=missing_fn)
        self.assertEqual(normal, [(first, 3)])
        self.assertEqual(frontier, [(first, 4), (second, 9)])

    def test_seed_rejects_malformed_coordinates(self):
        target = tile(0.1)
        tag = E.tile_tag(target, 64)
        malformed = [tag, "#3", tag + "#x", tag + "#-1",
                     tag + "#64", tag + "#3#4"]
        for value in malformed:
            with self.subTest(value=value), self.assertRaises(ValueError):
                H.validate_over_budget_seeds(
                    [value], [target], missing_fn=lambda item, parts: [3])

    def test_seed_rejects_coordinate_outside_recovery_scope(self):
        target, outside = tile(0.1), tile(0.2)
        with self.assertRaisesRegex(ValueError, "outside"):
            H.validate_over_budget_seeds(
                [E.cell_key(outside, 3, 64)], [target],
                missing_fn=lambda item, parts: [3])

    def test_seed_rejects_already_passed_coordinate(self):
        target = tile(0.1)
        with self.assertRaisesRegex(ValueError, "usable shard"):
            H.validate_over_budget_seeds(
                [E.cell_key(target, 3, 64)], [target],
                missing_fn=lambda item, parts: [4])

    def test_exact_completed_seed_is_idempotent_only_when_enabled(self):
        target = tile(0.1)
        key = E.cell_key(target, 3, 64)
        seeded, deferred = H.validate_frontier_directives(
            [key], [], [target], missing_fn=lambda item, parts: [4],
            allow_exact_completed=True,
            completed_fn=lambda item, cell: item is target and cell == 3)
        self.assertEqual(seeded, set())
        self.assertEqual(deferred, set())

    def test_noncanonical_completed_seed_still_rejected(self):
        target = tile(0.1)
        with self.assertRaisesRegex(ValueError, "usable shard"):
            H.validate_frontier_directives(
                [E.cell_key(target, 3, 64)], [], [target],
                missing_fn=lambda item, parts: [4],
                allow_exact_completed=True,
                completed_fn=lambda item, cell: False)

    def test_deferred_root_has_distinct_provenance_and_enters_frontier(self):
        target = tile(0.1)
        missing_fn = lambda item, parts: [3, 4]
        seeded, deferred = H.validate_frontier_directives(
            [], [E.cell_key(target, 4, 64)], [target],
            missing_fn=missing_fn)
        self.assertEqual(seeded, set())
        self.assertEqual(deferred, {E.cell_key(target, 4, 64)})
        normal, frontier = H.partition_work(
            [target], seeded | deferred, missing_fn=missing_fn)
        self.assertEqual(normal, [(target, 3)])
        self.assertEqual(frontier, [(target, 4)])

    def test_root_cannot_be_seeded_and_deferred(self):
        target = tile(0.1)
        key = E.cell_key(target, 4, 64)
        with self.assertRaisesRegex(ValueError, "both"):
            H.validate_frontier_directives(
                [key], [key], [target],
                missing_fn=lambda item, parts: [4])

    def test_old_state_loads_with_no_deferred_roots(self):
        old = H.fresh_state()
        del old["deferred_frontier"]
        original = H.load_json
        try:
            H.load_json = lambda path: old
            loaded = H.load_state()
            self.assertEqual(loaded["deferred_frontier"], {})
        finally:
            H.load_json = original

    @mock.patch.object(H.E, "missing_cells", return_value=[])
    @mock.patch.object(H, "save_log")
    @mock.patch.object(H.subprocess, "run")
    def test_frontier_uses_requested_checker(
            self, run, save_log, missing):
        run.return_value = SimpleNamespace(returncode=0, stdout="green")
        args = SimpleNamespace(
            depth=8, max_depth=20, frontier_workers=4,
            frontier_budget=600000, full_width=0.001,
            dfs_checker="minimum_dfs_check_parallel.py")
        outcome = H.run_frontier(tile(0.1), 3, args)
        self.assertEqual(outcome["kind"], "passed")
        argv = run.call_args.args[0]
        self.assertEqual(argv[argv.index("--checker") + 1],
                         "minimum_dfs_check_parallel.py")


if __name__ == "__main__":
    unittest.main()
