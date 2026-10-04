"""Structural tests for escalation/adaptive journal reconciliation."""
import copy
import json
import os
import sys
import tempfile
import unittest
from unittest import mock

import minimum_atlas_reconcile as R


def tile_result(b0=0.1, b1=0.2, y0=-0.9, y1=-0.85):
    return {"kind": "passed", "id": "tile", "b0": b0, "b1": b1,
            "y0": y0, "y1": y1, "regularity": "regularity.json",
            "boundaries": {face: face + ".json" for face in R.F.FACES}}


class ReconcileTests(unittest.TestCase):
    def test_driver_recognition_does_not_depend_on_launch_directory(self):
        parent = os.path.dirname(R.HERE)
        relative = os.path.join(os.path.basename(R.HERE),
                                "minimum_atlas_adaptive.py")
        self.assertTrue(R.is_this_adaptive_driver(
            R.HERE, ["python3", "-u", "minimum_atlas_adaptive.py"]))
        self.assertTrue(R.is_this_adaptive_driver(
            parent, ["python3", "-u", relative]))
        self.assertFalse(R.is_this_adaptive_driver(
            parent, ["python3", "minimum_atlas_adaptive.py"]))

    @mock.patch.object(R.A, "load_state")
    @mock.patch.object(R, "live_adaptive_drivers", return_value=[123])
    def test_apply_refuses_before_loading_a_live_journal(self, live, load):
        with mock.patch.object(sys, "argv", ["reconcile", "--apply"]), \
                mock.patch("builtins.print"):
            self.assertEqual(R.main(), 1)
        load.assert_not_called()

    def test_atomic_grid_tile_to_rect(self):
        self.assertEqual(R.rect_of(tile_result()), (0, 1, 0, 1))

    def test_coarse_and_off_grid_tiles_rejected(self):
        with self.assertRaisesRegex(ValueError, "not one atomic tile"):
            R.rect_of(tile_result(b1=0.3))
        with self.assertRaisesRegex(ValueError, "off the adaptive grid"):
            R.rect_of(tile_result(b0=0.11))

    @mock.patch.object(R.C, "validate_box")
    def test_only_active_hard_leaf_is_adoptable(self, validate):
        item = tile_result()
        escalation = {"tiles": {"fixed-key": item}}
        for state, expected in (
                ({"passed": {}, "hard_leaves": {"0:1,0:1": {}}}, "adopt"),
                ({"passed": {"0:1,0:1": {}}, "hard_leaves": {}},
                 "already-passed"),
                ({"passed": {}, "hard_leaves": {}},
                 "not-an-active-hard-leaf")):
            self.assertEqual(R.candidates(state, escalation)[0][3], expected)
        self.assertEqual(validate.call_count, 3)

    def test_recovery_journal_is_accepted_by_format_dispatch(self):
        item = tile_result()
        state = {"format": R.H.FORMAT, "strip": R.F.STRIP, "n_parts": 64,
                 "cells": {}, "over_budget": {}, "frontiers": {},
                 "tiles": {"0.1,0.2,-0.9,-0.85": item}}
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "recovery.json")
            with open(path, "w") as handle:
                json.dump(state, handle)
            self.assertEqual(R.load_reconciliation_state(path), state)

    def test_malformed_recovery_journals_are_rejected(self):
        item = tile_result()
        base = {"format": R.H.FORMAT, "strip": R.F.STRIP, "n_parts": 64,
                "tiles": {"0.1,0.2,-0.9,-0.85": item}}
        bad_states = []
        for mutate in (
                lambda state: state.update(n_parts=512),
                lambda state: state.update(strip={"delta": 0.001}),
                lambda state: state.update(tiles=[]),
                lambda state: state["tiles"].update(
                    {"0.1,0.2,-0.9,-0.85": {"kind": "unknown"}}),
                lambda state: state["tiles"][
                    "0.1,0.2,-0.9,-0.85"]["boundaries"].pop("gap_hi"),
                lambda state: state.update(
                    tiles={"wrong-key": tile_result()})):
            state = copy.deepcopy(base)
            mutate(state)
            bad_states.append(state)
        for state in bad_states:
            with self.assertRaises(ValueError):
                R.validate_recovery_state(state)


if __name__ == "__main__":
    unittest.main()
