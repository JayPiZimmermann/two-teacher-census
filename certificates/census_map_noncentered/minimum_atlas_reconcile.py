"""Adopt replayed escalation tiles into the selected adaptive atlas journal.

``minimum_shard_escalate.py`` and ``minimum_hard_recover.py`` deliberately
never write the adaptive journal, so they are safe beside a live refinement
run.  Once that run has stopped, this tool validates the stored artifact
specifications and replaces matching hard atomic leaves by passed leaves.  The
final atlas checker still independently replays every adopted forest and
boundary certificate.

Dry-run is the default.  ``--apply`` refuses to run while an adaptive driver is
alive in this directory.  ``--escalation-state`` selects either the default
64-cell journal or one refined-partition journal; each passed entry still
names a single full forest, so an atlas box never mixes shard specifications.
"""
import argparse
from decimal import Decimal
import glob
import json
import os

import minimum_atlas_adaptive as A
import minimum_atlas_check as C
import minimum_atlas_sweep as F
import minimum_hard_recover as H
import minimum_shard_escalate as E


HERE = os.path.dirname(os.path.abspath(__file__))
ADAPTIVE_DRIVER = os.path.join(HERE, "minimum_atlas_adaptive.py")


def is_this_adaptive_driver(cwd, argv):
    """Recognize this directory's driver, wherever it was launched from."""
    for arg in argv:
        if os.path.basename(arg) != "minimum_atlas_adaptive.py":
            continue
        path = arg if os.path.isabs(arg) else os.path.join(cwd, arg)
        if os.path.realpath(path) == os.path.realpath(ADAPTIVE_DRIVER):
            return True
    return False


def live_adaptive_drivers():
    found = []
    for proc in glob.glob("/proc/[0-9]*"):
        pid = int(os.path.basename(proc))
        try:
            cwd = os.path.realpath(os.readlink("/proc/%d/cwd" % pid))
            with open("/proc/%d/cmdline" % pid, "rb") as handle:
                argv = [part.decode() for part in handle.read().split(b"\0")
                        if part]
        except (OSError, UnicodeDecodeError):
            continue
        if not F.process_is_alive(pid):
            continue
        if is_this_adaptive_driver(cwd, argv):
            found.append(pid)
    return found


def index_of(values, value, name):
    target = Decimal(str(value))
    try:
        return values.index(target)
    except ValueError as exc:
        raise ValueError("escalation tile %s=%s is off the adaptive grid"
                         % (name, value)) from exc


def rect_of(item):
    rect = (index_of(A.BETA, item["b0"], "b0"),
            index_of(A.BETA, item["b1"], "b1"),
            index_of(A.Y, item["y0"], "y0"),
            index_of(A.Y, item["y1"], "y1"))
    if A.atomic_count(rect) != 1:
        raise ValueError("escalation result is not one atomic tile: %s"
                         % (rect,))
    return rect


def atlas_box(rect, item):
    tile = A.tile_of(rect)
    return {"id": item["id"],
            "b0": float(tile["b0"]), "b1": float(tile["b1"]),
            "y0": float(tile["y0"]), "y1": float(tile["y1"]),
            "regularity": item["regularity"],
            "boundaries": item["boundaries"]}


def validate_recovery_state(state):
    """Validate the narrow recovery-journal surface used by reconciliation."""
    if state.get("format") != H.FORMAT or state.get("strip") != F.STRIP \
            or state.get("n_parts") != 64:
        raise ValueError("hard-recovery journal has incompatible metadata")
    tiles = state.get("tiles")
    if not isinstance(tiles, dict):
        raise ValueError("hard-recovery journal has no tile mapping")
    for key, item in tiles.items():
        if not isinstance(key, str) or not isinstance(item, dict):
            raise ValueError("hard-recovery tile entry is not a mapping")
        kind = item.get("kind")
        if kind not in ("passed", "error"):
            raise ValueError("hard-recovery tile has unknown status")
        if kind != "passed":
            continue
        required = {"id", "b0", "b1", "y0", "y1", "regularity",
                    "boundaries"}
        if not required <= set(item) \
                or not isinstance(item["id"], str) \
                or not isinstance(item["regularity"], str) \
                or not isinstance(item["boundaries"], dict) \
                or set(item["boundaries"]) != set(F.FACES) \
                or not all(isinstance(path, str)
                           for path in item["boundaries"].values()):
            raise ValueError("hard-recovery passed tile has invalid schema")
        try:
            expected_key = ",".join(
                str(Decimal(str(item[name])))
                for name in ("b0", "b1", "y0", "y1"))
        except Exception as exc:
            raise ValueError("hard-recovery tile has invalid coordinates") \
                from exc
        if key != expected_key:
            raise ValueError("hard-recovery tile key disagrees with its box")
    return state


def load_reconciliation_state(path):
    """Load either an escalation journal or the separate recovery journal."""
    with open(path) as handle:
        state = json.load(handle)
    if not isinstance(state, dict):
        raise ValueError("reconciliation journal is not a mapping")
    if state.get("format") == H.FORMAT:
        return validate_recovery_state(state)
    # Preserve the escalation loader as the single validator for both n=64
    # and refined-partition escalation formats.
    return E.load_state_path(path)


def candidates(adaptive, escalation):
    out = []
    for _, item in sorted(escalation.get("tiles", {}).items()):
        if item.get("kind") != "passed":
            continue
        rect = rect_of(item)
        key = A.rect_key(rect)
        C.validate_box(atlas_box(rect, item), F.STRIP, "selected")
        if key in adaptive["passed"]:
            status = "already-passed"
        elif key in adaptive["hard_leaves"]:
            status = "adopt"
        else:
            status = "not-an-active-hard-leaf"
        out.append((key, rect, item, status))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument(
        "--escalation-state", default=E.STATE_PATH,
        help="64/refined escalation journal or minimum_hard_recover_state.json")
    args = parser.parse_args()
    # Check before reading the journal.  Checking only immediately before the
    # write admits a destructive race: a live driver can write its final
    # result and exit after this process loaded an older snapshot, after which
    # the stale snapshot would pass the liveness check and overwrite it.
    if args.apply:
        drivers = live_adaptive_drivers()
        if drivers:
            print("RECONCILIATION REFUSED: adaptive driver still alive: %s"
                  % ", ".join(map(str, drivers)))
            return 1
    try:
        adaptive = A.load_state()
        escalation_path = args.escalation_state
        if not os.path.isabs(escalation_path):
            escalation_path = os.path.join(HERE, escalation_path)
        escalation = load_reconciliation_state(escalation_path)
        items = candidates(adaptive, escalation)
    except (KeyError, ValueError, OSError, RuntimeError) as exc:
        print("RECONCILIATION INVALID:", exc)
        return 1
    for key, _, item, status in items:
        print("%-24s %-12s %s" % (key, status, item["regularity"]))
    adopt = [item for item in items if item[3] == "adopt"]
    conflicts = [item for item in items
                 if item[3] == "not-an-active-hard-leaf"]
    if conflicts:
        print("RECONCILIATION REFUSED: %d passed escalation tiles are not "
              "active hard leaves" % len(conflicts))
        return 1
    if not args.apply:
        print("dry run: %d tiles ready to adopt" % len(adopt))
        return 0
    # Repeat the check in case a driver was started while the artifacts were
    # being validated.  Cooperative operators should still start only one
    # journal writer at a time; these two checks close the exit-during-read
    # race that matters for the normal takeover sequence.
    drivers = live_adaptive_drivers()
    if drivers:
        print("RECONCILIATION REFUSED: adaptive driver still alive: %s"
              % ", ".join(map(str, drivers)))
        return 1
    for key, rect, item, _ in adopt:
        adaptive["passed"][key] = A.pass_record(
            rect, item, "escalation-journal")
        adaptive["hard_leaves"].pop(key, None)
    A.write_state(adaptive)
    print("adopted %d tiles; coverage is now %d/189" %
          (len(adopt), A.coverage(adaptive)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
