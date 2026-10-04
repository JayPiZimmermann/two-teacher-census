"""DELTA TABLE for the 2026-08-08 `sep_res` Jacobian sign fix.

The fix (see `sep_jacobian_audit.py` for what went wrong and how it was
caught) changes the Krawczyk verdicts that back every CERTIFIED per-family
statement of the noncentered separated stratum.  It cannot change the float
LOCATOR, which never touches the Jacobian.  So the question this script
answers, item by item, is: which PUBLISHED numbers move?

It re-runs `census_cert.certified_census` -- with the fixed Jacobian -- at
every point behind a published artifact, and diffs against the stored file:

  faces_0.json        the 28 face witnesses (CERTIFICATE.md section 7,
                      website appendix-c per-face censuses)
  wall_pos.json       the mixed-collar pointwise witness scan (section 6.3)
  ywallp_0.json       the mixed-sector candidate-transition scan (section 6.3)
  bandfix_0.json      the 158 antipodal band points (section 8), which also
                      feed the shipped-classifier agreement count

Compared per point: the census STRING, `n_separated`, the Schur types, and
the `certified` flag.  Any difference is printed as `was X, now Y`.

Usage:  python3 jacobian_fix_delta.py [faces|wall|ywall|band|all]
Output: jacobian_fix_delta.json
"""
import json
import os
import signal
import sys
import time

import census_cert as X

# Some band points near beta = pi make `certified_census` grind (the torque
# root solver stalls where the lens degenerates).  A per-point alarm keeps the
# audit finite and records the point as TIMEOUT rather than silently dropping
# it -- a skipped point is exactly the kind of hole a delta table must not
# have.
PER_POINT_SECONDS = 180


class _Timeout(Exception):
    pass


def _alarm(signum, frame):                                       # noqa
    raise _Timeout()

HERE = os.path.dirname(os.path.abspath(__file__))


def recompute(beta, y, N=160):
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(PER_POINT_SECONDS)
    try:
        r = X.certified_census(beta, y, N=N)
    finally:
        signal.alarm(0)
    return {
        "census": r["census"],
        "n_separated": r["n_separated"],
        "sep_schur": [s["schur"] for s in r["separated"]],
        "sep_D": [round(s["D_mid"], 6) for s in r["separated"]],
        "certified": r["certified_modulo_locator_completeness"],
        "n_failed_cert": len(r["locator_failed_certification"]),
        "n_candidates": r["locator_candidates"],
    }


def diff_one(stored, got, keys):
    out = {}
    for k, sk in keys:
        if sk not in stored:
            continue
        a, b = stored[sk], got[k]
        if isinstance(a, list) and isinstance(b, list):
            same = (len(a) == len(b)
                    and all(_close(x, z) for x, z in zip(a, b)))
        else:
            same = (a == b)
        if not same:
            out[k] = {"was": a, "now": b}
    return out


def _close(a, b):
    if isinstance(a, float) and isinstance(b, float):
        return abs(a - b) < 1e-5
    return a == b


KEYS = [("census", "census"), ("n_separated", "n_separated"),
        ("sep_schur", "sep_schur"), ("certified", "certified")]
KEYS_YWALL = [("census", "census"), ("n_separated", "n_sep"),
              ("certified", "certified")]


def run_set(name, points, keys, N=160):
    rows = []
    t0 = time.time()
    for i, (label, beta, y, stored) in enumerate(points):
        try:
            got = recompute(beta, y, N=N)
            d = diff_one(stored, got, keys)
        except _Timeout:
            got, d = {"timeout": PER_POINT_SECONDS}, {"TIMEOUT": True}
        except Exception as e:                                   # noqa
            got, d = {"error": repr(e)}, {"error": repr(e)}
        rows.append({"set": name, "label": label, "beta": beta, "y": y,
                     "delta": d, "now": got})
        print("  %-22s beta=%-6s y=%-8s %s" % (
            label, beta, y, "UNCHANGED" if not d else "CHANGED %s" % d),
            flush=True)
    n_ch = sum(1 for r in rows if r["delta"])
    print("%s: %d points, %d CHANGED, %.0fs"
          % (name, len(rows), n_ch, time.time() - t0), flush=True)
    return rows, n_ch


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    slot = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    nslots = int(sys.argv[3]) if len(sys.argv) > 3 else 1

    def mine(seq):
        return [z for n, z in enumerate(seq) if n % nslots == slot]
    out = {}
    allrows = []

    if which in ("faces", "all"):
        f = json.load(open(os.path.join(HERE, "faces_0.json")))
        pts = [("%s%s" % (w["face"], "*" if w.get("is_rep") else ""),
                w["beta"], w["y"], w) for w in f]
        rows, n = run_set("faces_28", mine(pts), KEYS)
        out["faces_28"] = {"n_points": len(rows), "n_changed": n}
        allrows += rows

    # `collar_wall` and `fold_wall` below are retained only as stable legacy
    # set keys in the committed delta JSON. They do not certify a wall; the
    # artifact's `legacy_set_names` field records their current meanings.
    if which in ("wall", "all"):
        w = json.load(open(os.path.join(HERE, "wall_pos.json")))["points"]
        pts = [("collar", z["beta"], z["y"], z) for z in w]
        rows, n = run_set("collar_wall", mine(pts), KEYS)
        out["collar_wall"] = {"n_points": len(rows), "n_changed": n}
        allrows += rows

    if which in ("ywall", "all"):
        w = json.load(open(os.path.join(HERE, "ywallp_0.json")))
        pts = [("foldwall", z["beta"], z["y"], z) for z in w]
        rows, n = run_set("fold_wall", mine(pts), KEYS_YWALL)
        out["fold_wall"] = {"n_points": len(rows), "n_changed": n}
        allrows += rows

    if which in ("band", "all"):
        b = json.load(open(os.path.join(HERE, "bandfix_0.json")))
        pts = [("band", z["beta"], z["y"], z) for z in b]
        rows, n = run_set("band_158", mine(pts), KEYS)
        out["band_158"] = {"n_points": len(rows), "n_changed": n}
        allrows += rows

    out["rows"] = allrows
    out["total_changed"] = sum(1 for r in allrows if r["delta"])
    out["total_points"] = len(allrows)
    with open(os.path.join(HERE, "jacobian_fix_delta_%d.json" % slot),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nTOTAL: %d points, %d CHANGED"
          % (out["total_points"], out["total_changed"]))


if __name__ == "__main__":
    main()
