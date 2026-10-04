"""Re-run the antipodal-band certified census AFTER the sep_res Jacobian sign
fix, at exactly the 158 points of the pre-fix scan (bandfix_0.json), in the
same record format, so that compare_band.py can be re-run against the shipped
classifier.  Identical to band_scan.py except that (a) the points are read
from the stored file rather than regenerated, and (b) a per-point alarm keeps
the run finite and records a TIMEOUT instead of silently dropping the point.

Usage: python3 band_rescan_fixed.py <slot> <nslots> [tag]  -> band<tag>_<slot>.json
"""
import json, os, signal, sys, time
import census_cert as X

HERE = os.path.dirname(os.path.abspath(__file__))
PER_POINT = 240


class _T(Exception):
    pass


def _a(s, f):                                                    # noqa
    raise _T()


slot, nslots = int(sys.argv[1]), int(sys.argv[2])
TAG = sys.argv[3] if len(sys.argv) > 3 else "fix2"
pre = json.load(open(os.path.join(HERE, "bandfix_0.json")))
pts = [(p["beta"], p["y"]) for p in pre]
mine = [p for n, p in enumerate(pts) if n % nslots == slot]

out = []
t0 = time.time()
for beta, y in mine:
    signal.signal(signal.SIGALRM, _a)
    signal.alarm(PER_POINT)
    try:
        r = X.certified_census(beta, y, N=180)
    except _T:
        out.append({"beta": beta, "y": y, "timeout": PER_POINT})
    except Exception as e:                                       # noqa
        out.append({"beta": beta, "y": y, "error": str(e)})
    else:
        out.append({
            "beta": beta, "y": y, "census": r["census"],
            "n_torque_roots": r["n_torque_roots"],
            "torque_t": [q["t_mid"] for q in r["torque_roots"]],
            "n_separated": r["n_separated"],
            "n_spurious_sep": sum(1 for s in r["separated"]
                                  if s["schur"] == "spurious"),
            "sep_D": [round(s["D_mid"], 6) for s in r["separated"]],
            "sep_schur": [s["schur"] for s in r["separated"]],
            "certified": r["certified_modulo_locator_completeness"],
            "torque_undecided": len(r["torque_undecided"]),
            "locator_failed": len(r["locator_failed_certification"]),
        })
    finally:
        signal.alarm(0)
    with open(os.path.join(HERE, "band%s_%d.json" % (TAG, slot)), "w") as fh:
        json.dump(out, fh)
print("slot %d: %d points, %.0fs" % (slot, len(out), time.time() - t0))
