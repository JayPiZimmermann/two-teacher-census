"""Driver for TASK 2: certified covering of the beta axis with per-beta root
counts of Wtau and Wpot, plus enclosures of the torque-lens branch in (beta,y).

Usage: python3 run_coincident.py <interval_index 0..3> <n_boxes>
Intervals: 0=(0,beta*_1) 1=(beta*_1,pi) 2=(pi,beta*_2) 3=(beta*_2,2pi)
Writes coincident_<i>.json
"""
import json
import os
import sys
import time

import mpmath
from mpmath import iv

import noncentered as J
import certify_coincident as C

HERE = os.path.dirname(os.path.abspath(__file__))
J.set_prec(C.PREC)

u_star, B1, B2 = C.certify_beta_star()
PIv = J.PI_IV()
BOUNDS = [(iv.mpf(0), B1), (B1, PIv), (PIv, B2), (B2, 2 * PIv)]
NAMES = ["(0, beta*_1)", "(beta*_1, pi)", "(pi, beta*_2)", "(beta*_2, 2pi)"]
MINW_BETA = mpmath.mpf("1e-9")


def y_of_root_box(pc, mbox):
    d = pc["d"]
    p = mbox - d / 2
    q = mbox + d / 2
    et = 1 if pc["nt"] % 2 == 0 else -1
    eu = 1 if pc["nu"] % 2 == 0 else -1
    ht = (et * p) * iv.sin(p)
    hu = (eu * q) * iv.sin(q)
    s0, s1 = -hu, ht
    return s0, s1, J.y_from_s(s0, s1)


def certify_box(blo, bhi):
    B = iv.mpf([blo, bhi])
    r = C.certify_beta_box(B)
    if not r.get("ok"):
        return r
    # lens roots -> y
    ys = []
    for pc, rec in zip(J.pieces_of(B), r["pieces"]):
        for (x0, x1) in rec.get("Wtau_root_boxes", []):
            mb = C.refine_root("tau", pc["d"], iv.sin(B), x0, x1)
            s0, s1, y = y_of_root_box(pc, mb)
            tb = J.t_from_m(pc, mb)
            ys.append({
                "piece_t": [float(J.mid(pc["t_lo"])), float(J.mid(pc["t_hi"]))],
                "t": [mpmath.nstr(J.lo(tb), 22), mpmath.nstr(J.hi(tb), 22)],
                "m": [mpmath.nstr(J.lo(mb), 22), mpmath.nstr(J.hi(mb), 22)],
                "y": None if y is None else [mpmath.nstr(J.lo(y), 22),
                                             mpmath.nstr(J.hi(y), 22)],
                "y_width": None if y is None else J.width(y),
            })
    r["lens_roots"] = ys
    return r


def run(i, nboxes):
    lo_, hi_ = BOUNDS[i]
    a, b = J.hi(lo_), J.lo(hi_)
    boxes = []
    for k in range(nboxes):
        boxes.append((a + (b - a) * k / nboxes, a + (b - a) * (k + 1) / nboxes))
    done = []
    uncovered = []
    stack = list(reversed(boxes))
    t0 = time.time()
    while stack:
        x0, x1 = stack.pop()
        try:
            r = certify_box(x0, x1)
        except Exception as e:            # noqa
            r = {"ok": False, "reason": "exception: %s" % e}
        if r.get("ok"):
            done.append({
                "beta": [mpmath.nstr(x0, 22), mpmath.nstr(x1, 22)],
                "width": float(x1 - x0),
                "Wtau_total": r["Wtau_total"], "Wpot_total": r["Wpot_total"],
                "Wtau_interior": r["Wtau_interior"],
                "Wpot_interior": r["Wpot_interior"],
                "lens_roots": r["lens_roots"],
            })
        else:
            if x1 - x0 < MINW_BETA:
                uncovered.append({"beta": [mpmath.nstr(x0, 22), mpmath.nstr(x1, 22)],
                                  "width": float(x1 - x0),
                                  "reason": r.get("reason") or r.get("unresolved")})
            else:
                m = (x0 + x1) / 2
                stack.append((m, x1))
                stack.append((x0, m))
    cov = sum(mpmath.mpf(d["width"]) for d in done)
    out = {
        "interval": NAMES[i],
        "beta_range": [mpmath.nstr(a, 25), mpmath.nstr(b, 25)],
        "requested_boxes": nboxes,
        "certified_boxes": len(done),
        "measure_of_range": float(b - a),
        "measure_covered": float(cov),
        "measure_uncovered_inside_range": float((b - a) - cov),
        "collar_at_lower_end": float(J.hi(lo_) - J.lo(lo_)),
        "collar_at_upper_end": float(J.hi(hi_) - J.lo(hi_)),
        "uncovered_boxes": uncovered,
        "boxes": done,
        "seconds": time.time() - t0,
    }
    with open(os.path.join(HERE, "coincident_%d.json" % i), "w") as fh:
        json.dump(out, fh)
    print("interval %d %s: %d boxes, covered %.15f of %.15f, uncovered %d, %.1fs"
          % (i, NAMES[i], len(done), float(cov), float(b - a), len(uncovered),
             out["seconds"]))


if __name__ == "__main__":
    run(int(sys.argv[1]), int(sys.argv[2]))
