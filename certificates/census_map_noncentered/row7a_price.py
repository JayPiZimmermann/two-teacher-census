"""PROSPECTING (not evidence): price the row-7a exclusion sweep on a RESTRICTED
student band, before any certificate is emitted.

The full sweep over `(beta, y, s, D)` with `D` down to the census gap floor is
already measured as unreachable.  This probe asks the narrower question the
scope statement needs: on a band `D >= d0` where the SCHUR test can fire, how
many nodes does the exclusion actually cost, and does the bisection TERMINATE?

It reuses `foldfree_sweep.verdict` and `foldfree_sweep.child` unchanged, so
the verdicts priced here are the verdicts a certificate would carry.

Usage: python3 row7a_price.py <b0> <b1> <y0> <y1> <d0> <d1> <budget>
"""
import sys
import time

import foldfree_sweep as S


def main():
    b0, b1, y0, y1, d0, d1 = (float(sys.argv[i]) for i in range(1, 7))
    budget = int(sys.argv[7]) if len(sys.argv) > 7 else 200000
    two_pi = 6.283185307179587
    spec = {"b0": b0, "b1": b1, "y0": y0, "y1": y1, "seam": 0.137,
            "delta": d0, "dmax_pad": 0.05, "fitr": 0.05, "minw": 1e-9,
            "n_parts": 1, "order": "map_first"}
    root = (b0, b1, y0, y1, spec["seam"], spec["seam"] + two_pi, d0, d1)
    print("root beta=[%g,%g] y=[%g,%g] s=[%g,%g] D=[%g,%g] budget %d"
          % (root[0], root[1], root[2], root[3], root[4], root[5], root[6],
             root[7], budget), flush=True)
    stack = [(root, 0)]
    counts = {"MAP": 0, "SCHUR": 0, "FIT": 0, "UNDECIDED": 0}
    bydepth = {}
    n = 0
    live = 0.0
    t0 = time.time()
    maxdepth = 0
    while stack:
        bx, dep = stack.pop()
        n += 1
        maxdepth = max(maxdepth, dep)
        if n > budget:
            print("BUDGET EXHAUSTED: %d nodes, stack %d, alive volume so far "
                  "%.3e" % (n, len(stack), live), flush=True)
            break
        v = S.verdict(spec, *bx)
        if v is None:
            w = max(bx[1] - bx[0], bx[3] - bx[2], bx[5] - bx[4],
                    bx[7] - bx[6])
            if w < spec["minw"]:
                counts["UNDECIDED"] += 1
                live += ((bx[1] - bx[0]) * (bx[3] - bx[2]) * (bx[5] - bx[4])
                         * (bx[7] - bx[6]))
                continue
            bydepth[dep] = bydepth.get(dep, 0) + 1
            stack.append((S.child(bx, 1), dep + 1))
            stack.append((S.child(bx, 0), dep + 1))
            continue
        counts[v] += 1
        if n % 5000 == 0:
            print("  nodes %8d  map %7d schur %6d fit %5d undec %d  "
                  "stack %6d  depth %3d  %.0fs"
                  % (n, counts["MAP"], counts["SCHUR"], counts["FIT"],
                     counts["UNDECIDED"], len(stack), dep, time.time() - t0),
                  flush=True)
    secs = time.time() - t0
    print("nodes %d  map %d schur %d fit %d undecided %d  maxdepth %d  "
          "%.1f s  (%.1f ms/node)"
          % (n, counts["MAP"], counts["SCHUR"], counts["FIT"],
             counts["UNDECIDED"], maxdepth, secs, 1000.0 * secs / max(n, 1)))
    print("SPLIT COUNT BY DEPTH (survivors per level is the sweep's exponent):")
    ks = sorted(bydepth)
    for k in ks:
        if k % 4 == 0 or k == ks[-1]:
            print("  depth %3d : %8d splits" % (k, bydepth[k]))
    if not stack and not counts["UNDECIDED"]:
        print("COMPLETE: every leaf discharged")
        return 0
    print("INCOMPLETE")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
