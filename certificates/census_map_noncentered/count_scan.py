"""
The SEPARATED-FAMILY COUNT scan (the counting-slack measurement of the
count-matching route).

MEASUREMENT, not certification: for a dense grid of teachers per proposed
arrangement label, run the certificate's own float locator (census_cert.
locate_separated, mirrored here in numpy for speed and VALIDATED against the
reference at the face representatives) and record the number of distinct
separated families found.  The per-face max/min is the target any provable
upper bound must MATCH at the witnesses for count-matching completeness
(found >= bound  =>  found = all).

Validity domain (SKILL 1v): float64 Newton on the well-conditioned radial
clearing; the reference locator's own instrument limits apply (families at
gap D < dmin invisible; near-fold pairs may need denser seeding - checked
here by doubling N at every point where the count changes between N and 2N).

Output: separated_count_scan.json
"""
import json
import math
import sys

import numpy as np

import noncentered as J
import census_cert as C

TP = 2 * math.pi
PI = math.pi


# --------------------------------------------------------------------------
# numpy mirror of the widgets.js kernel atoms (float locator only)
# --------------------------------------------------------------------------

def v_phiJ(t):
    x = np.mod(t, TP)
    x = np.where(x > PI, TP - x, x)
    return (PI - x) * np.cos(x) + np.sin(x)


def v_HJ(t):
    x = np.mod(t, TP)
    return np.where(x <= PI, (PI - x) * np.sin(x), (x - PI) * np.sin(x))


def v_res(beta, s0, s1, t0, t1):
    """(G0, G1) of census_cert._f_res, vectorized; nan where massDet ~ 0."""
    D = t0 - t1
    phiD = v_phiJ(D)
    hD = v_HJ(D)
    M = PI * PI - phiD * phiD
    bad = np.abs(M) < 1e-12
    Msafe = np.where(bad, 1.0, M)
    P0 = s0 * v_phiJ(t0) + s1 * v_phiJ(t0 - beta)
    P1 = s0 * v_phiJ(t1) + s1 * v_phiJ(t1 - beta)
    A0 = -(s0 * v_HJ(t0) + s1 * v_HJ(t0 - beta))
    A1 = -(s0 * v_HJ(t1) + s1 * v_HJ(t1 - beta))
    c0 = (PI * P0 - phiD * P1) / Msafe
    c1 = (PI * P1 - phiD * P0) / Msafe
    G0 = c1 * hD + A0
    G1 = c0 * hD - A1
    G0 = np.where(bad, np.nan, G0)
    G1 = np.where(bad, np.nan, G1)
    return G0, G1


def locate_np(beta, s0, s1, N=260, dmin=1e-4, dvals=None):
    """numpy mirror of census_cert.locate_separated: same seeds, same Newton
    (finite-difference Jacobian, h = 1e-7, tol 1e-13, 60 iterations), same
    folding, dedup and dmin filter."""
    if dvals is None:
        dvals = C.d_grid(dmin)
    dvals = np.asarray(dvals, dtype=float)
    t1s = TP * (np.arange(N) + 0.5) / N
    T1, DV = np.meshgrid(t1s, dvals, indexing="ij")
    x0 = (T1 + DV).ravel()
    x1 = T1.ravel()
    alive = np.ones(x0.shape, dtype=bool)
    done = np.zeros(x0.shape, dtype=bool)
    h = 1e-7
    for _ in range(60):
        idx = alive & ~done
        if not idx.any():
            break
        f0, f1 = v_res(beta, s0, s1, x0[idx], x1[idx])
        bad = np.isnan(f0) | np.isnan(f1)
        conv = (np.abs(f0) + np.abs(f1)) < 1e-13
        fa0, fa1 = v_res(beta, s0, s1, x0[idx] + h, x1[idx])
        fb0, fb1 = v_res(beta, s0, s1, x0[idx], x1[idx] + h)
        bad |= np.isnan(fa0) | np.isnan(fb0) | np.isnan(fa1) | np.isnan(fb1)
        j00 = (fa0 - f0) / h
        j01 = (fb0 - f0) / h
        j10 = (fa1 - f1) / h
        j11 = (fb1 - f1) / h
        det = j00 * j11 - j01 * j10
        bad |= np.abs(det) < 1e-14
        detsafe = np.where(np.abs(det) < 1e-300, 1.0, det)
        dx0 = (-f0 * j11 + f1 * j01) / detsafe
        dx1 = (-f1 * j00 + f0 * j10) / detsafe
        sub_alive = ~bad & ~conv
        ii = np.where(idx)[0]
        x0[ii[sub_alive]] += dx0[sub_alive]
        x1[ii[sub_alive]] += dx1[sub_alive]
        alive[ii[bad]] = False
        done[ii[conv & ~bad]] = True
    found = []
    a_all, b_all = x0[done], x1[done]
    for a, b in zip(a_all, b_all):
        d = a - b
        if d < 0:
            a, b, d = b, a, -d
        d = d % TP
        b = b % TP
        if d > PI:
            d = TP - d
            b = (b + TP - d) % TP
        if d < dmin:
            continue
        if any(abs(b - q[0]) < 1e-6 and abs(d - q[1]) < 1e-6 for q in found):
            continue
        found.append((b, d))
    found.sort()
    return found


def count_families(beta, y, N=260, dmin=1e-4):
    """Distinct separated families (exact fit dropped), as the census counts
    them: unordered pair mod 2pi, dedup 1e-5, both-students-on-lattice
    dropped."""
    psi = (y + 1) * PI / 2
    s0, s1 = math.sin(psi), math.cos(psi)
    cands = locate_np(beta, s0, s1, N=N, dmin=dmin)

    def on_lattice(t):
        for c in (0.0, beta % TP):
            d = (t - c) % TP
            if min(d, TP - d) < 1e-5:
                return True
        return False

    uniq = []
    for (b, d) in cands:
        th0, th1 = (b + d) % TP, b
        k = tuple(sorted((round(th0, 7), round(th1, 7))))
        if any(abs(k[0] - u[0]) < 1e-5 and abs(k[1] - u[1]) < 1e-5
               for u in uniq):
            continue
        uniq.append(k)
    fams = [k for k in uniq
            if not (on_lattice(k[0]) and on_lattice(k[1]))]
    return fams


# --------------------------------------------------------------------------
# validation against the reference locator (SKILL 1v)
# --------------------------------------------------------------------------

def validate(points):
    """The numpy locator must reproduce the reference family list at the
    given (beta, y) points, to 1e-6 in (th1, D)."""
    out = []
    for (beta, y) in points:
        psi = (y + 1) * PI / 2
        s0, s1 = math.sin(psi), math.cos(psi)
        ref = C.locate_separated(beta, s0, s1, N=64, dmin=1e-4)
        got = locate_np(beta, s0, s1, N=64, dmin=1e-4)
        ok = len(ref) == len(got) and all(
            abs(r[0] - g[0]) < 1e-6 and abs(r[1] - g[1]) < 1e-6
            for r, g in zip(ref, got))
        out.append({"beta": beta, "y": y, "n_ref": len(ref),
                    "n_np": len(got), "match": bool(ok)})
    return out


# --------------------------------------------------------------------------
# the per-face grids
# --------------------------------------------------------------------------

def linspace(a, b, n):
    return [a + (b - a) * (i + 0.5) / n for i in range(n)]


def face_grids():
    """Finite grids keyed by the proposed labels of CERTIFICATE.md section 7.

    The grids sample the fundamental domain `beta in (0,pi)`. Their placement
    relative to sampled candidate-transition brackets does not certify global
    faces, walls, or constancy between teachers.
    """
    g = {}
    g["F1"] = [(b, y) for b in linspace(0.1, 3.13, 16)
               for y in linspace(0.25, 0.75, 8)]
    g["F1_near_collar"] = [(b, y) for b in linspace(0.1, 3.13, 12)
                           for y in [0.03, 0.06, 0.10, 0.90, 0.94, 0.97]]
    g["collar"] = [(b, y) for b in [0.05, 0.10, 0.15, 0.30, 0.50, 0.75]
                   for y in [0.002, 0.004, 0.008]]
    g["F2"] = [(b, y) for b in linspace(0.1, 2.5, 12)
               for y in [-0.01, -0.02, -0.03, -0.05]]
    g["F3"] = [(b, y) for b in [3.05, 3.08, 3.10, 3.12, 3.13]
               for y in [-0.05, -0.07, -0.09, -0.11]]
    g["F4"] = [(b, y) for b in linspace(0.1, 2.2, 10)
               for y in linspace(-0.75, -0.25, 6)]
    g["F5"] = [(b, y) for b in linspace(2.3, 3.13, 10)
               for y in linspace(-0.7, -0.3, 6)]
    g["F6"] = [(b, y) for b in linspace(0.1, 2.5, 12)
               for y in [-0.99, -0.98, -0.97, -0.95]]
    g["F7"] = [(b, y) for b in [3.05, 3.08, 3.10, 3.12, 3.13]
               for y in [-0.95, -0.93, -0.91, -0.89]]
    # Candidate-transition neighbourhoods: sample across pointwise
    # count-change brackets. The scan does not certify that a fold pair is born
    # or that an intervening wall exists.
    g["wall_mixed"] = [(b, y) for b in [0.3, 0.75, 1.5, 2.2, 3.0, 3.13]
                       for y in [-0.06, -0.08, -0.10, -0.12, -0.14, -0.16]]
    return g


def main():
    reps = [(0.75, 0.40), (0.75, -0.02), (3.10, -0.06), (0.75, -0.40),
            (2.60, -0.50), (0.75, -0.98), (3.10, -0.94), (0.15, 0.004)]
    val = validate(reps)
    print("validation:", json.dumps(val, indent=1))
    if not all(v["match"] for v in val):
        print("VALIDATION FAILED - not scanning")
        json.dump({"validation": val, "scan": None},
                  open("separated_count_scan.json", "w"), indent=1)
        sys.exit(1)

    grids = face_grids()
    scan = {}
    for face, pts in grids.items():
        rows = []
        for (b, y) in pts:
            fams = count_families(b, y, N=260)
            n = len(fams)
            n2 = n
            # seed-doubling check at every point (not only at maxima):
            # a count that CHANGES under doubling is an instrument limit
            fams2 = count_families(b, y, N=520)
            n2 = len(fams2)
            rows.append({"beta": round(b, 6), "y": round(y, 6),
                         "n": n, "n_2N": n2,
                         "D": [round(f[0] - f[1] if f[0] > f[1] else 0, 6)
                               for f in []],
                         "stable": bool(n == n2)})
            print(face, round(b, 3), round(y, 3), n, n2, flush=True)
        ns = [r["n_2N"] for r in rows]
        scan[face] = {"points": rows, "min": min(ns), "max": max(ns),
                      "all_stable": all(r["stable"] for r in rows)}
        print(face, "min", min(ns), "max", max(ns))
    json.dump({"validation": val, "scan": scan},
              open("separated_count_scan.json", "w"), indent=1)
    print("wrote separated_count_scan.json")


if __name__ == "__main__":
    main()
