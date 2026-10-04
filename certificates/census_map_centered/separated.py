"""
TASK 4: is any SEPARATED critical family of the CENTERED model with MIXED-SIGN
teacher masses a local minimum?

Model.  Teacher = two atoms, directions (0, beta), signed masses (s0, s1).
Students = two atoms, directions (th0, th1), signed weights (c0, c1).
Loss (up to a teacher-only constant)

    L(c,th) = 1/2 [ sum_ij c_i c_j phi(th_i-th_j)
                    - 2 sum_ij c_i s_j phi(th_i-beta_j) ]

with the CENTERED kernel phi = phiCos, phi(0) = kap = pi/2, phi' = -H.

Criticality (separated: th0 != th1 mod pi, c0 c1 != 0):
    radial:  kap c0 + phi(D) c1 = P(th0),   phi(D) c0 + kap c1 = P(th1)
    torque:  A(th_i) = 0, i.e. c0 = Trq(th1)/H(D), c1 = -Trq(th0)/H(D)
with D = th0-th1, P(t) = s0 phi(t)+s1 phi(t-beta), Trq(t) = -(s0 H(t)+s1 H(t-beta)).

Second variation.  Hessian in (c0,c1,th0,th1) is [[Gram, B],[B^T, A]] with
    Gram = [[kap, phi(D)],[phi(D), kap]]          (PD iff D != 0 mod pi)
    B    = [[0, c1 H(D)],[-c0 H(D), 0]]
    A00  = c0 Tau(th0) - c0 c1 H'(D),  A11 = c1 Tau(th1) - c0 c1 H'(D),
    A01  = c0 c1 H'(D),   Tau(t) = P(t) - 2 W(t),  W(t)=s0|sin t|+s1|sin(t-beta)|
The point is a strict local minimum iff the Schur complement A - B^T Gram^-1 B
is positive definite, i.e. (multiplying by massDet = kap^2-phi(D)^2 > 0)

    T00 = massDet*A00 - kap*H(D)^2*c0^2 > 0   and   det T > 0.

This is exactly widgets.js `schurLabel`, which the shipped classification uses
for the NONCENTERED separated stratum; the centered mixed-sign stratum lists no
separated row at all, which is the gap this script probes.

Method.  float64 numpy scan to LOCATE candidates, float64 Newton to polish,
then EVERY surviving candidate is re-solved and re-typed with mpmath at 60
decimal digits.  No verdict is taken from float64.  Configurations with a
near-singular student Gram (|D mod pi| small) are screened out and reported
separately -- that is where float64 manufactures nonsense.
"""
import sys
import json
import math
import numpy as np
import mpmath

PI = math.pi
KAP = PI / 2

# ---------------------------------------------------------------- float layer


def fmod(x, m=PI):
    return np.mod(x, m)


def f_phi(t):
    x = fmod(t)
    return (PI / 2 - x) * np.cos(x) + np.sin(x)


def f_H(t):
    x = fmod(t)
    return (PI / 2 - x) * np.sin(x)


def f_dH(t):
    x = fmod(t)
    return (PI / 2 - x) * np.cos(x) - np.sin(x)


def f_P(b, s0, s1, t):
    return s0 * f_phi(t) + s1 * f_phi(t - b)


def f_W(b, s0, s1, t):
    return s0 * np.abs(np.sin(t)) + s1 * np.abs(np.sin(t - b))


def f_Trq(b, s0, s1, t):
    return -(s0 * f_H(t) + s1 * f_H(t - b))


def residual(b, s0, s1, t0, t1):
    """Criticality residual in the WELL-CONDITIONED direction.

    Masses come from the RADIAL rows (divide by massDet = kap^2 - phi(D)^2,
    which vanishes only at coincidence D = 0 mod pi), and the two TORQUE rows
    are the residual.  The earlier arrangement -- masses from the torque rows,
    dividing by H(D) -- is singular on the whole quarter locus D = pi/2 mod pi,
    which is exactly where an extra separated family lives; that division is why
    the shipped separatedScan needs its |H(D)| > 5e-3 guard and cannot see it.
    """
    D = t0 - t1
    fd = f_phi(D)
    det = KAP * KAP - fd * fd
    P0 = f_P(b, s0, s1, t0)
    P1 = f_P(b, s0, s1, t1)
    with np.errstate(divide="ignore", invalid="ignore"):
        c0 = (KAP * P0 - fd * P1) / det
        c1 = (KAP * P1 - fd * P0) / det
        Hd = f_H(D)
        r0 = s0 * f_H(t0) + s1 * f_H(t0 - b) - c1 * Hd
        r1 = s0 * f_H(t1) + s1 * f_H(t1 - b) + c0 * Hd
    return r0, r1, c0, c1


def scan_teacher(b, s0, s1, N=240, hmin=2e-3, rtol=0.35):
    """Grid + VECTORISED Newton; returns deduped float candidates.

    Only cells that are local minima of the residual norm are polished, and the
    polish runs on all of them at once, which is what makes a full (beta,y)
    sweep affordable.  Everything here is float64 and is a LOCATOR only: every
    survivor is re-solved and re-typed with mpmath at 60 digits downstream.
    """
    g = (np.arange(N) + 0.5) * PI / N
    T0, T1 = np.meshgrid(g, g, indexing="ij")
    r0, r1, _, _ = residual(b, s0, s1, T0, T1)
    Dm = np.mod(T0 - T1, PI)
    gap = np.minimum(Dm, PI - Dm)
    bad = ~np.isfinite(r0) | ~np.isfinite(r1) | (gap < hmin)
    score = np.where(bad, np.inf, np.abs(r0) + np.abs(r1))
    S = np.pad(score, 1, constant_values=np.inf)
    loc = np.ones_like(score, dtype=bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            loc &= score <= S[1 + di:1 + di + N, 1 + dj:1 + dj + N]
    # two passes: local minima of the residual norm, plus EVERY cell under a
    # tight absolute threshold, so a narrow basin cannot fall between samples
    idx = np.argwhere((loc & (score < rtol)) | (score < 0.01))
    if not len(idx):
        return []
    x0 = T0[idx[:, 0], idx[:, 1]].copy()
    x1 = T1[idx[:, 0], idx[:, 1]].copy()
    live = np.ones(len(x0), dtype=bool)
    h = 1e-7
    for _ in range(70):
        f0, f1, _, _ = residual(b, s0, s1, x0, x1)
        a0, a1, _, _ = residual(b, s0, s1, x0 + h, x1)
        b0, b1, _, _ = residual(b, s0, s1, x0, x1 + h)
        j00, j01 = (a0 - f0) / h, (b0 - f0) / h
        j10, j11 = (a1 - f1) / h, (b1 - f1) / h
        det = j00 * j11 - j01 * j10
        ok = np.isfinite(det) & (np.abs(det) > 1e-13) & np.isfinite(f0) & np.isfinite(f1)
        d0 = np.where(ok, (-f0 * j11 + f1 * j01) / np.where(ok, det, 1), 0.0)
        d1 = np.where(ok, (-f1 * j00 + f0 * j10) / np.where(ok, det, 1), 0.0)
        step = np.clip(np.abs(d0) + np.abs(d1), 1e-300, None)
        damp = np.minimum(1.0, 0.5 / step)
        x0 = x0 + 0.9 * damp * d0
        x1 = x1 + 0.9 * damp * d1
        live &= ok
    f0, f1, cc0, cc1 = residual(b, s0, s1, x0, x1)
    good = live & np.isfinite(f0) & np.isfinite(f1) & (np.abs(f0) + np.abs(f1) < 1e-9)
    out = []
    for k in np.flatnonzero(good):
        t0, t1 = float(np.mod(x0[k], PI)), float(np.mod(x1[k], PI))
        c0, c1 = float(cc0[k]), float(cc1[k])
        dm = math.fmod(t0 - t1, PI)
        if dm < 0:
            dm += PI
        if min(dm, PI - dm) < hmin:
            continue
        latd = lambda t: min(min(t % PI, PI - t % PI),
                             min((t - b) % PI, PI - (t - b) % PI))
        if latd(t0) < 1e-6 and latd(t1) < 1e-6:
            continue                                   # the exact fit
        if t1 < t0:
            t0, t1, c0, c1 = t1, t0, c1, c0
        if not any(abs(q[0] - t0) < 1e-6 and abs(q[1] - t1) < 1e-6 for q in out):
            out.append((t0, t1, c0, c1))
    return out


# ------------------------------------------------------------ mpmath layer

def mp_setup(dps=60):
    mpmath.mp.dps = dps


def m_phi(t):
    x = mpmath.fmod(t, mpmath.pi)
    if x < 0:
        x += mpmath.pi
    return (mpmath.pi / 2 - x) * mpmath.cos(x) + mpmath.sin(x)


def m_H(t):
    x = mpmath.fmod(t, mpmath.pi)
    if x < 0:
        x += mpmath.pi
    return (mpmath.pi / 2 - x) * mpmath.sin(x)


def m_dH(t):
    x = mpmath.fmod(t, mpmath.pi)
    if x < 0:
        x += mpmath.pi
    return (mpmath.pi / 2 - x) * mpmath.cos(x) - mpmath.sin(x)


def m_res(b, s0, s1, t0, t1):
    """mpmath mirror of the well-conditioned residual (masses from the radial
    rows, torque rows as the residual)."""
    kap = mpmath.pi / 2
    D = t0 - t1
    fd = m_phi(D)
    Hd = m_H(D)
    det = kap * kap - fd * fd
    P0 = s0 * m_phi(t0) + s1 * m_phi(t0 - b)
    P1 = s0 * m_phi(t1) + s1 * m_phi(t1 - b)
    c0 = (kap * P0 - fd * P1) / det
    c1 = (kap * P1 - fd * P0) / det
    r0 = s0 * m_H(t0) + s1 * m_H(t0 - b) - c1 * Hd
    r1 = s0 * m_H(t1) + s1 * m_H(t1 - b) + c0 * Hd
    return [r0, r1], c0, c1


def m_solve(b, s0, s1, t0, t1, iters=80):
    x = [mpmath.mpf(t0), mpmath.mpf(t1)]
    h = mpmath.mpf(10) ** (-mpmath.mp.dps // 2)
    for _ in range(iters):
        f, _, _ = m_res(b, s0, s1, x[0], x[1])
        fa, _, _ = m_res(b, s0, s1, x[0] + h, x[1])
        fb, _, _ = m_res(b, s0, s1, x[0], x[1] + h)
        J = mpmath.matrix([[(fa[0] - f[0]) / h, (fb[0] - f[0]) / h],
                           [(fa[1] - f[1]) / h, (fb[1] - f[1]) / h]])
        det = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
        if abs(det) < mpmath.mpf(10) ** (-mpmath.mp.dps + 5):
            break
        d0 = (-f[0] * J[1, 1] + f[1] * J[0, 1]) / det
        d1 = (-f[1] * J[0, 0] + f[0] * J[1, 0]) / det
        x = [x[0] + d0, x[1] + d1]
        if abs(d0) + abs(d1) < mpmath.mpf(10) ** (-mpmath.mp.dps + 8):
            break
    f, c0, c1 = m_res(b, s0, s1, x[0], x[1])
    return x, c0, c1, abs(f[0]) + abs(f[1])


def m_type(b, s0, s1, t0, t1, c0, c1):
    """Reduced (Schur) angular Hessian and the full 4x4 Hessian spectrum."""
    kap = mpmath.pi / 2
    D = t0 - t1
    Hd, Hp, Fd = m_H(D), m_dH(D), m_phi(D)
    P = lambda t: s0 * m_phi(t) + s1 * m_phi(t - b)
    Wg = lambda t: s0 * abs(mpmath.sin(t)) + s1 * abs(mpmath.sin(t - b))
    Tau = lambda t: P(t) - 2 * Wg(t)
    A00 = -c0 * c1 * Hp + c0 * Tau(t0)
    A11 = -c0 * c1 * Hp + c1 * Tau(t1)
    A01 = c0 * c1 * Hp
    massDet = kap * kap - Fd * Fd
    T00 = massDet * A00 - kap * Hd * Hd * c0 * c0
    T11 = massDet * A11 - kap * Hd * Hd * c1 * c1
    T01 = massDet * A01 - Fd * Hd * Hd * c0 * c1
    detT = T00 * T11 - T01 * T01
    Hess = mpmath.matrix([[kap, Fd, 0, c1 * Hd],
                          [Fd, kap, -c0 * Hd, 0],
                          [0, -c0 * Hd, A00, A01],
                          [c1 * Hd, 0, A01, A11]])
    try:
        ev = mpmath.eigsy(Hess, eigvals_only=True)
        evs = sorted([mpmath.mpf(ev[i]) for i in range(4)])
    except Exception:
        evs = None
    return {"massDet": massDet, "T00": T00, "T11": T11, "T01": T01,
            "detT": detT, "Hd": Hd, "eig": evs}


def run(nb=48, ny=32, mixed=True, dps=60, out=None):
    mp_setup(dps)
    results = []
    stats = {"teachers": 0, "candidates": 0, "verified": 0,
             "minima": 0, "psd_boundary": 0, "screened_near_coincident": 0,
             "float_flips": 0}
    for ib in range(nb):
        beta = PI * (ib + 0.5) / nb
        for iy in range(ny):
            y = (iy + 0.5) / ny if mixed else -1 + (iy + 0.5) / ny
            psi = (y + 1) * PI / 2
            s0, s1 = math.sin(psi), math.cos(psi)
            stats["teachers"] += 1
            cands = scan_teacher(beta, s0, s1)
            stats["candidates"] += len(cands)
            for (t0, t1, c0, c1) in cands:
                mb = mpmath.mpf(beta)
                mpsi = (mpmath.mpf(y) + 1) * mpmath.pi / 2
                ms0, ms1 = mpmath.sin(mpsi), mpmath.cos(mpsi)
                x, C0, C1, res = m_solve(mb, ms0, ms1, t0, t1)
                if res > mpmath.mpf(10) ** (-dps + 12):
                    continue
                Dm = mpmath.fmod(x[0] - x[1], mpmath.pi)
                if Dm < 0:
                    Dm += mpmath.pi
                gap = min(Dm, mpmath.pi - Dm)
                if gap < mpmath.mpf("1e-4"):
                    stats["screened_near_coincident"] += 1
                    continue
                stats["verified"] += 1
                ty = m_type(mb, ms0, ms1, x[0], x[1], C0, C1)
                is_min = (ty["T00"] > 0 and ty["detT"] > 0)
                near = (abs(ty["T00"]) < mpmath.mpf("1e-30")
                        or abs(ty["detT"]) < mpmath.mpf("1e-30"))
                # float64 comparison, to record where float64 would disagree
                fty = float(ty["T00"]), float(ty["detT"])
                if is_min:
                    stats["minima"] += 1
                if near:
                    stats["psd_boundary"] += 1
                if is_min or near:
                    results.append({
                        "beta": float(beta), "y": float(y),
                        "s0": float(s0), "s1": float(s1),
                        "t0": mpmath.nstr(x[0], 30), "t1": mpmath.nstr(x[1], 30),
                        "c0": mpmath.nstr(C0, 30), "c1": mpmath.nstr(C1, 30),
                        "gap": mpmath.nstr(gap, 20),
                        "T00": mpmath.nstr(ty["T00"], 25),
                        "T11": mpmath.nstr(ty["T11"], 25),
                        "detT": mpmath.nstr(ty["detT"], 25),
                        "massDet": mpmath.nstr(ty["massDet"], 20),
                        "eig": [mpmath.nstr(e, 20) for e in ty["eig"]] if ty["eig"] else None,
                        "is_min": bool(is_min), "near_boundary": bool(near),
                        "residual": mpmath.nstr(res, 10)})
        print("  beta row %d/%d done  (verified %d, minima %d)"
              % (ib + 1, nb, stats["verified"], stats["minima"]), flush=True)
    payload = {"mixed_sign": mixed, "nb": nb, "ny": ny, "dps": dps,
               "stats": stats, "hits": results}
    if out:
        with open(out, "w") as fh:
            json.dump(payload, fh, indent=1)
    return payload


if __name__ == "__main__":
    mixed = "--positive" not in sys.argv
    nb = int(sys.argv[sys.argv.index("--nb") + 1]) if "--nb" in sys.argv else 48
    ny = int(sys.argv[sys.argv.index("--ny") + 1]) if "--ny" in sys.argv else 32
    out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else None
    p = run(nb=nb, ny=ny, mixed=mixed, out=out)
    print(json.dumps(p["stats"], indent=1))
    print("hits:", len(p["hits"]))
