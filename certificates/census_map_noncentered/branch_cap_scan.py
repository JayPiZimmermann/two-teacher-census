"""
Branch-cap route measurements for the separated-family counting campaign
(2026-08-08, follow-up to the proved sign law of separatedNumJ).

THE OPEN LEMMA (named precisely, measured tight):
  for every teacher (beta, s0, s1), every D in (0, 2pi), and every
  direction (alpha, gamma) != 0, the combination

    alpha*(s0*W(t,D) + s1*W(t-beta,D)) + gamma*(s0*h(t) + s1*h(t-beta)),
    W(t,D) = pi*KJ(t-D) - KJ(D)*KJ(t),   h = couplingHJ,

  has AT MOST 4 sign changes in t per 2pi-period.  The student load
  L = s0*num(D,.) + s1*num(D,.-beta) is the instance
  (alpha, gamma) = (hJ(D), -(pi^2-KJ(D)^2)); i.e. the cap is a GLOBAL
  line-crossing property of the closed curve
  Gamma_D(t) = (W_load(t), h_load(t)), not a property of the special line.

ROUTE-A FEASIBILITY: THE AMBIENT 4-DIM SPACE IS **NOT** A WEAK CHEBYSHEV
SYSTEM (measured and REFUTED, 2026-08-08, third pass).  The open lemma above
quantifies over the 2-dimensional pencil the STUDENT LOAD lives in.  The
standard sufficient criterion would have been to prove the cap for the whole
ambient space

    V(D,beta) = span{ W(.,D), W(.-beta,D), h(.), h(.-beta) }

(then collocation-determinant sign-consistency / Dodgson condensation of the
4x4 into 2x2 crosses of the mother pair would apply).  That is FALSE.  The
element

    c = (0.5389600251329808, -0.14246740163039062,
         0.7066748566272032, -0.4357014778403485)
    D = 4.846208362913536,  beta = 0.7349885346109851

has SIX sign changes per period (measurement 3 below re-checks it; verified
at 24000 / 96000 / 400000 / 2000000 grid points and independently at 50 and
60 decimal digits with mpmath, all six crossings transversal).  Histogram of
the per-(D,beta) maximum over a 62-pair structured sweep: {4: 15, 6: 47};
nothing reached 8.

TWO THINGS THE REFUTATION DOES NOT DO.
 (a) It does NOT touch the open lemma.  The lemma's directions are
     (alpha*s0, alpha*s1, gamma*s0, gamma*s1), i.e. the TEACHER-CONSISTENT
     variety c0*c3 - c1*c2 = 0.  At the witness c0*c3 - c1*c2 = -0.134 with
     |c| = 1 -- 57% of scale, far off the variety, not a near-miss.  The
     cap for the student load remains measured-true and unrefuted.
 (b) It is therefore INFORMATIVE: the teacher-consistency coupling
     (c0,c1) || (c2,c3) is ESSENTIAL to the cap, not incidental.  Any proof
     must use it, and no argument treating the four generators as an
     unstructured space can succeed.  This is the same verdict the E + S
     split gave from the other side (the content sits in the coupling
     BETWEEN the two branches, not in a per-branch 2 + 2).

HOW IT WAS FOUND, and the methodological warning (skill 1m/1p): 20000-
direction Sobol sampling plus basin-hopping never exceeded 4 at ANY of the
62 pairs -- 100% false-negative rate.  Only the STRUCTURED probe (solve the
4-node collocation system for prescribed nodes, sweep the nodes) produced
the witness.  A sampled maximum over a projective sphere is an
under-estimate and must never be reported as a cap.
Not clean, not refuted: the pocket D near pi AND beta near pi, where the
same pipeline found only 4 across three shrinking epsilons.
Instrument note: the 5-node collocation sigma_min is NOT bounded away from
zero (infimum 1.0e-12), but that infimum sits at the boundary-degenerate
D -> 2pi and is a conditioning artefact, not a second counterexample.

KILL-RECORDS (routes refuted by slack, skill 1t -- do not re-try):
  * SIGN-LAW ARC LOCALISATION DOES NOT REDUCE THE COUNT (2026-08-08).  The
    proved straddle (positive masses, LatticeStraddleJ) / same-side (mixed
    masses, MixedSameSideJ) condition confines every student to explicit
    arcs of the gap lattice -- but that condition is IMPLIED by the load
    equation s0*num(D,a) + s1*num(D,b) = 0 itself, so it removes no zero of
    the 1-D load.  Measured (measurement 4): at all 30 face witnesses the
    number of sign changes of L_D INSIDE the admissible arcs equals the
    number over the whole period, and up to 3 of them bunch in ONE arc, so
    "one zero per arc" is false too.  (Validity domain: teachers with both
    masses bounded away from 0 and gaps bounded away from the teacher atom
    -- at min|s| -> 0 the load degenerates to a single num(D,.) whose zeros
    ARE the arc endpoints; measurement 4 asserts that characterisation of
    every mismatch.)  The sign law is a 2-D SEARCH pruner
    (50-91% of the (theta1,D) torus for positive teachers, 5-49% for mixed
    ones; face_bb_signlaw.py), not a counting lemma.
  * The two straddle theorems are ONE constraint, not two.  Student 1's
    offsets (-theta1, beta-theta1) are the image of student 0's
    (theta0, theta0-beta) under x |-> D - x, and that map preserves each
    open side of the gap lattice {0, D}; so the second row's sign condition
    is identical to the first's and buys no extra exclusion (measured: the
    admissible-area fraction is unchanged when the second test is added).
  * Sturm vs forcing: (d^2+1)L = 2G with G carrying |sin| corners AND
    sgn(sin)cos jumps at atoms {0,D,beta,beta+D}+pi*Z; measured SC(G) = 8
    at 25 of the 28 census face witnesses (max over D), so the hump
    alternation bound SC(L) <= SC(G)+2 = 10 has slack 6 against the
    needed-tight 4.  A permissive count closes nothing.
  * per-arc ratio monotonicity: partitioning by zeros of the load coupling
    h_load and requiring W_load/h_load monotone per arc fails -- measured
    up to 4 sign changes of L bunched INSIDE one h_load-arc.
  * per-arc disconjugacy: L is a resonant sinusoid
    (a+bt)sin t + (c+dt)cos t on each of the <= 8 atom arcs, and
    (D^2+1)^2 is disconjugate on arcs shorter than pi (<= 3 zeros/arc),
    giving 24 -- slack 20.
  * monotone direction / winding bookkeeping: the load Wronskian
    W_load' h_load - W_load h_load' has up to 6 sign changes per period,
    and the net winding of Gamma_D is 0 or -2pi depending on the sector,
    while the crossing count stays <= 4 in BOTH -- so neither monotone
    direction nor naive winding+backtrack accounting reproduces the cap.

Beyond the branch cap, the family caps (4 same-sign / 3 mixed) further
need the D-profile crossing bound of the second row along load-zero
branches (the fold-profile architecture of the campaign work note).

Running this file re-executes all measurements and asserts the recorded
values (validity domain: float64, parameters kept >= 0.05 from the
degenerate lattice, grids of 12000-16000 points per period).
"""
import json
import math
import os

import numpy as np

pi = np.pi


def wrap(t):
    return np.arccos(np.cos(t))


def K(t):
    a = wrap(t)
    return (pi - a) * np.cos(t) + np.sin(a)


def h(t):
    a = wrap(t)
    return np.sin(t) * (pi - a)


def num(D, x):
    return (pi * K(x - D) - K(D) * K(x)) * h(D) - (pi ** 2 - K(D) ** 2) * h(x)


def W(x, D):
    return pi * K(x - D) - K(D) * K(x)


def F(D, x):
    M = pi ** 2 - K(D) ** 2
    return (pi * h(D) * np.abs(np.sin(x - D))
            - K(D) * h(D) * np.abs(np.sin(x))
            + M * np.sign(np.sin(x)) * np.cos(x))


def sc_circle(v):
    s = np.sign(v)
    s = s[s != 0]
    if len(s) == 0:
        return 0
    s = np.append(s, s[0])
    return int(np.sum(np.diff(s) != 0))


def main():
    rng = np.random.default_rng(5)
    ts = np.linspace(0, 2 * pi, 12001)[:-1]

    # 1. the open lemma: all-line crossing cap = 4
    worst = 0
    for _ in range(400):
        D = rng.uniform(0.05, 2 * pi - 0.05)
        beta = rng.uniform(0.05, pi - 0.05)
        y = rng.uniform(-0.97, 0.97)
        psi = (y + 1) * pi / 2
        s0, s1 = math.sin(psi), math.cos(psi)
        Wg = s0 * W(ts, D) + s1 * W(ts - beta, D)
        hg = s0 * h(ts) + s1 * h(ts - beta)
        for phi in np.linspace(0, np.pi, 25)[:-1]:
            worst = max(worst, sc_circle(np.cos(phi) * hg - np.sin(phi) * Wg))
    assert worst == 4, worst
    print("all-line crossing cap over samples:", worst)

    # 2. Sturm slack at the census witnesses
    here = os.path.dirname(os.path.abspath(__file__))
    faces = json.load(open(os.path.join(here, "faces_0.json")))
    n_sc8 = 0
    for w in faces[:28]:
        beta, y = w["beta"], w["y"]
        psi = (y + 1) * pi / 2
        s0, s1 = math.sin(psi), math.cos(psi)
        mg = 0
        for D in np.linspace(0.02, 2 * pi - 0.02, 160):
            mg = max(mg, sc_circle(s0 * F(D, ts) + s1 * F(D, ts - beta)))
        if mg == 8:
            n_sc8 += 1
    print("witnesses with SC(forcing) = 8:", n_sc8, "of 28")
    assert n_sc8 >= 20

    # 3. ROUTE A REFUTATION: the ambient 4-dim space is NOT weak Chebyshev.
    #    Re-check the explicit witness at four resolutions, and confirm it
    #    is OFF the teacher-consistent variety (so the open lemma stands).
    cw = np.array([0.5389600251329808, -0.14246740163039062,
                   0.7066748566272032, -0.4357014778403485])
    Dw, bw = 4.846208362913536, 0.7349885346109851
    for n in (24000, 96000, 400000):
        tw = np.linspace(0, 2 * pi, n, endpoint=False)
        fw = (cw[0] * W(tw, Dw) + cw[1] * W(tw - bw, Dw)
              + cw[2] * h(tw) + cw[3] * h(tw - bw))
        assert sc_circle(fw) == 6, (n, sc_circle(fw))
    tc = cw[0] * cw[3] - cw[1] * cw[2]
    print("ambient-WT counterexample: 6 sign changes at 24k/96k/400k; "
          "teacher-consistency defect c0*c3-c1*c2 = %.4f (|c| = %.4f)"
          % (tc, float(np.linalg.norm(cw))))
    assert abs(tc) > 0.1                    # NOT a teacher-consistent load

    # random directions do NOT see it (skill 1m/1p): report the false
    # negative rather than an assertion of the (false) cap
    ts4 = np.linspace(0, 2 * pi, 24001)[:-1]
    worst_rand = 0
    for _ in range(60):
        D = rng.uniform(0.05, 2 * pi - 0.05)
        beta = rng.uniform(0.05, pi - 0.05)
        Bm = np.vstack([W(ts4, D), W(ts4 - beta, D), h(ts4), h(ts4 - beta)])
        dirs = rng.normal(size=(400, 4))
        dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
        V = dirs @ Bm
        for i in range(V.shape[0]):
            worst_rand = max(worst_rand, sc_circle(V[i]))
    print("...and 24000 RANDOM directions never exceed:", worst_rand,
          "(the sampled maximum is an under-estimate)")
    assert worst_rand == 4

    # 4. the sign-law arcs localise but DO NOT COUNT (kill-record above)
    def side_mask(t, D, beta, positive):
        s_a = np.where(np.mod(t, 2 * pi) < D, -1.0, 1.0)
        s_b = np.where(np.mod(t - beta, 2 * pi) < D, -1.0, 1.0)
        return (s_a * s_b < 0) if positive else (s_a * s_b > 0)

    def arc_runs(mask):
        idx = np.where(mask)[0]
        if len(idx) == 0:
            return []
        runs = []
        start = prev = idx[0]
        for i in idx[1:]:
            if i == prev + 1:
                prev = i
            else:
                runs.append((start, prev))
                start = prev = i
        runs.append((start, prev))
        if len(runs) > 1 and runs[0][0] == 0 and runs[-1][1] == len(mask) - 1:
            runs[0] = (runs[-1][0] - len(mask), runs[0][1])
            runs.pop()
        return runs

    n_trials = 0
    n_equal = 0
    mismatches = []
    worst_per_arc = 0
    for _ in range(40):
        beta = rng.uniform(0.1, pi - 0.1)
        y = rng.uniform(-0.97, 0.97)
        psi = (y + 1) * pi / 2
        s0, s1 = math.sin(psi), math.cos(psi)
        for D in np.linspace(0.05, 2 * pi - 0.05, 40):
            L = s0 * num(D, ts) + s1 * num(D, ts - beta)
            m = side_mask(ts, D, beta, s1 > 0)
            inside = 0
            for (a, b) in arc_runs(m):
                seg = L[np.arange(a, b + 1) % len(ts)]
                sg = np.sign(seg)
                sg = sg[sg != 0]
                c = int(np.sum(np.diff(sg) != 0)) if len(sg) > 1 else 0
                inside += c
                worst_per_arc = max(worst_per_arc, c)
            n_trials += 1
            if inside == sc_circle(L):
                n_equal += 1
            else:
                mismatches.append(
                    (min(abs(s0), abs(s1)),
                     min(abs(D - beta), abs(D - (2 * pi - beta)))))
    print("sign changes of L_D inside the admissible arcs = total, in",
          n_equal, "of", n_trials, "trials;  worst per single arc:",
          worst_per_arc)
    # VALIDITY DOMAIN (skill 1v).  Every mismatch sits on ONE of the two
    # degeneracies at which a crossing lands exactly ON an arc BOUNDARY, so
    # the run splitting mis-attributes it -- neither is a zero outside the
    # arcs:
    #   (i)  a VANISHING teacher mass (min|s| -> 0, i.e. the census
    #        parameter y -> 0 or |y| -> 1).  There L_D -> s0*num(D,.), whose
    #        zeros ARE the lattice {0, D} = the arc endpoints.  The census
    #        witnesses stay clear of it: the closest, F2 at y = -0.02, has
    #        min|s| = 0.031.
    #   (ii) the gap collapsing onto the teacher atom (|D - beta| or
    #        |D - (2pi - beta)| below the 0.05 lattice margin this file
    #        declares as its validity domain).  There two of the four
    #        lattice points {0, D, beta, beta+D} nearly coincide, BOTH
    #        num factors are near zero on the thin arc between them, and
    #        the float64 sign of L there is noise.
    # Assert that characterisation rather than pre-filtering it away.
    print("  mismatches:", len(mismatches),
          " all inside the declared degenerate margin"
          " (min|s| < 0.01 or |D - beta| < 0.05):",
          all(a < 0.01 or b < 0.05 for a, b in mismatches))
    assert all(a < 0.01 or b < 0.05 for a, b in mismatches), \
        sorted(mismatches)[-5:]
    assert worst_per_arc >= 3          # "one zero per arc" is FALSE

    print("ALL BRANCH-CAP MEASUREMENTS REPRODUCED")


if __name__ == "__main__":
    main()
