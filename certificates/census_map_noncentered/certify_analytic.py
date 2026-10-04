r"""ANALYTIC CLOSURE of the noncentered coincident-stratum root counts.

Certifies, for EVERY beta in (0, pi) -- including the two outer intervals whose
adaptive box coverings never completed -- the interior root counts of the
torque and potential Wronskians, by exact piecewise inequalities on the sinc-
reduced forms instead of a box covering whose margin must beat the box width.
By the proved mirror symmetry (kernel parities, CERTIFICATE.md section 2.4)
the counts transfer to (pi, 2pi).

## The four smooth pieces and their reduced forms

For beta in (0, pi) the t-period [0, 2pi] splits at the kink lattice
{0, beta, pi, beta+pi} into four pieces.  With p, q the half-branch
coordinates of t and t - beta (noncentered.py), Wtau = |pq| fhat and
Wpot = |pq| ghat where fhat = sin d + d sinc(p) sinc(q),
ghat = -sin d + d sinc(p) sinc(q), d in {beta, beta - 2pi}, sin d = sin beta
(validate.py section (c), sympy-exact on all four sign quadrants).  Interior
roots are exactly the interior zeros of fhat / ghat.

  L  (t in (0, beta)):        d = beta - 2pi;  sinc p sinc q =
        S(t) = sin t sin(beta-t) / ((pi-t)(pi-beta+t)) >= 0
  A1 (t in (beta, pi)):       d = beta;  sinc p sinc q = sinc(x) sinc(x+beta),
        x = pi - t in (0, pi-beta)
  A2 (t in (pi, beta+pi)):    d = beta;  sinc p sinc q = sinc(a) sinc(b),
        a = t - pi, b = beta - a, a + b = beta
  A3 (t in (beta+pi, 2pi)):   d = beta;  sinc p sinc q = sinc(x) sinc(x+beta),
        x = t - beta - pi in (0, pi-beta)

## The exact inequalities

(S1)  sin x - x cos x > 0 on (0, pi]  [proved below; equivalently sinc is
      strictly decreasing on (0, pi]].
(S2)  sinc(a) sinc(b) > sinc(a+b) for a, b > 0, a + b <= pi  [the identity
      (a+b) sin a sin b - ab sin(a+b)
        = a sin a (sin b - b cos b) + b sin b (sin a - a cos a),
      sympy-exact, plus (S1)].

Wtau: on A1, A2, A3 the sinc product is >= 0 (all arguments in [-pi, pi]), so
      fhat = sin beta + beta * (sinc p sinc q) >= sin beta > 0: NO interior
      roots.  All interior torque roots live on L.
Wpot: on L,  ghat = -sin beta - (2pi - beta) S <= -sin beta < 0.
      on A1/A3,  sinc(x) sinc(x+beta) < sinc(beta) for x > 0 [(S1): the second
      factor < sinc(beta), the first <= 1], so
      ghat = -sin beta + beta sinc(x)sinc(x+beta) < -sin beta + sin beta = 0.
      on A2,  sinc(a) sinc(b) > sinc(beta) by (S2), so ghat > 0.
      NO interior roots anywhere: Wpot count = 2 (the lattice roots), for
      every beta in (0, pi).

Wtau on L: with u = beta/2, v = pi - u, t = u + sigma,
      S = (sin^2 u - sin^2 sigma) / (v^2 - sigma^2) =: h(sigma), even in sigma,
      and fhat = 0  <=>  F(sigma) := (2pi - beta) h(sigma) - sin beta = 0.
(Q)   h is STRICTLY decreasing in |sigma| on (0, u): the derivative numerator
      is -[sin(2 sigma)(v^2 - sigma^2) - 2 sigma (sin^2 u - sin^2 sigma)], and
        Q(u, sigma) := sinc(2 sigma)(v^2 - sigma^2) - sin(u+sigma) sin(u-sigma)
      is > 0 on {0 <= sigma <= u <= pi/2} \ {corner (pi/2, pi/2)}:
      * sigma <= SIGMA0 = 1.171: certified by adaptive interval arithmetic
        below (margins are O(1) there);
      * sigma > SIGMA0 (the corner region; y = pi/2 - sigma <= 0.4,
        x = pi/2 - u <= y): the analytic chain
          sin 2y >= 2y (1 - (2/3) y^2),   pi - 2y <= pi,
          (x+y)(pi+x-y) >= (x+y)(pi - 0.4),
          sin(x+y) sin(y-x) <= (x+y) y,
        gives Q >= (x+y) y [2 (1 - (2/3) 0.16)(1 - 0.4/pi) - 1] >= 0.55 (x+y) y,
        strictly positive off the corner.
(F)   F(0) at the piece midpoint: F(0) = (2 sin u / v)(sin u - v cos u)
      [sympy-exact], and g(u) = sin u - (pi - u) cos u is strictly increasing
      on [0, pi/2] (g' = 2 cos u + (pi - u) sin u > 0, both terms >= 0 and not
      simultaneously 0), g(0) = -pi < 0, with unique zero u* = beta*_1 / 2
      (tan u* = pi - u*, certified enclosure in certify_beta_star).  At the
      endpoints F -> -sin beta < 0.  Hence on L:
        beta < beta*_1 :  F < 0 everywhere,           0 interior roots;
        beta = beta*_1 :  F(0) = 0, F < 0 elsewhere,  1 interior (double) root;
        beta > beta*_1 :  F(0) > 0, one strict crossing per side,
                                                       2 interior roots.

## Verdict

For every beta in (0, pi):   Wpot_total = 2, and
    Wtau_total = 2 on (0, beta*_1), 3 at beta*_1, 4 on (beta*_1, pi).
Mirrored by the exact symmetry to (pi, 2pi).  The remaining columns beta = 0
(all determinants identically zero) and beta = pi (sin beta = 0: all three
reduce to the common factor d * sinc p sinc q, zero set exactly the lattice
{0, pi}) are the exact degenerate cases of CERTIFICATE.md sections 3.1/6.

This SUPERSEDES the incomplete adaptive coverings of the outer intervals and
closes the six beta collars: no box covering is needed anywhere on the
coincident stratum.

Writes analytic_counts.json.
"""
import json
import os
import time

import mpmath
import sympy as sp
from mpmath import iv

import noncentered as J
from certify_coincident import sinc_iv

J.set_prec(120)
PIv = J.PI_IV()
OUT = {"object": "analytic closure of the coincident-stratum root counts",
       "checks": {}}

t0 = time.time()

# ---------------------------------------------------------------- sympy
u_, s_, a_, b_, t_, be_ = sp.symbols("u sigma a b t beta", positive=True)
v_ = sp.pi - u_

# (1) the L-piece S in the (u, sigma) chart
S_t = sp.sin(t_) * sp.sin(be_ - t_) / ((sp.pi - t_) * (sp.pi - be_ + t_))
S_us = (sp.sin(u_) ** 2 - sp.sin(s_) ** 2) / (v_ ** 2 - s_ ** 2)
ok1 = sp.simplify(sp.expand_trig(
    S_t.subs({t_: u_ + s_, be_: 2 * u_}) - S_us)) == 0

# (2) derivative numerator of h(sigma)
num = sp.simplify(sp.together(sp.diff(S_us, s_)) * (v_ ** 2 - s_ ** 2) ** 2)
target = -(sp.sin(2 * s_) * (v_ ** 2 - s_ ** 2)
           - 2 * s_ * (sp.sin(u_) ** 2 - sp.sin(s_) ** 2))
ok2 = sp.simplify(sp.expand_trig(sp.expand(num - target))) == 0

# (3) sin(u+s) sin(u-s) = sin^2 u - sin^2 s
ok3 = sp.simplify(sp.expand_trig(
    sp.sin(u_ + s_) * sp.sin(u_ - s_) - (sp.sin(u_) ** 2 - sp.sin(s_) ** 2))) == 0

# (4) F(0) = (2 sin u / v)(sin u - v cos u), from F = (2pi-2u) h(0) - sin 2u
F0 = (2 * sp.pi - 2 * u_) * (sp.sin(u_) ** 2 / v_ ** 2) - sp.sin(2 * u_)
ok4 = sp.simplify(F0 - (2 * sp.sin(u_) / v_) * (sp.sin(u_) - v_ * sp.cos(u_))) == 0

# (5) the superadditivity identity behind (S2)
lhs = (a_ + b_) * sp.sin(a_) * sp.sin(b_) - a_ * b_ * sp.sin(a_ + b_)
rhs = a_ * sp.sin(a_) * (sp.sin(b_) - b_ * sp.cos(b_)) \
    + b_ * sp.sin(b_) * (sp.sin(a_) - a_ * sp.cos(a_))
ok5 = sp.simplify(sp.expand_trig(sp.expand(lhs - rhs))) == 0

# (6) the L-piece fhat/ghat forms against the quadrant reduction:
#     E * F(p, q) with p = pi - t, q = beta - t - pi, E = -1 equals
#     |pq| (sin beta - (2pi - beta) S(t)),  and the ghat analogue.
p_ = sp.pi - t_
q_ = be_ - t_ - sp.pi
Fpq = p_ * q_ * sp.sin(q_ - p_) + (q_ - p_) * sp.sin(p_) * sp.sin(q_)
Gpq = -p_ * q_ * sp.sin(q_ - p_) + (q_ - p_) * sp.sin(p_) * sp.sin(q_)
absPQ = p_ * (sp.pi + t_ - be_)          # |p q| = (pi - t)(pi - (beta - t))
ok6a = sp.simplify(sp.expand_trig(sp.expand(
    -Fpq - absPQ * (sp.sin(be_) - (2 * sp.pi - be_) * S_t)))) == 0
ok6b = sp.simplify(sp.expand_trig(sp.expand(
    -Gpq - absPQ * (-sp.sin(be_) - (2 * sp.pi - be_) * S_t)))) == 0

OUT["checks"]["sympy"] = {
    "S_in_(u,sigma)_chart": bool(ok1),
    "h_derivative_numerator": bool(ok2),
    "product_formula": bool(ok3),
    "F_at_midpoint": bool(ok4),
    "superadditivity_identity": bool(ok5),
    "L_piece_fhat_form": bool(ok6a),
    "L_piece_ghat_form": bool(ok6b),
}
assert all(OUT["checks"]["sympy"].values()), OUT["checks"]["sympy"]
print("sympy identities:", OUT["checks"]["sympy"])

# ------------------------------------------------- (S1) sin x - x cos x > 0
# on (0, 2.8]: sin x - x cos x >= x^3/3 - x^5/24 = x^3 (1/3 - x^2/24) > 0 for
# x^2 < 8 (series bounds sin x >= x - x^3/6 and cos x <= 1 - x^2/2 + x^4/24,
# both alternating-series enclosures, valid for x >= 0).
# on [2.8, pi]: cos x <= cos(2.8) < 0 (certified below), sin x >= 0, so
# sin x - x cos x >= -x cos(2.8) >= 2.8 |cos(2.8)| > 0.
c28 = iv.cos(J.I("2.8"))
OUT["checks"]["S1_cos28_negative"] = J.hi(c28) < 0
assert J.hi(c28) < 0

# ------------------------------------------------- (Q) interval certification
SIGMA0 = mpmath.mpf("1.171")
# corner-region constants:  pi/2 - 0.4 <= SIGMA0  and the margin constant
corner_gap = J.hi(PIv / 2 - J.I("0.4"))
OUT["checks"]["corner_covers"] = float(corner_gap) <= float(SIGMA0)
marg = 2 * (1 - sp.Rational(2, 3) * sp.Rational(4, 10) ** 2) \
    * (1 - sp.Rational(4, 10) / sp.pi) - 1
OUT["checks"]["corner_margin_constant"] = float(marg.evalf(30))
assert OUT["checks"]["corner_covers"] and marg.evalf(30) > 0.5


def Qval(U, Ssig):
    V = PIv - U
    return sinc_iv(2 * Ssig) * (V * V - Ssig * Ssig) \
        - (iv.sin(U) ** 2 - iv.sin(Ssig) ** 2)


half = J.lo(PIv / 2)
todo = [(mpmath.mpf(0), min(SIGMA0, half), mpmath.mpf(0), half)]
boxes = 0
min_margin = None
while todo:
    s0, s1, u0, u1 = todo.pop()
    if u1 <= s0:                       # entirely below the diagonal sigma <= u
        continue
    Q = Qval(iv.mpf([u0, u1]), iv.mpf([s0, s1]))
    boxes += 1
    lo = J.lo(Q)
    if lo > 0:
        if min_margin is None or lo < min_margin:
            min_margin = lo
        continue
    if (s1 - s0) < mpmath.mpf("1e-6") and (u1 - u0) < mpmath.mpf("1e-6"):
        raise SystemExit("Q certification FAILED at sigma=[%s,%s] u=[%s,%s]"
                         % (s0, s1, u0, u1))
    if s1 - s0 >= u1 - u0:
        m = (s0 + s1) / 2
        todo += [(s0, m, u0, u1), (m, s1, u0, u1)]
    else:
        m = (u0 + u1) / 2
        todo += [(s0, s1, u0, m), (s0, s1, m, u1)]
OUT["checks"]["Q_interval"] = {
    "region": "0 <= sigma <= min(u, %s), 0 <= u <= pi/2" % SIGMA0,
    "boxes": boxes,
    "min_certified_lower_bound": float(min_margin),
}
print("Q > 0 certified on the interval region: %d boxes, min margin %.3e"
      % (boxes, float(min_margin)))

# ------------------------------------------------- (F) g' > 0 on [0, pi/2]
# g'(u) = 2 cos u + (pi - u) sin u: both terms are >= 0 on [0, pi/2]
# (cos u >= 0, sin u >= 0, pi - u > 0), and they never vanish together:
# for u in [0, 1.4] the first term is >= 2 cos(1.4) > 0; for u in [1.4, pi/2]
# the second is >= (pi/2) sin(1.4) > 0.
c14 = iv.cos(J.I("1.4"))
s14 = iv.sin(J.I("1.4"))
OUT["checks"]["gprime"] = {
    "cos_1.4_positive": J.lo(c14) > 0,
    "pi_half_sin_1.4_positive": J.lo(PIv / 2 * s14) > 0,
}
assert J.lo(c14) > 0 and J.lo(PIv / 2 * s14) > 0

# ------------------------------------------------- cross-check vs coverings
# analytic counts vs the committed 309/305-box covering of (beta*_1, pi) and
# (pi, beta*_2): every covered box there certified (Wtau, Wpot) = (4, 2), and
# the analytic count on (beta*_1, pi) is (4, 2).  Spot-check the float sign
# structure at a few betas as an instrument sanity check (float only).
import math
sanity = []
for b in (0.05, 0.5, 1.0, 2.0, 2.3, 2.8, 3.1):
    n = 0
    for k in range(1, 40000):
        t1, t2 = b * k / 40000, b * (k + 1) / 40000
        if J.f_Wtau(b, t1) * J.f_Wtau(b, t2) < 0:
            n += 1
    pred = 0 if b < 2.2256696309587180 else 2
    sanity.append({"beta": b, "float_interior_sign_changes_on_L": n,
                   "analytic": pred, "agree": n == pred})
OUT["checks"]["float_sanity"] = sanity
assert all(z["agree"] for z in sanity)

OUT["verdict"] = {
    "Wpot_total": "2 for every beta in (0, 2pi) \\ {pi}",
    "Wtau_total": "2 on (0, beta*_1) u (beta*_2, 2pi); 3 at beta*_1, beta*_2 "
                  "(one interior double root); 4 on (beta*_1, pi) u (pi, beta*_2)",
    "scope": "every beta; the columns beta = 0 and beta = pi are the exact "
             "degenerate cases (all determinants identically zero at 0; "
             "common zero set = lattice {0, pi} at pi)",
    "supersedes": "the adaptive beta coverings (309/305 boxes) and the "
                  "incomplete outer-interval runs; no beta collar remains "
                  "on the coincident stratum",
}
OUT["seconds"] = time.time() - t0
HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, "analytic_counts.json"), "w") as fh:
    json.dump(OUT, fh, indent=1)
print("verdict:", json.dumps(OUT["verdict"], indent=1))
print("%.0fs; wrote analytic_counts.json" % OUT["seconds"])
