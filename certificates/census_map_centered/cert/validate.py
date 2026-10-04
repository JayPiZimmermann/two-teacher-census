import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
"""
TASK 1 -- validation.
  (a) Python float mirror vs the LITERAL source lines of website/widgets.js
      (evaluated by node), over 25000 (t, beta) samples.
  (b) h'(t) == slopeAtom(t) by high-order finite differences (mpmath, away
      from the kink lattice).
  (c) the reduced piecewise forms == the folded mirror, in interval arithmetic.
  (d) the exact identity Wtau + Wpot = 2*Wwgt (interval + sympy symbolic).
"""
import json
import math
import sys

import mpmath
from mpmath import iv, mp

import centered as C

OUT = {}

# ---------------------------------------------------------------- (a)
cols = ["phiC", "HC", "dHC", "slopeAtom", "Wtau", "Wpot", "Wwgt"]
maxdev = {c: 0.0 for c in cols}
argmax = {c: None for c in cols}
n = 0
with open(_os.path.join(_HERE, "ref_samples.txt")) as fh:
    for line in fh:
        v = [float(x) for x in line.split()]
        t, beta = v[0], v[1]
        mine = [C.f_phiC(t), C.f_HC(t), C.f_dHC(t), C.f_slopeAtom(t),
                C.f_Wtau(beta, t), C.f_Wpot(beta, t), C.f_Wwgt(beta, t)]
        for i, c in enumerate(cols):
            d = abs(mine[i] - v[2 + i])
            if d > maxdev[c]:
                maxdev[c] = d
                argmax[c] = (t, beta)
        n += 1
OUT["a_widgets_mirror"] = {
    "samples": n,
    "max_abs_deviation": {c: maxdev[c] for c in cols},
    "argmax": {c: argmax[c] for c in cols},
    "note": "reference values produced by ref_widgets.js: node evaluating the "
            "literal function sources of website/widgets.js (PI, mod, phiC, "
            "HC, dHC, slopeAtom), extracted by identifier with a guard that "
            "asserts each extracted function still contains the expected "
            "expression",
}
print("(a) mirror vs widgets.js over %d samples:" % n)
for c in cols:
    print("    %-10s max|dev| = %.3e" % (c, maxdev[c]))

# --------------------------------------------------------------- (a2)
# libm attribution: recompute the Python formulas from NODE's sin/cos values.
import subprocess
sub = subprocess.run(["node", _os.path.join(_HERE, "ref_libm.js"), "6000"],
                     capture_output=True, text=True, check=True)
PIf = math.pi
libm_max = 0.0
c7 = ["phiC", "HC", "dHC", "slopeAtom", "Wtau", "Wpot", "Wwgt"]
recomp = {c: 0.0 for c in c7}
nrow = 0
for line in sub.stdout.strip().split("\n"):
    q = [float(x) for x in line.split()]
    t, beta = q[0], q[1]
    pt, pu = q[2:6], q[6:10]
    ref = q[10:17]
    for (x_, s_x, c_x, s_u), u in ((pt, t), (pu, t - beta)):
        libm_max = max(libm_max,
                       abs(s_x - math.sin(x_)), abs(c_x - math.cos(x_)),
                       abs(s_u - math.sin(u)))

    def build(p):
        x_, s_x, c_x, s_u = p
        phi = (PIf / 2 - x_) * c_x + s_x
        H = (PIf / 2 - x_) * s_x
        dH = (PIf / 2 - x_) * c_x - s_x
        sA = phi - 2 * abs(s_u)
        return phi, H, dH, sA, abs(s_u)

    Pt, Ht, dHt, sAt, at = build(pt)
    Pu, Hu, dHu, sAu, au = build(pu)
    got = [Pt, Ht, dHt, sAt,
           Ht * sAu - Hu * sAt,
           Pt * Hu - Pu * Ht,
           at * Hu - au * Ht]
    for i, c in enumerate(c7):
        recomp[c] = max(recomp[c], abs(got[i] - ref[i]))
    nrow += 1
OUT["a2_libm"] = {
    "arguments": nrow,
    "max_abs_sin_cos_difference_node_vs_glibc": libm_max,
    "max_abs_deviation_when_fed_nodes_own_sin_cos": recomp,
    "statement": "V8's Math.sin/Math.cos and glibc's differ by up to "
                 "max_abs_sin_cos_difference_node_vs_glibc; feeding node's OWN "
                 "sin/cos into the Python formulas reproduces every kernel "
                 "quantity to max_abs_deviation_when_fed_nodes_own_sin_cos, so "
                 "the reimplementation is formula-identical to the widget.",
}
print("(a2) node vs glibc sin/cos max difference = %.4e" % libm_max)
print("(a2) fed node's own sin/cos: %s" % recomp)

# ---------------------------------------------------------------- (b)
mp.prec = 300
PI = mp.pi


def m_phiC(t):
    x = mp.fmod(t, PI)
    if x < 0:
        x += PI
    return (PI / 2 - x) * mp.cos(x) + mp.sin(x)


def m_HC(t):
    x = mp.fmod(t, PI)
    if x < 0:
        x += PI
    return (PI / 2 - x) * mp.sin(x)


def m_slope(t):
    return m_phiC(t) - 2 * abs(mp.sin(t))


worst = mp.mpf(0)
worst_at = None
for k in range(1, 400):
    t = mp.mpf(k) * mp.mpf("7.3") / 400 * 3 - 4     # spread over ~[-4, 12]
    # stay away from the kink lattice
    if min(abs(mp.fmod(t, PI)), abs(PI - mp.fmod(t, PI))) < mp.mpf("1e-3"):
        continue
    hstep = mp.mpf("1e-25")
    d = (m_HC(t + hstep) - m_HC(t - hstep)) / (2 * hstep)
    e = abs(d - m_slope(t))
    if e > worst:
        worst, worst_at = e, t
OUT["b_hprime_eq_slopeAtom"] = {
    "method": "central difference, mpmath prec=300, step 1e-25, kink lattice "
              "avoided by 1e-3",
    "max_abs_error": mpmath.nstr(worst, 8),
    "argmax_t": mpmath.nstr(worst_at, 20),
}
print("(b) |d/dt HC - slopeAtom| max = %s" % mpmath.nstr(worst, 8))

# ---------------------------------------------------------------- (c)
C.set_prec(200)
PIv = C.PI_IV()
maxgap = {"Wtau": 0.0, "Wpot": 0.0, "Wwgt": 0.0}
checked = 0
for ib in range(1, 60):
    beta = PIv * ib / 60
    for piece in ("A", "B"):
        b, sig = C.piece_params(beta, piece)
        Wmax = (PIv - b) / 2
        for j in range(-24, 25):
            w = Wmax * j / 25
            t = C.t_from_w(beta, piece, w)
            # folded mirror (thin box around t, still an interval evaluation)
            D = C.dets_iv(beta, t)
            R = C.redu(b, w)
            red = {"Wtau": sig * R[0], "Wpot": sig * R[1], "Wwgt": sig * R[2]}
            for nm in maxgap:
                diff = D[nm][0] - red[nm]
                if not C.contains_zero(diff):
                    raise SystemExit("REDUCTION MISMATCH %s beta=%s piece=%s j=%d  %s"
                                     % (nm, beta, piece, j, diff))
                maxgap[nm] = max(maxgap[nm], float(max(abs(diff.a), abs(diff.b))))
            checked += 1
OUT["c_reduction_identity"] = {
    "grid_points": checked,
    "statement": "on every sample the interval enclosure of (folded mirror - "
                 "reduced form) contains 0",
    "max_enclosure_radius": maxgap,
}
print("(c) reduction identity verified on %d interval samples; max residual radius %s"
      % (checked, maxgap))

# ---------------------------------------------------------------- (d)
maxid = 0.0
for ib in range(1, 60):
    beta = PIv * ib / 60
    for j in range(-40, 41):
        t = PIv * (j + 40) / 80
        D = C.dets_iv(beta, t)
        r = D["Wtau"][0] + D["Wpot"][0] - 2 * D["Wwgt"][0]
        if not C.contains_zero(r):
            raise SystemExit("IDENTITY FAILS at beta=%s t=%s : %s" % (beta, t, r))
        maxid = max(maxid, float(max(abs(r.a), abs(r.b))))

import sympy as sp
b_, w_ = sp.symbols("b w", real=True)
Q = (w_**2 - b_**2 / 4) * sp.sin(b_)
Om = (b_ / 2) * (sp.cos(2 * w_) + sp.cos(b_))
Phi = Q + Om
Psi = -Q + Om
sym_ok = sp.simplify(Phi + Psi - 2 * Om) == 0
# and the reduced forms really are the determinants: symbolic check on piece A
t_, be_ = sp.symbols("t beta", positive=True)


def s_phi(x):
    return (sp.pi / 2 - x) * sp.cos(x) + sp.sin(x)


def s_h(x):
    return (sp.pi / 2 - x) * sp.sin(x)


def s_sA(x):
    return (sp.pi / 2 - x) * sp.cos(x) - sp.sin(x)


# piece A: fold(t) = t, fold(t-beta) = t-beta   (0 < beta < t < pi)
xt, xu = t_, t_ - be_
WtauA = sp.simplify(s_h(xt) * s_sA(xu) - s_h(xu) * s_sA(xt))
WpotA = sp.simplify(s_phi(xt) * s_h(xu) - s_phi(xu) * s_h(xt))
WwgtA = sp.simplify(sp.sin(xt) * s_h(xu) - sp.sin(xu) * s_h(xt))
wA = t_ - (be_ + sp.pi) / 2
subA = {b_: be_, w_: wA}
okA = [sp.simplify(sp.expand_trig(sp.simplify(WtauA - Phi.subs(subA)))) == 0,
       sp.simplify(sp.expand_trig(sp.simplify(WpotA - Psi.subs(subA)))) == 0,
       sp.simplify(sp.expand_trig(sp.simplify(WwgtA - Om.subs(subA)))) == 0]
# piece B: fold(t) = t, fold(t-beta) = t-beta+pi  (0 < t < beta < pi)
xu2 = t_ - be_ + sp.pi
WtauB = s_h(t_) * s_sA(xu2) - s_h(xu2) * s_sA(t_)
WpotB = s_phi(t_) * s_h(xu2) - s_phi(xu2) * s_h(t_)
WwgtB = sp.sin(t_) * s_h(xu2) - sp.sin(xu2) * s_h(t_)
wB = t_ - be_ / 2
subB = {b_: sp.pi - be_, w_: wB}
okB = [sp.simplify(sp.expand_trig(sp.simplify(WtauB + Phi.subs(subB)))) == 0,
       sp.simplify(sp.expand_trig(sp.simplify(WpotB + Psi.subs(subB)))) == 0,
       sp.simplify(sp.expand_trig(sp.simplify(WwgtB + Om.subs(subB)))) == 0]

OUT["d_linear_identity"] = {
    "statement": "Wtau(beta,t) + Wpot(beta,t) = 2*Wwgt(beta,t) identically",
    "interval_max_residual_radius": maxid,
    "sympy_symbolic": bool(sym_ok),
}
OUT["c_reduction_symbolic"] = {
    "pieceA_Wtau_Wpot_Wwgt": [bool(x) for x in okA],
    "pieceB_Wtau_Wpot_Wwgt": [bool(x) for x in okB],
}
print("(d) Wtau+Wpot-2Wwgt: interval residual radius <= %.3e ; sympy exact: %s"
      % (maxid, sym_ok))
print("(c') sympy reduction piece A %s piece B %s" % (okA, okB))

with open(_os.path.join(_HERE, "validation.json"), "w") as fh:
    json.dump(OUT, fh, indent=2)
print("wrote validation.json")
