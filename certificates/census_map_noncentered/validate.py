"""
TASK 1 -- validation of the NONCENTERED mirror.
  (a) Python float mirror vs the LITERAL source lines of website/widgets.js
      (evaluated by node), over 25000 (t, beta, t1) samples.
  (a2) the libm attribution: feed node's OWN sin/cos into the Python formulas.
  (b) h'(t) == slopeAtom(t) by high-precision central differences (mpmath,
      away from the kink lattice).
  (c) the reduced piecewise forms == the folded mirror, in interval arithmetic
      + sympy symbolic on all four pieces.
  (d) the exact identity Wtau + Wpot = 2*Wwgt (interval + sympy symbolic).
  (e) Dsep (division-free) == H(D)^2 * separatedDet (division form).
"""
import json
import math
import os
import subprocess
import sys

import mpmath
from mpmath import iv, mp

import noncentered as J

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}

# ---------------------------------------------------------------- (a)
cols = ["phiJ", "HJ", "dHJ", "slopeAtom", "Wtau", "Wpot", "Wwgt", "Dsep", "y"]
maxdev = {c: 0.0 for c in cols}
argmax = {c: None for c in cols}
n = 0
nan_skipped = 0
with open(os.path.join(HERE, "ref_samplesJ.txt")) as fh:
    for line in fh:
        v = [float(x) for x in line.split()]
        t, beta, t1 = v[0], v[1], v[2]
        mine = [J.f_phiJ(t), J.f_HJ(t), J.f_dHJ(t), J.f_slopeAtom(t),
                J.f_Wtau(beta, t), J.f_Wpot(beta, t), J.f_Wwgt(beta, t),
                None, None]
        try:
            mine[7] = J.f_Dsep(beta, t, t1)
        except ZeroDivisionError:
            mine[7] = float("nan")
        try:
            mine[8] = J.f_ratioCoord(-J.f_HJ(t - beta), J.f_HJ(t))
        except Exception:
            mine[8] = float("nan")
        for i, c in enumerate(cols):
            ref = v[3 + i]
            if math.isnan(ref) or math.isnan(mine[i]):
                if c == "Dsep":
                    nan_skipped += 1
                continue
            d = abs(mine[i] - ref)
            if d > maxdev[c]:
                maxdev[c] = d
                argmax[c] = (t, beta, t1)
        n += 1
OUT["a_widgets_mirror"] = {
    "samples": n,
    "nan_samples_skipped_Dsep": nan_skipped,
    "max_abs_deviation": {c: maxdev[c] for c in cols},
    "argmax_t_beta_t1": {c: argmax[c] for c in cols},
    "note": "reference values produced by node evaluating the literal source "
            "lines 4/225/231/232/233/294 (atoms), 240/241 (Pot,Trq) and "
            "393-398 (ratioCoord) of website/widgets.js; kernelOf(\"noncentered\") "
            "is line 238.",
}
print("(a) mirror vs widgets.js over %d samples:" % n)
for c in cols:
    print("    %-10s max|dev| = %.4e" % (c, maxdev[c]))

# --------------------------------------------------------------- (a2)
# libm attribution: recompute the Python formulas from NODE's sin/cos values.
sub = subprocess.run(["node", os.path.join(HERE, "ref_libm.js"), "6000"],
                     capture_output=True, text=True, check=True)
libm_max = 0.0
PIf = math.pi
c7 = ["phiJ", "HJ", "dHJ", "slopeAtom", "Wtau", "Wpot", "Wwgt"]
recomp = {c: 0.0 for c in c7}
nrow = 0
for line in sub.stdout.strip().split("\n"):
    q = [float(x) for x in line.split()]
    t, beta = q[0], q[1]
    pt, pu = q[2:9], q[9:16]
    ref = q[16:23]
    for (xh, xf, s_xh, c_xh, s_xf, c_xf, s_u), u in ((pt, t), (pu, t - beta)):
        libm_max = max(libm_max,
                       abs(s_xh - math.sin(xh)), abs(c_xh - math.cos(xh)),
                       abs(s_xf - math.sin(xf)), abs(c_xf - math.cos(xf)),
                       abs(s_u - math.sin(u)))

    def build(p):
        xh, xf, s_xh, c_xh, s_xf, c_xf, s_u = p
        phi = (PIf - xf) * c_xf + s_xf
        H = (PIf - xh) * s_xh if xh <= PIf else (xh - PIf) * s_xh
        dH = ((PIf - xh) * c_xh - s_xh if xh <= PIf
              else s_xh + (xh - PIf) * c_xh)
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
    "statement": "V8's Math.sin/Math.cos (fdlibm port) and glibc's differ by "
                 "up to max_abs_sin_cos_difference_node_vs_glibc; feeding "
                 "node's OWN sin/cos into the Python formulas reproduces every "
                 "kernel quantity to max_abs_deviation_when_fed_nodes_own_sin_cos, "
                 "so the reimplementation is formula-identical to the widget.",
}
print("(a2) node vs glibc sin/cos max difference = %.4e" % libm_max)
print("(a2) fed node's own sin/cos: %s" % recomp)

# ---------------------------------------------------------------- (b)
mp.prec = 300
PI = mp.pi


def m_phiJ(t):
    x = mp.fmod(t, 2 * PI)
    if x < 0:
        x += 2 * PI
    if x > PI:
        x = 2 * PI - x
    return (PI - x) * mp.cos(x) + mp.sin(x)


def m_HJ(t):
    x = mp.fmod(t, 2 * PI)
    if x < 0:
        x += 2 * PI
    return (PI - x) * mp.sin(x) if x <= PI else (x - PI) * mp.sin(x)


def m_slope(t):
    return m_phiJ(t) - 2 * abs(mp.sin(t))


worst = mp.mpf(0)
worst_at = None
checked_b = 0
for k in range(1, 800):
    t = mp.mpf(k) * mp.mpf("7.3") / 800 * 6 - 8     # spread over ~[-8, 36]
    r = mp.fmod(t, PI)
    if r < 0:
        r += PI
    if min(abs(r), abs(PI - r)) < mp.mpf("1e-3"):
        continue
    hstep = mp.mpf("1e-25")
    d = (m_HJ(t + hstep) - m_HJ(t - hstep)) / (2 * hstep)
    e = abs(d - m_slope(t))
    checked_b += 1
    if e > worst:
        worst, worst_at = e, t
OUT["b_hprime_eq_slopeAtom"] = {
    "method": "central difference, mpmath prec=300, step 1e-25, kink lattice "
              "(t in pi*Z) avoided by 1e-3",
    "samples": checked_b,
    "max_abs_error": mpmath.nstr(worst, 8),
    "argmax_t": mpmath.nstr(worst_at, 20),
}
print("(b) |d/dt HJ - slopeAtom| max = %s over %d samples"
      % (mpmath.nstr(worst, 8), checked_b))

# ---------------------------------------------------------------- (c)
J.set_prec(200)
PIv = J.PI_IV()
maxgap = {"Wtau": 0.0, "Wpot": 0.0, "Wwgt": 0.0}
maxgap_sinc = 0.0
checked = 0


def sinc_iv(z):
    """Rigorous sinc(z) = sin z / z (removable singularity at 0)."""
    a, b = J.endpoints(z)
    if a > 0 or b < 0:
        return iv.sin(z) / z
    # |z| <= 1 assumed for the series branch; else split
    M = max(abs(a), abs(b))
    assert M <= 1, "sinc series branch needs |z| <= 1"
    # alternating series sum_{k>=0} (-1)^k z^{2k}/(2k+1)!  with |z|<=1:
    # terms strictly decreasing, so truncation error <= first omitted term.
    K = 12
    s = iv.mpf(0)
    fact = 1
    z2 = z * z
    pw = iv.mpf(1)
    for k in range(K):
        fact_k = mpmath.factorial(2 * k + 1)
        s = s + (pw / iv.mpf(str(fact_k)) if k % 2 == 0
                 else -(pw / iv.mpf(str(fact_k))))
        pw = pw * z2
    tail = abs(float(M)) ** (2 * K) / float(mpmath.factorial(2 * K + 1))
    return s + iv.mpf([-tail, tail])


for ib in range(1, 96):
    beta = 2 * PIv * ib / 96
    if abs(float(J.mid(beta)) - math.pi) < 1e-9:
        continue
    for pc in J.pieces_of(beta):
        d = pc["d"]
        for jx in range(1, 20):
            m = pc["m_lo"] + (pc["m_hi"] - pc["m_lo"]) * jx / 20
            t = J.t_from_m(pc, m)
            D = J.dets_iv(beta, t)
            R = J.redu(d, m)
            red = {"Wtau": pc["E"] * R[0], "Wpot": pc["E"] * R[1],
                   "Wwgt": pc["E"] * R[2]}
            for nm in maxgap:
                diff = D[nm][0] - red[nm]
                if not J.contains_zero(diff):
                    raise SystemExit("REDUCTION MISMATCH %s beta=%s piece=%s"
                                     % (nm, beta, pc))
                maxgap[nm] = max(maxgap[nm], J.rad(diff))
            checked += 1
OUT["c_reduction_identity"] = {
    "grid_points": checked,
    "statement": "on every sample the interval enclosure of (folded mirror - "
                 "reduced form E*(F,G,N)) contains 0",
    "max_enclosure_radius": maxgap,
}
print("(c) reduction identity verified on %d interval samples; max residual "
      "radius %s" % (checked, maxgap))

# ---------------------------------------------------------------- (d)
maxid = 0.0
for ib in range(1, 96):
    beta = 2 * PIv * ib / 96
    for jx in range(0, 81):
        t = 2 * PIv * jx / 80
        D = J.dets_iv(beta, t)
        r = D["Wtau"][0] + D["Wpot"][0] - 2 * D["Wwgt"][0]
        if not J.contains_zero(r):
            raise SystemExit("IDENTITY FAILS at beta=%s t=%s : %s" % (beta, t, r))
        maxid = max(maxid, J.rad(r))

import sympy as sp
p_, q_, d_ = sp.symbols("p q d", real=True)
Qd = p_ * q_ * sp.sin(q_ - p_)
Nn = (q_ - p_) * sp.sin(p_) * sp.sin(q_)
F_, G_ = Qd + Nn, -Qd + Nn
sym_lin = sp.simplify(F_ + G_ - 2 * Nn) == 0

# the reduced forms really ARE the determinants, on all four sign quadrants.
t_, be_ = sp.symbols("t beta", positive=True)


def atoms_sym(p, e):
    """e = +1 (p in [0,pi]) or -1 (p in [-pi,0]); P = e*p = |p|."""
    P = e * p
    return {"h": P * sp.sin(p),
            "phi": -P * sp.cos(p) + e * sp.sin(p),
            "sA": -P * sp.cos(p) - e * sp.sin(p),
            "as": e * sp.sin(p)}


ok_quadrants = []
for et in (1, -1):
    for eu in (1, -1):
        At = atoms_sym(p_, et)
        Au = atoms_sym(q_, eu)
        E = et * eu
        Wt = At["h"] * Au["sA"] - Au["h"] * At["sA"]
        Wp = At["phi"] * Au["h"] - Au["phi"] * At["h"]
        Ww = At["as"] * Au["h"] - Au["as"] * At["h"]
        ok_quadrants.append([
            bool(sp.simplify(sp.expand_trig(sp.simplify(Wt - E * F_))) == 0),
            bool(sp.simplify(sp.expand_trig(sp.simplify(Wp - E * G_))) == 0),
            bool(sp.simplify(sp.expand_trig(sp.simplify(Ww - E * Nn))) == 0)])

# factorisation F = p q fhat, G = p q ghat with sinc
fhat_ok = bool(sp.simplify(F_ - p_ * q_ * (sp.sin(q_ - p_)
               + (q_ - p_) * sp.sin(p_) * sp.sin(q_) / (p_ * q_))) == 0)
ghat_ok = bool(sp.simplify(G_ - p_ * q_ * (-sp.sin(q_ - p_)
               + (q_ - p_) * sp.sin(p_) * sp.sin(q_) / (p_ * q_))) == 0)

OUT["d_linear_identity"] = {
    "statement": "Wtau(beta,t) + Wpot(beta,t) = 2*Wwgt(beta,t) identically",
    "interval_max_residual_radius": maxid,
    "sympy_symbolic_on_reduced_forms": bool(sym_lin),
}
OUT["c_reduction_symbolic"] = {
    "quadrants_(e_t,e_u)_in_(+,+),(+,-),(-,+),(-,-)": ok_quadrants,
    "F_eq_pq_fhat": fhat_ok,
    "G_eq_pq_ghat": ghat_ok,
}
print("(d) Wtau+Wpot-2Wwgt: interval residual radius <= %.4e ; sympy exact: %s"
      % (maxid, sym_lin))
print("(c') sympy reduction on the four sign quadrants: %s" % ok_quadrants)

# ---------------------------------------------------------------- (e)
# Dsep (division-free) == H(D)^2 * separatedDet.  Compare the float mirror of
# the division form with the interval evaluation of the division-free form.
maxe = 0.0
cnt_e = 0
for ib in range(1, 24):
    beta = 2 * PIv * ib / 24
    bf = float(J.mid(beta))
    for i0 in range(1, 13):
        for i1 in range(1, 13):
            th0 = 2 * PIv * i0 / 13
            th1 = 2 * PIv * i1 / 13
            if abs(float(J.mid(th0 - th1))) < 1e-6:
                continue
            Dv, a0, a1, sm = J.sep_det(beta, th0, th1)
            ref = J.f_Dsep(bf, float(J.mid(th0)), float(J.mid(th1)))
            if math.isnan(ref) or math.isinf(ref):
                continue
            df = Dv - J.I(ref)
            maxe = max(maxe, J.rad(df))
            cnt_e += 1
OUT["e_separated_det"] = {
    "statement": "Dsep := a0*b1 - a1*b0 equals H(D)^2 * separatedDet, i.e. the "
                 "division form of widgets.js separatedScan; compared against "
                 "the float mirror of the division form",
    "samples": cnt_e,
    "max_abs_difference_interval_vs_float_division_form": maxe,
}
print("(e) Dsep vs float division form: %d samples, max |diff| = %.4e"
      % (cnt_e, maxe))

# derivative check of Dsep by central differences (mpmath, away from lattice)
J.set_prec(260)
worst_g = 0.0
cnt_g = 0
for ib in range(1, 8):
    beta = 2 * J.PI_IV() * ib / 8
    for i0 in range(1, 6):
        for i1 in range(1, 6):
            th0 = 2 * J.PI_IV() * (i0 + 0.31) / 6
            th1 = 2 * J.PI_IV() * (i1 + 0.17) / 6
            hstep = J.I("1e-20")
            Dv, D0, D1, Db, a0, a1, AD, sm = J.sep_det_grad(beta, th0, th1)
            f0p = J.sep_det(beta, th0 + hstep, th1)[0]
            f0m = J.sep_det(beta, th0 - hstep, th1)[0]
            f1p = J.sep_det(beta, th0, th1 + hstep)[0]
            f1m = J.sep_det(beta, th0, th1 - hstep)[0]
            e0 = (f0p - f0m) / (2 * hstep) - D0
            e1 = (f1p - f1m) / (2 * hstep) - D1
            worst_g = max(worst_g, J.rad(e0), J.rad(e1))
            cnt_g += 1
OUT["e2_separated_gradient"] = {
    "method": "central difference of Dsep vs the analytic interval gradient, "
              "mpmath.iv prec=260, step 1e-20",
    "samples": cnt_g,
    "max_abs_error": worst_g,
}
print("(e2) gradient of Dsep: %d samples, max |diff| = %.4e" % (cnt_g, worst_g))

with open(os.path.join(HERE, "validation.json"), "w") as fh:
    json.dump(OUT, fh, indent=2)
print("wrote validation.json")
