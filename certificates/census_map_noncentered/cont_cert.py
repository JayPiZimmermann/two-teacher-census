"""Certified grid-point tracking of selected census-map zero families.

The engine and the emitter.  `cont_check.py` is the replay checker; the two
share every deterministic rule through this module, exactly as
`foldfree_sweep.py` / `foldfree_check.py` do, so that a replay RE-DERIVES each
verdict instead of reading it.

--------------------------------------------------------------------------
WHAT ONE STEP IS
--------------------------------------------------------------------------
A step carries one tracked zero of the census angle map from a teacher `p` at
which it is already certified to a neighbouring teacher `p'`:

  1. PREDICT.  A DAMPED Newton corrector for the row-form residual
     `(E0, E1)` in the `(s, D)` chart (`th0 = s + D`, `th1 = s`), run at the
     WORKING PRECISION on thin intervals, starting from the midpoint of the
     enclosure certified at `p`.  This is a guess and nothing else: it enters
     no verdict, and a bad guess can only make the certification fail.
  2. CERTIFY.  A 2x2 Krawczyk test (`census_cert._krawczyk_sD`) on the box of
     radius `r` about the predicted point, `r` halved from `r0` until either
     `K(Z) subset int(Z)` -- EXACTLY ONE zero of the map in `Z`, enclosed --
     or `r < rmin`, which is an UNDECIDED verdict.
  3. READ THE TYPE.  On the Krawczyk-refined enclosure: the interval Schur
     type of the census (`census_cert.schur_label`, the `widgets.js` label in
     interval arithmetic), and the sign of the ANGLE-MAP JACOBIAN DETERMINANT.

The Jacobian determinant is the one the Lean schema consumes.  In the `(s,D)`
chart the row form's own derivative gives it with no new arithmetic:

      det_theta(E0,E1) = dE0/dD * dE1/ds  -  dE0/ds * dE1/dD

(the chain rule `d/dth0 = d/dD`, `d/dth1 = d/ds - d/dD` has determinant 1),
and the tree's `separatedAngleJacDetJ` is exactly `-det_theta(E0,E1)` because
the tree's second equation is the negative of this file's `E1`
(`generalJAngleEq1` carries `couplingHJ (th1 - th0)`, and `couplingHJ` is odd).
So a SIGN-DEFINITE enclosure of `det_theta` over the solution box is
`separatedAngleJacDetJ != 0` at the enclosed zero -- the schema's `hfold`, at
that zero.  `cont_gate.py` is the falsifier for both identities.

--------------------------------------------------------------------------
WHY A BISECTION TREE, AND OVER WHAT
--------------------------------------------------------------------------
The steps of a walk along one parameter line are the LEAVES OF A BISECTION
TREE OVER THAT LINE'S PARAMETER INTERVAL, and preorder traversal is exactly
walk order.  A node is a parameter interval `[u0, u1]` in the line's
normalised coordinate; `1` = internal, meaning the step was refused and the
walk goes `u0 -> m -> u1` through the midpoint `m = (u0+u1)/2`; `0` = leaf,
meaning the walk stepped straight to `u1` and the three following bits are the
verdict there.  The midpoint rule is determined by the node, so a reader that
knows the root cells reconstructs every parameter exactly and no coordinate
appears in the file -- the same discipline as `face_bb_signlaw.py`.

Root cells are the GRID intervals, so every grid parameter is the right
endpoint of some leaf and therefore carries a verdict.  Bisection adds
certified teachers between grid points; it never removes one.

Lines, and their order in the bitstream, for each tracked family in turn:

    ANCHOR       a degenerate leaf: the certification at (beta[0], y[0])
                 from the family's declared float seed
    SPINE        a walk in `y` at `beta = beta[0]`, root cells the y-grid
                 intervals, anchored at the ANCHOR
    TOOTH j      a walk in `beta` at `y = y[j]`, root cells the beta-grid
                 intervals, anchored at the SPINE's certified state at y[j]

so the comb covers the whole (nb+1) x (ny+1) grid with one certification per
teacher and per family.

--------------------------------------------------------------------------
WHAT A VALID REPLAY ESTABLISHES -- AND WHAT IT DOES NOT
--------------------------------------------------------------------------
ESTABLISHED, for the CONCRETE census map and for the teachers the spec names
and no others: at every teacher of the regenerated grid, and for each tracked
family, there is EXACTLY ONE zero of the census angle map in the certified box
(Krawczyk), the boxes of distinct families are pairwise disjoint, and at that
zero `separatedAngleJacDetJ != 0` and the interval Schur type is the reported
one.

NOT ESTABLISHED: the `hfold` hypothesis of
`censusAngleMapJ_zeroCount_eq_witness`.  That hypothesis quantifies over ALL
zeros in the compact student region at ALL teachers of the face.  A
continuation certifies fold-freeness ALONG THE SHEETS IT TRACKS.  It cannot by
itself exclude a zero that is not a continuation of a seeded family, nor the
birth of a new pair between two grid teachers.  There is a bootstrap -- if
fold-freeness held at every zero the count would be constant, so every zero
would be a continuation of one at the witness -- and it is CIRCULAR until
degenerate zeros are excluded independently.  Nothing in this file excludes
them, and no claim here should be read as if it did.

Not established either: anything between grid teachers.  A step is a statement
about its endpoint, not about the segment.

``CERTIFICATE.md`` (status 2026-08-24) is the authoritative claim ledger for
these replay artifacts.

Run:
    python3 cont_cert.py <specname>          # emit  cont_<tag>.json
    python3 cont_check.py cont_<tag>.json    # replay, exit 0 = valid+complete
"""
import base64
import json
import os
import sys
import time
from fractions import Fraction

import mpmath
from mpmath import iv

import census_cert as X
import noncentered as J

HERE = os.path.dirname(os.path.abspath(__file__))

# verdict codes -- three bits, as in ivcert
CERT_SADDLE = 0        # unique zero, jacdet sign-definite, Schur type `saddle`
CERT_SPURIOUS = 1      # unique zero, jacdet sign-definite, Schur type `spurious`
CERT_FOLDFREE = 2      # unique zero, jacdet sign-definite, Schur type undecided
FOLD_CANDIDATE = 3     # unique zero but jacdet enclosure STRADDLES zero
UNDECIDED_JUMP = 4     # certified, but the enclosure is further from the
                       # previous one than `maxjump`: the walk did not stay on
                       # its sheet, so the step is REFUSED and the parameter
                       # interval is bisected
UNDECIDED_MINRAD = 6   # Krawczyk never certified down to rmin
UNDECIDED_OTHER = 7
NAME = {0: "CERT_SADDLE", 1: "CERT_SPURIOUS", 2: "CERT_FOLDFREE",
        3: "FOLD_CANDIDATE", 4: "UNDECIDED_JUMP", 5: "RESERVED5",
        6: "UNDECIDED_MINRAD", 7: "UNDECIDED_OTHER"}
ACCEPTING = (CERT_SADDLE, CERT_SPURIOUS, CERT_FOLDFREE)


# --------------------------------------------------------------------------
# grids -- regenerated from the spec, never read from the file
# --------------------------------------------------------------------------

def _mpf(v):
    return mpmath.mpf(repr(v)) if not isinstance(v, mpmath.mpf) else v


def grid(spec, axis):
    """The (n+1) grid values of one axis, as exact mpf from the spec."""
    lo, hi, n = _mpf(spec[axis + "0"]), _mpf(spec[axis + "1"]), spec["n" + axis]
    return [lo + (hi - lo) * mpmath.mpf(k) / mpmath.mpf(n) for k in range(n + 1)]


def param_at(spec, axis, u):
    """The parameter value at normalised coordinate `u` (an exact Fraction)."""
    lo, hi = _mpf(spec[axis + "0"]), _mpf(spec[axis + "1"])
    return lo + (hi - lo) * mpmath.mpf(u.numerator) / mpmath.mpf(u.denominator)


def teacher(beta_m, y_m):
    """The census teacher at (beta, y): masses (sin psi, cos psi)."""
    B = iv.mpf([beta_m, beta_m])
    Y = iv.mpf([y_m, y_m])
    s0, s1 = X.masses_at(Y)
    return X.Teacher(B, s0, s1)


# --------------------------------------------------------------------------
# one step
# --------------------------------------------------------------------------

def thin(m):
    return iv.mpf([m, m])


def _res_mid(Tc, s, D):
    r = X.sep_res_rowform_jac(Tc, thin(s) + thin(D), thin(s), None,
                              Dm=thin(D))
    return [J.mid(z) for z in r]


def newton_correct(Tc, s, D, iters, tol, trust, delta):
    """DAMPED Newton corrector for the row-form residual in the `(s, D)`
    chart, at the working precision.  A PREDICTOR: its output is a guess,
    never a verdict, and a bad guess can only make a certification fail.

    Two dampings, both needed, both measured (2026-08-09):

    * a TRUST RADIUS on the `(s, D)` step, plus backtracking until the
      residual `|E0| + |E1|` does not increase.  The undamped step is
      `J^-1 E`, and `det J` at the exact fit falls like a high power of
      `beta` as the fit merges into the coincidence stratum -- measured
      `-5.5e-5` at `beta = 0.10`, `-3.3e-3` at `0.20`, `-0.61` at `0.50`,
      `-20.8` at `1.00`.  A full step from the exact fit at `beta = 0.15` to
      `beta = 0.20` therefore has length `3.2` and lands on ANOTHER family --
      which is precisely what the first F4 run did, and what the
      disjointness check caught.
    * a GAP FLOOR `delta`, refusing any iterate with `|D| < delta`.  The row
      form `E = massDet(D) * G` vanishes IDENTICALLY on `D = 0`
      (`massDet(0) = 0` and `h(0) = 0`), so the coincidence line is a
      SPURIOUS zero manifold of `E` and an undamped corrector is attracted to
      it.  `delta` is the same validity-domain object as the enumeration's
      gap floor, not a tuning knob.
    """
    for _ in range(iters):
        E0, E1, a, b, c, d = _res_mid(Tc, s, D)
        n0 = abs(E0) + abs(E1)
        det = a * d - b * c
        if det == 0:
            break
        ds = (-E0 * d + E1 * b) / det
        dD = (-E1 * a + E0 * c) / det
        L = abs(ds) + abs(dD)
        if L > trust:
            ds, dD = ds * trust / L, dD * trust / L
        t = mpmath.mpf(1)
        accepted = False
        for _ in range(40):
            Dn = D + t * dD
            if abs(Dn) >= delta:
                q = _res_mid(Tc, s + t * ds, Dn)
                if abs(q[0]) + abs(q[1]) <= n0:
                    accepted = True
                    break
            t = t / 2
        if not accepted:
            break
        s, D = s + t * ds, D + t * dD
        if t * (abs(ds) + abs(dD)) < tol:
            break
    return s, D


def _branches(Tc, S, DD):
    """Single-half-branch data for the box, or None if a kink is crossed."""
    bb = {}
    T0 = S + DD
    for key, Xb in (("n0", T0), ("n1", S), ("n0b", T0 - Tc.B),
                    ("n1b", S - Tc.B), ("nD", DD)):
        ps = J.branch_pieces(Xb)
        if len(ps) != 1:
            return None
        bb[key] = ps[0][1]
    return bb


def jac_det_theta(Tc, T0, T1, br, Dm):
    """`det_theta(E0,E1) = dE0/dD * dE1/ds - dE0/ds * dE1/dD`, the angle-map
    Jacobian determinant up to the sign of the tree's second equation.  Gated
    by `cont_gate.py` against central differences AND against a transcription
    of the tree's own `separatedAngleJacDetJ` closed forms."""
    r = X.sep_res_rowform_jac(Tc, T0, T1, br, Dm=Dm)
    return r[3] * r[4] - r[2] * r[5]


def certify(spec, Tc, s, D):
    """Krawczyk certification at a predicted point.  Returns
    (code, S, DD, schur, jsign, radius)."""
    r = _mpf(spec["r0"])
    rmin = _mpf(spec["rmin"])
    while r >= rmin:
        S = iv.mpf([s - r, s + r])
        DD = iv.mpf([D - r, D + r])
        br = _branches(Tc, S, DD)
        k, N0, N1 = X._krawczyk_sD(Tc, S, DD, br)
        if k == "unique":
            for _ in range(spec.get("refine", 80)):
                kk, M0, M1 = X._krawczyk_sD(Tc, N0, N1, br)
                if kk != "unique":
                    break
                if J.width(M0) >= J.width(N0) and J.width(M1) >= J.width(N1):
                    break
                N0, N1 = M0, M1
            R1, DDf = N0, N1
            R0 = R1 + DDf
            rr = X.sep_res(Tc, R0, R1, br, Dm=DDf)
            if rr is None:
                return UNDECIDED_OTHER, None, None, None, 0, float(r)
            lab, _sg = X.schur_label(Tc, R0, R1, rr[2], rr[3])
            jsign = X.sgn(jac_det_theta(Tc, R0, R1, br, DDf))
            if jsign == 0:
                return FOLD_CANDIDATE, R1, DDf, lab, 0, float(r)
            if lab == "saddle":
                code = CERT_SADDLE
            elif lab == "spurious":
                code = CERT_SPURIOUS
            else:
                code = CERT_FOLDFREE
            return code, R1, DDf, lab, jsign, float(r)
        r = r / 2
    return UNDECIDED_MINRAD, None, None, None, 0, None


def predict(spec, Tc, s_prev, D_prev):
    return newton_correct(Tc, s_prev, D_prev, spec.get("newton_iters", 40),
                          mpmath.mpf(10) ** (-spec.get("newton_tol_exp", 40)),
                          _mpf(spec["trust"]), _mpf(spec["delta"]))


def step(spec, beta_m, y_m, s_prev, D_prev):
    """PREDICT, CERTIFY, then apply the SHEET GUARD at the teacher
    `(beta_m, y_m)`, starting from the midpoint of the previous certified
    enclosure.

    The guard is part of the deterministic rule, not a diagnostic: a
    certified box further than `maxjump` from the previous one is a step that
    left its sheet, and the verdict is `UNDECIDED_JUMP`, which makes the
    caller BISECT the parameter interval and try again on half the step.  The
    box it refuses is a perfectly good enclosure of a perfectly good zero --
    what fails is the reading of the sequence as ONE tracked family, and that
    reading is the whole point."""
    Tc = teacher(beta_m, y_m)
    s, D = predict(spec, Tc, s_prev, D_prev)
    code, S, DD, lab, jsign, rad = certify(spec, Tc, s, D)
    if code in ACCEPTING:
        disp = abs(J.mid(S) - s_prev) + abs(J.mid(DD) - D_prev)
        if disp > _mpf(spec["maxjump"]):
            return UNDECIDED_JUMP, S, DD, lab, jsign, rad
    return code, S, DD, lab, jsign, rad


# --------------------------------------------------------------------------
# the walk over one parameter line, as a bisection tree
# --------------------------------------------------------------------------

class Bits(object):
    """A bit sink (emitter) or source (checker)."""

    def __init__(self, bits=None):
        self.out = [] if bits is None else None
        self.inp = bits
        self.pos = 0

    def put(self, b):
        self.out.append(b)

    def putn(self, v, n):
        for k in range(n - 1, -1, -1):
            self.out.append((v >> k) & 1)

    def get(self):
        if self.inp is None or self.pos >= len(self.inp):
            raise IndexError("bitstream exhausted")
        b = self.inp[self.pos]
        self.pos += 1
        return b

    def getn(self, n):
        v = 0
        for _ in range(n):
            v = (v << 1) | self.get()
        return v


class Walk(object):
    """One parameter line.  `axis` is 'b' or 'y'; `fixed` is the other
    coordinate's value; `state` is the running (s, D) midpoint pair."""

    def __init__(self, spec, axis, fixed, state, bits, emit, stats, record):
        self.spec = spec
        self.axis = axis
        self.fixed = fixed
        self.state = state
        self.bits = bits
        self.emit = emit
        self.stats = stats
        self.record = record          # callback(u, S, DD) at grid coordinates
        self.min_u = Fraction(1, spec.get("min_split", 64) * spec["n" + axis])

    def teacher_at(self, u):
        p = param_at(self.spec, self.axis, u)
        return (p, self.fixed) if self.axis == "b" else (self.fixed, p)

    def leaf_here(self, u):
        """Certify at `u` from the current state; emit/consume the verdict."""
        beta_m, y_m = self.teacher_at(u)
        if self.emit:
            code, S, DD, lab, jsign, rad = step(self.spec, beta_m, y_m,
                                                self.state[0], self.state[1])
            self.bits.putn(code, 3)
        else:
            code = self.bits.getn(3)
            if code not in ACCEPTING:
                self.stats["undec"] += 1
                self.stats["undec_where"].append((self.axis, str(u), NAME[code]))
                return False
            got, S, DD, lab, jsign, rad = step(self.spec, beta_m, y_m,
                                               self.state[0], self.state[1])
            if got != code:
                self.stats["bad"] += 1
                if len(self.stats["bad_where"]) < 20:
                    self.stats["bad_where"].append(
                        (self.axis, str(u), NAME[code], NAME[got]))
                return False
        if code not in ACCEPTING:
            self.stats["undec"] += 1
            self.stats["undec_where"].append((self.axis, str(u), NAME[code]))
            return False
        self.accept(u, code, S, DD)
        return True

    def accept(self, u, code, S, DD):
        """Book a certified leaf and move the walk's state onto it.  The
        DISPLACEMENT is recorded as a diagnostic only: soundness of what the
        certificate claims (one zero per box, boxes disjoint, Jacobian
        nonzero) does not depend on the walk staying on one sheet, but the
        reading of the boxes AS A TRACKED FAMILY does, and a jump would show
        up here as a large step."""
        self.stats["ok"] += 1
        self.stats["types"][NAME[code]] = \
            self.stats["types"].get(NAME[code], 0) + 1
        new = (J.mid(S), J.mid(DD))
        d = float(abs(new[0] - self.state[0]) + abs(new[1] - self.state[1]))
        if d > self.stats["max_disp"]:
            self.stats["max_disp"] = d
        self.state = new
        self.record(u, S, DD)

    def node(self, u0, u1):
        """A node of the tree over [u0, u1]: step from u0 to u1, bisecting the
        parameter interval when the step is refused."""
        if self.emit:
            saved = self.state
            beta_m, y_m = self.teacher_at(u1)
            code, S, DD, lab, jsign, rad = step(self.spec, beta_m, y_m,
                                                saved[0], saved[1])
            if code in ACCEPTING or (u1 - u0) <= self.min_u:
                self.bits.put(0)
                self.bits.putn(code, 3)
                if code not in ACCEPTING:
                    self.stats["undec"] += 1
                    self.stats["undec_where"].append(
                        (self.axis, str(u1), NAME[code]))
                    return False
                self.accept(u1, code, S, DD)
                return True
            self.bits.put(1)
            self.stats["splits"] += 1
            m = (u0 + u1) / 2
            return self.node(u0, m) and self.node(m, u1)
        b = self.bits.get()
        if b == 1:
            if (u1 - u0) <= self.min_u:
                self.stats["bad"] += 1
                self.stats["bad_where"].append(
                    (self.axis, str(u1), "SPLIT BELOW min_split", ""))
                return False
            self.stats["splits"] += 1
            m = (u0 + u1) / 2
            return self.node(u0, m) and self.node(m, u1)
        return self.leaf_here(u1)


# --------------------------------------------------------------------------
# the comb over one family
# --------------------------------------------------------------------------

def run_family(spec, seed, bits, emit, stats):
    """Anchor, spine (in y at beta[0]), then one tooth in beta per y-grid
    value.  Returns {(i, j): (S, DD)} over the grid, or None on failure."""
    nb, ny = spec["nb"], spec["ny"]
    bg, yg = grid(spec, "b"), grid(spec, "y")
    boxes = {}

    # --- ANCHOR: a degenerate leaf at (beta[0], y[0]) from the declared seed
    state = (_mpf(seed[0]), _mpf(seed[1]))
    if emit:
        code, S, DD, lab, jsign, rad = step(spec, bg[0], yg[0],
                                            state[0], state[1])
        bits.put(0)
        bits.putn(code, 3)
    else:
        if bits.get() != 0:
            stats["bad"] += 1
            stats["bad_where"].append(("anchor", "0", "INTERNAL NODE", ""))
            return None
        code = bits.getn(3)
        got, S, DD, lab, jsign, rad = step(spec, bg[0], yg[0],
                                           state[0], state[1])
        if code not in ACCEPTING:
            stats["undec"] += 1
            stats["undec_where"].append(("anchor", "0", NAME[code]))
            return None
        if got != code:
            stats["bad"] += 1
            stats["bad_where"].append(("anchor", "0", NAME[code], NAME[got]))
            return None
    if code not in ACCEPTING:
        stats["undec"] += 1
        stats["undec_where"].append(("anchor", "0", NAME[code]))
        return None
    stats["ok"] += 1
    stats["types"][NAME[code]] = stats["types"].get(NAME[code], 0) + 1
    boxes[(0, 0)] = (S, DD)
    spine_state = {0: (J.mid(S), J.mid(DD))}

    # --- SPINE: walk in y at beta[0]
    def rec_spine(u, S, DD):
        num = u * ny
        if num.denominator == 1:
            j = int(num)
            boxes[(0, j)] = (S, DD)
            spine_state[j] = (J.mid(S), J.mid(DD))

    w = Walk(spec, "y", bg[0], (J.mid(S), J.mid(DD)), bits, emit, stats,
             rec_spine)
    for j in range(ny):
        if not w.node(Fraction(j, ny), Fraction(j + 1, ny)):
            return None

    # --- TEETH: walk in beta at each y[j]
    for j in range(ny + 1):
        if j not in spine_state:
            stats["bad"] += 1
            stats["bad_where"].append(("spine", str(j), "NO STATE", ""))
            return None

        def rec_tooth(u, S, DD, j=j):
            num = u * nb
            if num.denominator == 1:
                boxes[(int(num), j)] = (S, DD)

        w = Walk(spec, "b", yg[j], spine_state[j], bits, emit, stats, rec_tooth)
        for i in range(nb):
            if not w.node(Fraction(i, nb), Fraction(i + 1, nb)):
                return None
    return boxes


# --------------------------------------------------------------------------
# disjointness of the tracked families
# --------------------------------------------------------------------------

def disjointness(spec, per_family):
    """At every grid teacher, the tracked families' enclosures must be
    pairwise disjoint AS UNORDERED STUDENT PAIRS.  Two boxes are disjoint if
    the `s` boxes are disjoint modulo 2pi, or the `D` boxes are disjoint; the
    swap `(s, D) -> (s + D, 2pi - D)` is applied to the second before the
    comparison as well, so a family and the mirror copy of another are not
    read as distinct."""
    nb, ny = spec["nb"], spec["ny"]
    tp = 2 * J.PI_IV()
    worst = None
    for i in range(nb + 1):
        for j in range(ny + 1):
            got = [(f, d[(i, j)]) for f, d in enumerate(per_family)
                   if (i, j) in d]
            for a in range(len(got)):
                for b in range(a + 1, len(got)):
                    (Sa, Da), (Sb, Db) = got[a][1], got[b][1]
                    seps = []
                    for (S2, D2) in ((Sb, Db), (Sb + Db, tp - Db)):
                        ds = _sep_mod(Sa, S2, tp)
                        dd = _sep_mod(Da, D2, tp)
                        seps.append(max(ds, dd))
                    sep = min(seps)
                    if worst is None or sep < worst[0]:
                        worst = (sep, i, j, got[a][0], got[b][0])
                    if sep <= 0:
                        return False, (i, j, got[a][0], got[b][0]), worst
    return True, None, worst


def _sep(A, B):
    a0, a1 = J.endpoints(A)
    b0, b1 = J.endpoints(B)
    return float(max(a0 - b1, b0 - a1))


def _sep_mod(A, B, tp):
    """Separation of two boxes on the circle of circumference 2pi: the box `B`
    has a copy at every `2pi` shift, so the boxes are disjoint only if `A`
    misses EVERY copy.  Hence the MINIMUM over the shifts -- taking the maximum
    would read `A` as disjoint from `B` whenever it misses the shifted copy,
    which is always."""
    worst = None
    for k in (-1, 0, 1):
        sh = iv.mpf(k) * tp
        v = _sep(A, B + sh)
        worst = v if worst is None else min(worst, v)
    return worst


# --------------------------------------------------------------------------
# emit
# --------------------------------------------------------------------------

def new_stats():
    return {"ok": 0, "bad": 0, "undec": 0, "splits": 0, "types": {},
            "max_disp": 0.0, "bad_where": [], "undec_where": []}


def emit_certificate(spec, tag):
    J.set_prec(spec.get("prec", 160))
    bits = Bits()
    stats = new_stats()
    per_family = []
    t0 = time.time()
    for seed in spec["seeds"]:
        d = run_family(spec, seed, bits, True, stats)
        if d is None:
            break
        per_family.append(d)
    secs = time.time() - t0
    ok, clash, worst = disjointness(spec, per_family)
    raw = bytearray()
    for k in range(0, len(bits.out), 8):
        byte = 0
        for b in bits.out[k:k + 8]:
            byte = (byte << 1) | b
        byte <<= (8 - len(bits.out[k:k + 8])) % 8
        raw.append(byte)
    doc = {
        "object": "pointwise Krawczyk certification of tracked census-map "
                  "zero families at a finite teacher grid",
        "authoritative_ledger": {
            "path": "CERTIFICATE.md",
            "status_date": "2026-08-24",
        },
        "tag": tag,
        "spec": spec,
        "n_nodes": len(bits.out),
        "bits_b64": base64.b64encode(bytes(raw)).decode(),
        "emit": {
            "seconds": secs,
            "families_completed": len(per_family),
            "certified_steps": stats["ok"],
            "parameter_splits": stats["splits"],
            "undecided": stats["undec"],
            "verdicts": stats["types"],
            "ms_per_step": 1000.0 * secs / max(stats["ok"], 1),
            "max_step_displacement": stats["max_disp"],
            "disjoint": ok,
            "worst_separation": None if worst is None else worst[0],
            # where the worst separation is attained, and where the first
            # overlap is if there is one: an emitter that reports only the
            # BOOLEAN sends its reader into a full replay to find out which
            # two families met and at which teacher.
            "worst_where": None if worst is None
                           else {"grid_i": worst[1], "grid_j": worst[2],
                                 "family_a": worst[3], "family_b": worst[4]},
            "clash": None if clash is None
                     else {"grid_i": clash[0], "grid_j": clash[1],
                           "family_a": clash[2], "family_b": clash[3]},
        },
        "NOT_CLAIMED": "This certifies fold-freeness ALONG THE TRACKED SHEETS "
                       "at the grid teachers only.  It is NOT the `hfold` "
                       "hypothesis of censusAngleMapJ_zeroCount_eq_witness, "
                       "which quantifies over ALL zeros at ALL teachers of the "
                       "face; continuation cannot exclude a zero that is not a "
                       "continuation of a seeded family, nor the birth of a "
                       "new pair, and the count bootstrap is circular until "
                       "degenerate zeros are excluded independently.",
    }
    return doc, stats


SPECS = {
    # F4: mixed sector, middle band, outside the torque lens. The rectangle is
    # a regional certificate box inside the proposed F4 label; the sampled
    # candidate transitions are not globally enclosed.
    # Two tracked families: the single separated family of F4, and the exact
    # fit th = (beta, 0), which is a zero of the angle map at every teacher.
    "F4": {
        "face": "F4",
        # The rectangle.  F4 is `0 < beta < beta*_1 = 2.2257` (outside the
        # torque lens) and `-1-w_a(beta) < y < w_a(beta)`; a located family
        # count of 1 (as against 3 in F2 above and F6 below) is observed over
        # `y in [-0.09,-0.91]` at every beta probed. This rectangle stays away
        # from the sampled count-change brackets and the degenerate beta ends;
        # it is not a certificate for an entire global face.
        # The beta floor is not cosmetic: at the degenerate column `beta = 0`
        # the exact fit merges into the coincidence stratum and its own
        # Jacobian determinant falls like `beta^5.9`, so the PREDICTOR needs a
        # parameter step of that order.  The sign test does not -- the
        # determinant is enclosed to `~1e-40` against a value of `5.5e-5` at
        # `beta = 0.10` -- which is why the wall here is conditioning of the
        # walk and not reach of the enclosure.
        "b0": 0.10, "b1": 2.20, "nb": 42,        # beta step 0.05
        "y0": -0.10, "y1": -0.90, "ny": 16,      # y step 0.05
        # located at the anchor teacher (0.10, -0.10) and declared here: the
        # seeds are spec INPUT, and a wrong seed can only make a leaf fail.
        "seeds": [[3.176223, 3.120344],          # the separated family of F4
                  [0.0, 0.10]],                  # the exact fit, th = (beta, 0)
        "r0": 1e-3, "rmin": 1e-11, "min_split": 64,
        "trust": 0.08, "delta": 0.02, "maxjump": 0.35,
        "refine": 80, "newton_iters": 40, "newton_tol_exp": 40,
        "prec": 160,
    },
}

SPECS["F6"] = {
    "face": "F6",
    # F6: mixed sector, the band above `y = -1`, outside the torque lens.  It
    # is THIN -- a located family count of 3 (against 1 in F4 above) is
    # observed for `y` between about `-0.915` and `-1`, so the y-grid step
    # here is `0.02`, not the `0.05` F4 uses.  Four tracked families: the
    # three separated families of F6, and the exact fit.
    "b0": 0.30, "b1": 2.20, "nb": 38,        # beta step 0.05
    "y0": -0.94, "y1": -0.98, "ny": 2,       # y step 0.02
    "seeds": [[0.273683, 0.625841],
              [0.275053, 1.359207],
              [0.276300, 3.063624],
              [0.0, 0.30]],
    "r0": 1e-3, "rmin": 1e-11, "min_split": 64,
    "trust": 0.08, "delta": 0.02, "maxjump": 0.35,
    "refine": 80, "newton_iters": 40, "newton_tol_exp": 40,
    "prec": 160,
}

SPECS["F1a"] = {
    "face": "F1 (open end)",
    # F1 is the SAME-SIGN sector `0 < beta < pi, 0 < y < 1`, and it is the one
    # face the covering sweep never closed.  Continuation does not use the
    # sweep, so the obstruction there does not transfer; what it does have to
    # respect is where F1's own thresholded quantities VANISH, and both of the
    # rectangle's four walls are set by a vanishing, not by a preference.
    #
    #   * `y -> 0` and `y -> 1` (the two teacher-mass strata).  F1 has FOUR
    #     separated families and the smallest of their student gaps goes to
    #     zero LINEARLY at each stratum: measured `D_min / y` = 0.115 at
    #     beta = 0.30, 0.26 at 0.50, 0.53 at 0.75, 0.87 at 1.00 for y down to
    #     0.008 (`work/contF1/scan1.log`, and `wall_pos.json` from the earlier
    #     wall campaign, whose "locator_failed" rows are exactly this family
    #     failing to certify as its gap closes).  So NO positive gap floor
    #     covers a rectangle that reaches either stratum -- `delta` is a named
    #     boundary of the claim, not a knob -- and the collar is 0.10 wide
    #     here because 0.10 is where the worst-case gap (at beta = 0.30)
    #     is still 0.0124, i.e. 2.5x the floor used.
    #   * `beta -> pi` (the degenerate column).  The same measurement in the
    #     other coordinate: `D_min / (pi - beta)` = 0.66 at y = 0.5, flat over
    #     beta = 3.05 and 3.10, and the four families crowd together as well
    #     (pairwise separation 0.0182 at beta = 3.00 against 0.132 at 0.30).
    #     The rectangle stops at 3.00, leaving `pi - beta < 0.1416`.
    #   * `beta -> 0`.  As in F4, the exact fit merges into the coincidence
    #     stratum and the walk's conditioning, not the enclosure, is the wall.
    #
    # Four separated families plus the exact fit = FIVE tracked families, the
    # most of any face here.
    "b0": 0.30, "b1": 2.60, "nb": 46,        # beta step 0.05
    "y0": 0.10, "y1": 0.90, "ny": 16,        # y step 0.05
    # located at the anchor teacher (0.30, 0.10) and declared here.
    "seeds": [[3.033032, 0.191329],
              [3.080304, 0.012398],
              [3.220976, 3.013292],
              [4.579480, 1.652807],
              [0.0, 0.30]],                  # the exact fit, th = (beta, 0)
    "r0": 1e-3, "rmin": 1e-11, "min_split": 64,
    # `delta` is below the rectangle's own worst measured gap (0.0124) by
    # 2.5x.  `maxjump` is set by the rule the F3/F7 pair priced: BELOW THE
    # FAMILY SEPARATION on the rectangle it runs on.  On F1 that rule is what
    # forces the SPLIT at `beta = 2.60`.  The four families' worst pairwise
    # separation over the rectangle runs 0.132 at `beta = 0.30`, 0.088 at
    # 2.60, 0.050 at 2.80 and 0.018 at 3.00 -- a factor of seven -- and a
    # single guard is then either useless at the crowded end or ruinous at the
    # open one.  Measured, twice: at `maxjump = 0.15` the walk emitted 4942
    # leaves that all replay bit-exactly and still failed the end-to-end
    # disjointness check (families 2 and 3 on ONE zero at `(2.95, 0.30)`,
    # separation `-2.6e-44`); at 0.05 it emitted 9497 and failed the same
    # check at `(2.95, 0.65)`.  So `F1a` carries the open end of the face at a
    # guard of 0.04 and `F1b` the crowded end at 0.008, and their union is the
    # rectangle above.  The per-step verdict cannot see a sheet jump -- only
    # the end-to-end check can -- which is F4's lesson at a face where it
    # bites twice.
    "trust": 0.04, "delta": 0.005, "maxjump": 0.04,
    "refine": 80, "newton_iters": 40, "newton_tol_exp": 40,
    "prec": 160,
}

SPECS["F1b"] = {
    "face": "F1 (crowded end)",
    # The other half of the split described in `F1a`: `beta` from 2.60 to
    # 3.00, where F1's four families close from 0.088 to 0.018 apart and the
    # sheet guard has to come down with them.  Same `y` range, same seeds
    # RELOCATED at this rectangle's own anchor `(2.60, 0.10)`, same gap floor.
    # `beta` stops at 3.00: at 3.10 the four gaps are 0.0056 / 0.0095 / 0.0397
    # / 0.0573 at `y = 0.05` and the closest pair is 0.0028 apart, which is
    # the degenerate column asserting itself and not a tractability question.
    "b0": 2.60, "b1": 3.00, "nb": 8,         # beta step 0.05
    "y0": 0.10, "y1": 0.90, "ny": 16,        # y step 0.05
    # located at the anchor teacher (2.60, 0.10) and declared here.
    "seeds": [[2.995352, 0.222828],
              [3.022747, 0.107640],
              [5.471261, 0.779964],
              [5.789020, 0.452341],
              [0.0, 2.60]],                  # the exact fit, th = (beta, 0)
    "r0": 1e-3, "rmin": 1e-11, "min_split": 64,
    "trust": 0.04, "delta": 0.005, "maxjump": 0.008,
    "refine": 80, "newton_iters": 40, "newton_tol_exp": 40,
    "prec": 160,
}

SPECS["F2"] = {
    "face": "F2",
    # F2: mixed sector, the band below `y = 0`, outside the torque lens -- the
    # `y -> -1-y` teacher-swap mirror of F6, and certified here in its own
    # right rather than inherited from it.  Three separated families and the
    # exact fit.
    #
    # The proposed band is `(w_a(beta),0)`. Its sampled count-change bracket
    # and the exact `y=0` stratum motivate the chosen regional box:
    #   * at `w_a` the located count drops 3 -> 1, between `y = -0.06` and
    #     `-0.08` at every beta from 0.20 to 2.00 and between `-0.08` and
    #     `-0.09` at 2.20 (`work/contF1/probe_F2.log`), so the rectangle stops
    #     at `-0.06`;
    #   * at the `y = 0` stratum the smallest family gap vanishes LINEARLY,
    #     `D_min/|y|` = 9.9 at `beta = 0.20` falling to 5.6 at `2.20`, so at
    #     `y = -0.01` that gap is 0.056 -- 2.8x the floor used, and no
    #     positive floor survives the stratum itself.
    # The beta range stops at 2.20 because `beta*_1 = 2.2257` is where the
    # torque lens begins and this regional box ends.
    "b0": 0.30, "b1": 2.20, "nb": 38,        # beta step 0.05
    "y0": -0.01, "y1": -0.06, "ny": 5,       # y step 0.01
    # located at the anchor teacher (0.30, -0.01) and declared here.
    "seeds": [[3.237426, 3.049914],
              [4.639756, 1.647731],
              [6.187720, 0.100262],
              [0.0, 0.30]],                  # the exact fit, th = (beta, 0)
    "r0": 1e-3, "rmin": 1e-11, "min_split": 64,
    "trust": 0.08, "delta": 0.02, "maxjump": 0.35,
    "refine": 80, "newton_iters": 40, "newton_tol_exp": 40,
    "prec": 160,
}

SPECS["F3"] = {
    "face": "F3",
    # F3: mixed sector, the band below `y = 0`, INSIDE the torque lens.  Three
    # separated families and the exact fit.  This is the face the README's
    # determinant-route table names as unreachable (`D = 0.0315, 0.0597,
    # 0.0961` at its `(3.08, -0.08)` witness, against `D* ~ 0.40`), so a
    # continuation certificate is the only route that reaches it at all.
    #
    # BOTH of the rectangle's limits are a vanishing measured in its own
    # coordinate (`work/contF1/probe_F3.log`):
    #   * `beta -> pi`: the smallest family gap is `D_min ~ c(beta)*|y|` with
    #     `c` = 1.52, 0.89, 0.58, 0.39, 0.26 at `beta` = 2.90, 3.00, 3.05,
    #     3.08, 3.10 -- i.e. `c -> 0` at the degenerate column, and the three
    #     families crowd together with it (pairwise separation 0.0161 at
    #     `(3.10, -0.10)`).  The rectangle stops at 3.10, where the worst gap
    #     is 0.0052 and the floor used is 0.002.
    #   * `y -> 0`: the same product vanishes at the stratum, which is why the
    #     rectangle starts at `y = -0.02` and not at 0.
    # All four of F3's replayed witnesses -- (3.10,-0.06), (3.10,-0.10),
    # (3.08,-0.08), (3.05,-0.09) -- are grid teachers of this rectangle.
    #
    # THE SHEET GUARD IS SET BELOW THE FAMILY SEPARATION, AND THE RULE THAT
    # SETS IT IS MEASURED, not chosen.  Three runs on the F3/F7 mirror pair
    # priced it.  At `maxjump = 0.05` F3 passed the end-to-end disjointness
    # check and F7 did NOT -- two tracked walks landed on ONE zero at
    # `(2.97, -0.90)`, separation `-7.6e-44` -- and the two faces are exact
    # mirrors (`y -> -1-y` composed with `theta -> beta - theta` preserves the
    # gap `D` exactly), so F3's pass at that setting was luck and not
    # evidence.  At `0.02` F7 failed again, at `(3.10, -0.90)`.  The reason is
    # a ratio and not a threshold: over one grid step of `beta = 0.01` near
    # the degenerate column the families MOVE about as far as they are APART
    # (`D` runs 0.0602, 0.0402, 0.0271 at `beta` = 3.05, 3.08, 3.10, while the
    # closest pair there is 0.0126 apart), so no guard at or above the step
    # displacement can separate them.  The guard here is 0.006, below the
    # 0.0126 separation; the walk then REFUSES the 0.01 step and bisects the
    # parameter interval until the displacement fits, which is exactly what
    # the bisection is for and costs only wall clock.
    "b0": 2.90, "b1": 3.10, "nb": 20,        # beta step 0.01
    "y0": -0.02, "y1": -0.10, "ny": 8,       # y step 0.01
    # located at the anchor teacher (2.90, -0.02) and declared here.
    "seeds": [[5.923946, 0.362417],
              [6.056445, 0.230870],
              [6.265864, 0.030401],
              [0.0, 2.90]],                  # the exact fit, th = (beta, 0)
    "r0": 1e-3, "rmin": 1e-11, "min_split": 64,
    "trust": 0.04, "delta": 0.002, "maxjump": 0.006,
    "refine": 80, "newton_iters": 40, "newton_tol_exp": 40,
    "prec": 160,
}

SPECS["F5a"] = {
    "face": "F5a regional box (outside the sampled candidate-transition band)",
    # F5: mixed sector, the middle band, INSIDE the torque lens.  This spec is
    # a DIAGNOSTIC as much as a certificate, and the reason is in the locator's
    # own output.  Over `beta` in [2.30, 3.10] and `y` in [-0.80, -0.20] the
    # located family count is 1 EVERYWHERE except exactly on `y = -1/2`, where
    # it is 3 (`work/contF1/probe_F5.log`, 56 teachers).  `y = -1/2` is the
    # EQUAL-MASS line -- `(s0,s1) = (sin psi, cos psi)` with `psi = (y+1)pi/2`,
    # so `psi = pi/4` -- and the two extra families there are a mirror pair
    # (identical gap `D`) beside an exactly ANTIPODAL one (`D = pi`).  Three of
    # F5's four recorded witnesses sit on that line and the fourth, `(3.00,
    # -0.40)`, has a count of 1. The named witnesses share a census string only
    # because the extra returned families are all `saddle` and the string lists
    # traps; this is not constancy over a global face.
    #
    # The continuation diagnostic bypasses the locator: seed the three
    # families ON the line and walk OFF it. Certification at later sampled
    # teachers establishes those individual zeros. A refused step does not by
    # itself prove degeneracy, a fold, or the location of a global boundary.
    # The measurements below are evidence for a candidate transition and
    # motivate treating F5a/F5b as separate regional certificate boxes; they
    # do not prove that the proposed global F5 label is one face or two.
    # Walking the three equal-mass-line families upward in `y` by certified
    # Krawczyk steps of 0.002, TWO of them stop at the SAME `y` at every beta
    # tried while the third walks on to `y = -0.30` untroubled
    # (`work/contF1/f5_foldscan.log`, regenerated 2026-08-10):
    #
    #     beta   last certified y   D of the two that meet   det_theta there
    #     2.30   -0.4980            3.0722 / 2.8027          -19.09 / +19.81
    #     2.60   -0.4580            2.5376 / 2.4975          -3.072 / +2.973
    #     2.90   -0.4240            2.3819 / 2.1120          -15.45 / +10.71
    #     3.00   -0.4160            2.2372 / 2.1506          -3.897 / +3.424
    #     3.10   -0.4140            2.3113 / 2.0100          -14.47 / +8.956
    #
    # The two tracked zeros approach in `D` with Jacobian determinants of
    # opposite sign and decreasing magnitude, behavior consistent with a
    # saddle-node. At `beta = 2.60` the
    # antipodal branch's determinant runs -117.5, -105.1, -86.2, -61.8, -23.2,
    # -3.07 at `y = -0.498, -0.490, -0.480, -0.470, -0.460, -0.458`, i.e.
    # linear extrapolation of `det^2` in `y` has a zero at `y = -0.45796`
    # (stable to 1e-5
    # whether fitted on the last two steps or the last four), just ABOVE the
    # last step that certifies.  This is the campaign's `|f| / |f'|` diagnostic
    # and is consistent with an approaching degeneracy. No accepted enclosure
    # in this diagnostic reaches past it.
    #
    # THE STOPPING HEIGHT IS INSTRUMENT-DEPENDENT AND THE STEPS ARE NOT, so the
    # instrument is stated with the numbers.  `f5_foldscan.py` takes the
    # certify settings from `F5b` below and DISABLES the sheet guard: `maxjump`
    # is a tracking device whose design rule (guard below the tracked families'
    # pairwise separation) is unsatisfiable at a fold, where that separation
    # goes to zero at a fold candidate. With the guard off, every walk above
    # stops on UNDECIDED_MINRAD: Krawczyk certification failed at the minimum
    # radius, leaving a possible degeneracy undecided.
    # With `F5b`'s own `maxjump = 0.10` the walk instead stops one or two steps
    # early at four of these five betas on UNDECIDED_JUMP (a zero WAS certified
    # there; only the sheet reading was refused), splitting the two branches by
    # a step at `beta = 3.00` (-0.4180 / -0.4200) and `3.10` (-0.4160 / -0.4180)
    # where guard-off they stop together, and is refused at its FIRST step at
    # `beta = 2.30`.  Superseded by this run and kept so the discrepancy is
    # traceable: the four-row version of the table above, which came from a
    # pre-split spec whose guard was non-binding, agrees exactly; its
    # `beta = 3.00` companion row `-0.418` and the extrapolation `y = -0.4585`
    # were the guarded instrument's, and `-0.4585` sits BELOW a height at which
    # both families are certified, which is what exposed the mixture.
    #
    # Thus the sampled evidence supports a candidate transition band inside
    # the proposed F5 label and its `y -> -1-y` mirror. The recorded witnesses
    # have family counts 3 near the equal-mass line and 1 at `(3.00,-0.40)`.
    # No global curve, inside/outside decomposition, or failure of `hfold` on a
    # complete F5 face follows from these finite walks.
    #
    # `F5a` and `F5b` are two separately replayed regional boxes chosen on
    # opposite sides of the sampled candidate band. `F5a` stops at
    # `beta=3.00`, `y=-0.40`; extending it requires new regional certificates.
    "b0": 2.40, "b1": 3.00, "nb": 12,        # beta step 0.05
    "y0": -0.40, "y1": -0.15, "ny": 5,       # y step 0.05
    # located at the anchor teacher (2.40, -0.40) and declared here.
    "seeds": [[4.776583, 1.843106],
              [0.0, 2.40]],                  # the exact fit, th = (beta, 0)
    "r0": 1e-3, "rmin": 1e-11, "min_split": 64,
    "trust": 0.04, "delta": 0.02, "maxjump": 0.10,
    "refine": 80, "newton_iters": 40, "newton_tol_exp": 40,
    "prec": 160,
}

SPECS["F5b"] = {
    "face": "F5b regional box (near the equal-mass line)",
    # This independently replayed rectangle surrounds the equal-mass line
    # `y=-1/2`, where the sampled family count is 3. It contains the
    # `(2.60,-0.50)` and `(2.80,-0.50)` witnesses. Attempts to widen it toward
    # the sampled transition return `FOLD_CANDIDATE`; that code is a refusal to
    # certify, not a proof that the refused box contains a fold.
    "b0": 2.60, "b1": 3.10, "nb": 10,        # beta step 0.05
    "y0": -0.52, "y1": -0.48, "ny": 8,       # y step 0.005
    # located at the anchor teacher (2.60, -0.52) and declared here.
    "seeds": [[2.154207, 1.752884],
              [4.481948, 2.970075],
              [4.832063, 2.056670],
              [0.0, 2.60]],                  # the exact fit, th = (beta, 0)
    "r0": 1e-3, "rmin": 1e-11, "min_split": 64,
    "trust": 0.04, "delta": 0.02, "maxjump": 0.10,
    "refine": 80, "newton_iters": 40, "newton_tol_exp": 40,
    "prec": 160,
}

SPECS["F7"] = {
    "face": "F7",
    # F7: mixed sector, the band above `y = -1`, INSIDE the torque lens -- the
    # `y -> -1-y` mirror of F3, certified in its own right.  Everything F3's
    # comment records about the two vanishings holds here with `|y|` replaced
    # by `|1 + y|`, and the same `maxjump` honesty note applies.  All four of
    # F7's replayed witnesses -- (3.10,-0.94), (3.10,-0.90), (3.08,-0.92),
    # (3.05,-0.91) -- are grid teachers of this rectangle.
    #
    # The rectangle is walked from `y = -0.98` toward `-0.90`, i.e. the EXACT
    # mirror of F3's direction, and that is not cosmetic.  Anchored at the
    # other end (`y = -0.90`) the three families' `s` coordinates lie within
    # 0.025 of each other while their gaps are 0.165 / 0.212 / 0.395, and the
    # walk hopped between two of them at `beta = 2.97` -- caught, as on F4,
    # only by the end-to-end disjointness check (`worst separation
    # -7.5e-44`, i.e. two tracked walks on ONE zero).  Anchored at `-0.98`
    # the same families are separated by 0.13 in `s` and the walk holds.
    "b0": 2.90, "b1": 3.10, "nb": 20,        # beta step 0.01
    "y0": -0.98, "y1": -0.90, "ny": 8,       # y step 0.01
    # located at the anchor teacher (2.90, -0.98) and declared here.
    "seeds": [[2.886921, 0.030401],
              [2.895870, 0.230870],
              [2.896822, 0.362417],
              [0.0, 2.90]],                  # the exact fit, th = (beta, 0)
    "r0": 1e-3, "rmin": 1e-11, "min_split": 64,
    "trust": 0.04, "delta": 0.002, "maxjump": 0.006,
    "refine": 80, "newton_iters": 40, "newton_tol_exp": 40,
    "prec": 160,
}

# A 5 x 3 corner of the same rectangle, for the emit/replay round trip: it
# exercises the anchor, the spine, the teeth, both verdict codes and the
# disjointness check in a couple of seconds.
SPECS["F4smoke"] = dict(SPECS["F4"], face="F4 (smoke corner)",
                        b1=0.50, nb=8, y1=-0.20, ny=2)
SPECS["F1smoke"] = dict(SPECS["F1a"], face="F1 (smoke corner)",
                        b1=0.50, nb=4, y1=0.20, ny=2)


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "F4"
    spec = json.loads(json.dumps(SPECS[name]))
    doc, stats = emit_certificate(spec, name)
    out = os.path.join(HERE, "cont_%s.json" % name)
    with open(out, "w") as fh:
        json.dump(doc, fh, indent=1, sort_keys=True)
    e = doc["emit"]
    print("emitted     : %s" % os.path.basename(out))
    print("grid        : beta %d x y %d = %d teachers, %d tracked families"
          % (spec["nb"] + 1, spec["ny"] + 1,
             (spec["nb"] + 1) * (spec["ny"] + 1), len(spec["seeds"])))
    print("steps       : %d certified, %d parameter splits, %d undecided"
          % (e["certified_steps"], e["parameter_splits"], e["undecided"]))
    print("verdicts    : %s" % e["verdicts"])
    print("disjoint    : %s (worst separation %s at %s)"
          % (e["disjoint"], e["worst_separation"], e.get("worst_where")))
    if e.get("clash") is not None:
        c, bg, yg = e["clash"], grid(spec, "b"), grid(spec, "y")
        print("   OVERLAP   : families %d and %d at grid teacher (%d, %d) = "
              "(beta %.6g, y %.6g)"
              % (c["family_a"], c["family_b"], c["grid_i"], c["grid_j"],
                 float(bg[c["grid_i"]]), float(yg[c["grid_j"]])))
    print("sheet guard : maxjump %g, worst step displacement %.4g"
          % (spec["maxjump"], e["max_step_displacement"]))
    print("wall clock  : %.1f s  (%.1f ms per certified step)"
          % (e["seconds"], e["ms_per_step"]))
    print("nodes       : %d bits, %d bytes packed"
          % (doc["n_nodes"], len(base64.b64decode(doc["bits_b64"]))))
    for w in stats["undec_where"][:10]:
        print("   UNDECIDED %s" % (w,))
    for w in stats["bad_where"][:10]:
        print("   FAILED    %s" % (w,))
    return 0 if not stats["undec"] and not stats["bad"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
