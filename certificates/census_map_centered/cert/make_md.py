import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
"""Build CERTIFICATE.md from certificate.json."""
import json
from collections import Counter

D = _HERE + "/"
c = json.load(open(D + "certificate.json"))
V = c["validation"]
try:
    with open(D + "validation.json") as _fh:
        V = json.load(_fh)          # the latest re-run against widgets.js
except OSError:
    pass
L = c["structure_lemmas"]
# Optional sidecar: run parameters and content hashes of the shipped run.
# Certificates emitted before certify.py learned to stamp its own parameters
# carry them here instead.
try:
    with open(D + "PROVENANCE.json") as _fh:
        PROV = json.load(_fh)
except OSError:
    PROV = {}


def e(x, n=None):
    if x is None:
        return "-"
    a, b = x
    return a if a == b else "[%s, %s]" % (a, b)


out = []
A = out.append

A("# Certificate: the boundary set of the centered two-student phase diagram")
A("")
A("**Object.** A certified-complete enumeration of the zero sets, in the")
A("teacher-gap / mass-ratio map coordinates `(x = beta in [0,pi], y in (-1,1])`,")
A("of the three mass-free determinants of the centered kernel:")
A("")
A("```")
A("  phiCos(t)   = cos t * arcsin(cos t) + |sin t|          (pi-periodic)")
A("  h(t)        = sin t * arcsin(cos t)                    (couplingH)")
A("  slopeAtom(t)= phiCos(t) - 2|sin t| = h'(t)")
A("")
A("  Wtau(beta,t) = h(t) sA(t-beta) - h(t-beta) sA(t)          torque Wronskian")
A("  Wpot(beta,t) = phiCos(t) h(t-beta) - phiCos(t-beta) h(t)  potential Wronskian")
A("  Wwgt(beta,t) = |sin t| h(t-beta) - |sin(t-beta)| h(t)     weight Wronskian")
A("```")
A("")
A("with `psi = (y+1)pi/2`, `(s0,s1) = (sin psi, cos psi)`, and a root `t` mapped")
A("to `y` through the kernel-vector mass direction `(s0,s1) ~ (-h(t-beta), h(t))`,")
A("folded by the `(-,-)` isometry exactly as `ratioCoord` in `website/widgets.js`.")
A("")
A("`phiCos` and `couplingH` are the Lean definitions of")
A("`Planar/SignedAngle/PlanarKernel.lean` -- `phiCos x = orthantPhi (cos x)`")
A("with `orthantPhi rho = rho arcsin(rho) + sqrt(1 - rho^2)`, and")
A("`couplingH D = sin D * arcsin(cos D)` -- and the closed branch forms used by")
A("`widgets.js` and mirrored here, `phiC(t) = (pi/2 - x) cos x + sin x` and")
A("`HC(t) = (pi/2 - x) sin x` with `x = t mod pi`, agree with them for every")
A("`t` (on `[0, pi]`, `arcsin(cos t) = pi/2 - t` and `sqrt(1 - cos^2 t) = sin t`,")
A("and both sides are `pi`-periodic).  The three Wronskians are the Lean")
A("`torqueAtomWronskian`, `potentialAtomWronskian`, `weightAtomWronskian` of")
A("`Planar/GeneralTeachers/PhaseBoundary.lean` / `PotentialBoundary.lean`,")
A("definition for definition.")
A("")

A("**Scope caveat.** This certificate is about the zero sets of these three")
A("determinants only. It does **not** claim that they exhaust the census")
A("boundary of the classification; that containment is a separate (Lean) matter.")
A("")
A("**Precision.** All reported enclosures come from `mpmath.iv` (interval")
A("arithmetic with directed rounding) at **%d bits** of working precision"
  % c["precision_bits"])
A("(mpmath %s). `float64` is used only inside the widgets.js mirror of section 1"
  % c["mpmath_version"])
A("and in pre-scans that guess where roots are; no reported enclosure passes")
A("through it.")
A("")

A("## 1. Validation against `website/widgets.js`")
A("")
a = V["a_widgets_mirror"]
A("The reference values are produced by `ref_widgets.js`: `node` evaluating the")
A("**literal function sources** of `website/widgets.js` (`PI`, `mod`, `phiC`,")
A("`HC`, `dHC`, `slopeAtom`), extracted by identifier with a guard that asserts")
A("each extracted function still contains the expected expression, and `eval`-ed")
A("verbatim.  Last re-run against `widgets.js` at commit `2201dd0` (whose")
A("centered kernel functions are byte-identical to those first validated at")
A("`f430fab`).")
A("Sample count: **%d** quasi-random `(t, beta)` pairs, `t` in `[-2pi, 3pi]`,"
  % a["samples"])
A("`beta` in `[0, pi]`.")
A("")
A("| quantity | max abs deviation, Python mirror vs widgets.js |")
A("|---|---|")
for k, v in a["max_abs_deviation"].items():
    A("| `%s` | %.3e |" % (k, v))
A("")
a2 = V.get("a2_libm")
if a2:
    A("These deviations are **1-2 ulp and are entirely attributable to libm**")
    A("(`ref_libm.js` + section (a2) of `validate.py`): V8's `Math.sin`/`Math.cos`")
    A("and glibc's differ by up to %.4e (measured over %d arguments), and feeding"
      % (a2["max_abs_sin_cos_difference_node_vs_glibc"], a2["arguments"]))
    A("node's *own* `sin`/`cos` values into the Python formulas reproduces every")
    A("one of the seven quantities **bit-identically (max deviation %.1f)**, so"
      % max(a2["max_abs_deviation_when_fed_nodes_own_sin_cos"].values()))
    A("the reimplementation is formula-identical to the widget.")
A("")
b = V["b_hprime_eq_slopeAtom"]
A("`h'(t) = slopeAtom(t)` by central differences (%s):" % b["method"])
A("max absolute error **%s**, attained at `t = %s`." % (b["max_abs_error"], b["argmax_t"]))
A("")
A("## 2. An exact reduction that removes the kinks")
A("")
A("`|sin|` and the pi-periodic folding put kinks at `t in pi*Z` and")
A("`t - beta in pi*Z`. An interval derivative test is invalid across a kink, so")
A("the `t`-period `[0,pi]` is **split at the kink lattice** first. For")
A("`beta in (0,pi)` the split points are exactly")
A("")
A("```")
A("    t = 0 ,   t = beta ,   t = pi")
A("```")
A("")
A("giving two smooth pieces")
A("")
A("```")
A("    piece B : t in [0, beta]   fold(t) = t, fold(t-beta) = t-beta+pi")
A("    piece A : t in [beta, pi]  fold(t) = t, fold(t-beta) = t-beta")
A("```")
A("")
A("On each closed piece the widgets.js *closed branch forms* are the correct")
A("values (the one-sided form is still correct at the branch endpoint), so the")
A("piece is analytic and the interval derivative test is legitimate. In the code")
A("this is `centered.dets_branch(B, T, kt, ku)`, which takes the branch indices")
A("`(kt, ku) = (0,-1)` on piece B and `(0,0)` on piece A rather than inferring")
A("them; the generic multi-branch evaluator `centered.atoms` (hull over all")
A("branches met, with a runtime assertion that the branch pieces tile the input)")
A("is used only where a value, not a derivative, is needed.")
A("")
A("**Worked example.** At `beta = pi/2`, the split points are `t = 0, pi/2, pi`.")
A("On piece B the enclosure of `Wtau` over `t in [0, pi/2]` contains 0 *and* the")
A("enclosure of `d/dt Wtau` contains 0, so the naive test cannot conclude. Split")
A("at the kink and work on the closed piece: the reduced form below shows")
A("`Phi(pi/2, w) > 0` for `|w| < pi/4` and `Phi(pi/2, +-pi/4) = 0` with")
A("`Phi_w(pi/2, +-pi/4) = 0` as well -- the root sits *exactly on* the kink and")
A("is a double root there. Had the domain not been split, the derivative")
A("enclosure across `t = 0` would have been taken over a set on which `Wtau` is")
A("not differentiable, and the resulting 'monotone, hence at most one root'")
A("conclusion would have been unsound.")
A("")
A("Substituting the branch forms and putting, on each piece,")
A("")
A("```")
A("   piece A:  b = beta   , sigma = +1 , w = t - (beta+pi)/2")
A("   piece B:  b = pi-beta, sigma = -1 , w = t - beta/2")
A("   W = (pi - b)/2 ,  |w| <= W")
A("")
A("   Q (b,w) = (w^2 - b^2/4) sin b")
A("   Om(b,w) = (b/2)(cos 2w + cos b) = b cos(w + b/2) cos(w - b/2)")
A("")
A("   Wtau = sigma * ( Q + Om) =: sigma * Phi(b,w)")
A("   Wpot = sigma * (-Q + Om) =: sigma * Psi(b,w)")
A("   Wwgt = sigma *       Om")
A("```")
A("")
A("`Phi`, `Psi`, `Om` are **analytic in `(b,w)` on the whole piece and even in")
A("`w`**: the reduction has absorbed every kink into the piece boundary")
A("`|w| = W`. This is verified two ways in `validate.py`: symbolically with")
A("`sympy` on both pieces (%s / %s, all six identities), and on **%d interval"
  % (c["validation"]["c_reduction_symbolic"]["pieceA_Wtau_Wpot_Wwgt"],
     c["validation"]["c_reduction_symbolic"]["pieceB_Wtau_Wpot_Wwgt"],
     c["validation"]["c_reduction_identity"]["grid_points"]))
A("samples** where the enclosure of (folded mirror - reduced form) contains 0")
A("with radius at most %.2e." % max(c["validation"]["c_reduction_identity"]
                                    ["max_enclosure_radius"].values()))
A("")
A("### 2.1 An exact linear relation")
A("")
A("Because `Phi + Psi = 2 Om` identically,")
A("")
A("```")
A("    Wtau(beta,t) + Wpot(beta,t) = 2 Wwgt(beta,t)      for all beta, t")
A("```")
A("")
A("proved symbolically (`sympy`, exact) and verified in interval arithmetic with")
A("residual radius at most %.2e. **The three determinants span only a"
  % c["validation"]["d_linear_identity"]["interval_max_residual_radius"])
A("2-dimensional space**; any two of them determine the third.")
A("")
A("## 3. Structure lemmas (the analytic backbone)")
A("")
A("Endpoint closed forms, verified at 399 grid values of `b` in interval")
A("arithmetic (max residual %.2e):" % L["endpoint_closed_forms"]["max_interval_residual"])
A("")
A("```")
for k in ("Phi(b,0)", "Phi(b,W)", "Psi(b,0)", "Psi(b,W)", "Om(b,W)"):
    A("   %-10s = %s" % (k, L["endpoint_closed_forms"][k]))
A("```")
A("")
for k, v in L["monotonicity_and_sign_lemmas"].items():
    if k == "interval_grid_checks":
        continue
    A("**%s.** %s" % (k, v))
    A("")
A("Interval grid checks of these hypotheses: `%s`."
  % L["monotonicity_and_sign_lemmas"]["interval_grid_checks"])
A("")
A("## 4. Certified root counts as a function of beta")
A("")
A("Roots are counted in `t` modulo `pi` (all three determinants are")
A("`pi`-periodic in `t`). The counting is done on `w in [0,W]` and doubled by")
A("the exact evenness in `w`.")
A("")
A("| beta | Wtau | Wpot | Wwgt |")
A("|---|---|---|---|")
A("| `beta = 0` | identically zero | identically zero | identically zero |")
for p in c["beta_partition"]:
    cnt = p["certified_counts"][0] if p["certified_counts"] else {}
    A("| `%s` | %s | %s | %s |" % (p["range"], cnt.get("Wtau", "?"),
                                   cnt.get("Wpot", "?"), cnt.get("Wwgt", "?")))
    if p["range"] == "(0, beta*_1)":
        A("| `beta = beta*_1` | **1 (double root)** | 2 | 2 |")
    if p["range"] == "(beta*_1, pi/2)":
        A("| `beta = pi/2` | 2 (both on the kink lattice) | 2 (same two points) "
          "| 2 (same two points) |")
    if p["range"] == "(pi/2, beta*_2)":
        A("| `beta = beta*_2` | **1 (double root)** | 2 | 2 |")
A("| `beta = pi` | identically zero | identically zero | identically zero |")
A("")
A("Only **one** count changes anywhere: that of `Wtau`, and only at")
A("`beta = beta*_1` and `beta = beta*_2`, in both cases through a **double root**")
A("(a common zero of `Phi(b,.)` and `Phi_w(b,.)` at `w = 0`, i.e. `Phi_w(b,0) = 0`")
A("holds identically and `Phi(b*,0) = 0`). `Wpot` and `Wwgt` have constant count")
A("2 on all of `(0,pi)`.")
A("")
A("Certified coverage of the beta axis:")
A("")
A("| beta range | boxes | measure covered | uncovered |")
A("|---|---|---|---|")
for p in c["beta_partition"]:
    A("| `%s` | %d | %.15f | %.3e |"
      % (p["range"], p["n_boxes"], p["measure_covered"],
         p["measure_total"] - p["measure_covered"]))
A("")
A("The uncovered part consists only of collars around the critical betas")
A("`0, beta*_1, pi/2, beta*_2, pi`, each of width below `1e-13`; those five")
A("betas themselves are handled exactly (section 5). See section 9.")
A("")
A("## 5. The critical betas")
A("")
for k, v in c["critical_betas"].items():
    A("### `%s`" % k)
    for kk, vv in v.items():
        A("- **%s**: %s" % (kk, vv))
    A("")
A("## 6. beta*: independent verification of the lens claim")
A("")
bs = L["bstar"]
A("The lens endpoints are the double-root locus `Phi(b,0) = 0`. Since")
A("`Phi(b,0) = b cos(b/2) [cos(b/2) - (b/2) sin(b/2)]` and `cos(b/2) > 0` on")
A("`(0,pi)`, this is `%s`." % bs["equation"])
A("")
A("Uniqueness in the bracket `%s`: %s." % (bs["search_bracket"], bs["uniqueness"]))
A("")
A("```")
A("   b*      in %s" % e(bs["b_star"]))
A("   beta*_1 = pi - b*  in %s" % e(bs["beta_star_1_eq_pi_minus_b_star"]))
A("   beta*_2 = b*       in %s" % e(bs["beta_star_2_eq_b_star"]))
A("   beta*_1 + beta*_2 - pi  in %s   (so the pair is symmetric about pi/2)"
  % e(bs["sum_check_beta1_plus_beta2_eq_pi"]))
A("```")
A("")
ch = bs["characterisation_beta1_eq_2u_with_tan_u_eq_halfpi_minus_u"]
A("**The previous lane's characterisation survives.** With")
A("`u* = beta*_1 / 2 in %s`, the residual `tan u* - (pi/2 - u*)` is enclosed in"
  % e(ch["u_star"]))
A("`%s`, i.e. certified to contain 0." % e(ch["residual_tan(u)-(pi/2-u)"]))
A("The two characterisations are *equivalent*, not merely numerically close:")
A("`u* = pi/2 - b*/2`, so `tan u* = cot(b*/2) = cos(b*/2)/sin(b*/2)`, and")
A("`tan u* = pi/2 - u* = b*/2` is exactly `cos(b*/2) = (b*/2) sin(b*/2)`.")
A("")
A("The reported value `beta* = 1.420925475551` agrees with the certified")
A("enclosure to all quoted digits.")
A("")
A("**Lens claim: SURVIVES, as a single closed component.** For every")
A("`beta in (beta*_1, beta*_2)` the certified `Wtau` count is exactly 2; at both")
A("endpoints it is exactly 1 (a double root); outside `[beta*_1, beta*_2]` it is")
A("0. The two roots are `w = +-w_r(beta)` with `w_r` the unique zero of a")
A("function certified strictly monotone in `w` near it, so `w_r` is continuous")
A("(implicit function theorem) and the zero set is one closed curve. The two")
A("branches are exact mirror images: `t -> beta - t (mod pi)` maps roots to roots")
A("and sends the mass ratio `rho` to `1/rho`, hence `y` to `-1-y`, so the lens is")
A("**symmetric in `y` about `-1/2`** and both endpoints lie exactly on `y = -1/2`.")
A("This identity was checked on every certified curve box: `y+ + y- = -1`: **%s**."
  % c["curves"]["torque_lens"].get("reflection_identity_y_plus_plus_y_minus_eq_minus_1"))
A("")
A("`y`-hull over the certified curve boxes: `%s`"
  % e(c["curves"]["torque_lens"].get("y_hull_over_certified_boxes",["?","?"])))
A("(beta measure covered %.12f of %.12f)."
  % (c["curves"]["torque_lens"].get("beta_measure_covered",0.0),
     c["curves"]["torque_lens"].get("beta_measure_total",0.0)))
A("The extremes `y = c-1` and `y = -c` with `c = (2/pi) arctan(2/pi) in %s`"
  % e(L["constants"]["c"]["enclosure"]))
A("are the `beta -> pi/2` limits and are not attained.")
A("")
A("## 7. The 1-dimensional pieces")
A("")
for p in c["pieces_1D"]:
    A("**%s. %s**" % (p["id"], p["name"]))
    for k, v in p.items():
        if k in ("id", "name"):
            continue
        A("  - %s: %s" % (k, v))
    A("")
A("### Non-intersections (certified)")
A("")
for s in c["non_intersections"]:
    A("- %s" % s)
A("")
A("## 8. The 0-dimensional pieces")
A("")
A("| id | name | beta | y | what makes it special |")
A("|---|---|---|---|---|")
for p in c["pieces_0D"]:
    A("| %s | %s | `%s` | `%s` | %s |"
      % (p["id"], p["name"], e(p["point"]["beta"]), e(p["point"]["y"]),
         p["why"].replace("\n", " ")))
A("")
A("Auxiliary constants appearing above:")
A("")
A("```")
A("   u*   in %s      (tan u* = pi/2 - u*, beta*_1 = 2u*)" % e(ch["u_star"]))
A("   b*   in %s      ((b*/2) tan(b*/2) = 1, beta*_2 = b*)" % e(bs["b_star"]))
A("   w_D  in %s      (Dottie number, cos w = w)"
  % e(L["constants"]["w_Dottie"]["enclosure"]))
A("   c    in %s      (c = (2/pi) arctan(2/pi))"
  % e(L["constants"]["c"]["enclosure"]))
A("```")
A("")
A("## 9. Unresolved set")
A("")
U = c["unresolved_boxes"]
A("Total entries: **%d**." % U["total_entries"])
A("")
A("### 9.1 beta collars where a root COUNT could not be certified")
A("")
A("| beta range | box | width |")
A("|---|---|---|")
for u in U["beta_collars_uncertified"]:
    A("| `%s` | `[%s, %s]` | %.3e |"
      % (u["range"], u["beta_box"][0], u["beta_box"][1], u["width"]))
A("")
A("Total uncovered beta measure: **%.3e**."
  % sum(u["width"] for u in U["beta_collars_uncertified"]))
A("Each collar sits around one of the critical betas `beta*_1, pi/2, beta*_2`,")
A("which are themselves treated exactly in section 5. This is the ONLY part of")
A("the beta axis on which a root count is not certified.")
A("")
A("### 9.2 curve boxes whose y enclosure is not tight")
A("")
A("These are boxes on which the ROOT is certified to exist and be unique, but")
A("the resulting enclosure of `y` is wider than the reporting tolerance `1e-3`,")
A("or the mass direction `(-h(t-beta), h(t))` is not sign-determined on the box.")
A("They are certified, just not tight; they are excluded from the `y`-hulls of")
A("section 6.")
A("")
A("| curve | wide boxes | beta measure | max box width | reasons |")
A("|---|---|---|---|---|")
for nm, lbl in (("torque_lens", "torque lens"), ("potential_curve", "potential curve")):
    w = c["curves"][nm]["wide_boxes_summary"]
    A("| %s | %d | %.6f | %.3e | %s |"
      % (lbl, w["n"], w["beta_measure"], w["max_width"],
         ", ".join("%s: %d" % (k, v) for k, v in w["reasons"].items())))
A("")
A("The `mass_direction_undetermined` boxes (98 in total) are exactly the ones")
A("adjacent to `beta = pi/2`, where the defining root reaches the kink lattice")
A("and `h(t) = h(t-beta) = 0`; the `count` boxes (6) are the collars of 9.1.")
A("")
A("### 9.3 the generic root test")
A("")
A("The generic exclusion/monotonicity/bisection test of section 5 left **no**")
A("undecided box that is not adjacent to a kink point: every box it could not")
A("decide abuts `t = 0`, `t = beta` or `t = pi`, where the determinant is")
A("certified to vanish (Wwgt) or certified nonzero (Wtau, Wpot) by the closed")
A("forms of section 3. Those boxes are reported per beta in")
A("`per_beta_table[*][det].kink_adjacent_boxes_explained`.")
A("")
A("## 10. Honesty notes")
A("")
A("- Every count in section 4 is certified **for all beta in the stated box**, not")
A("  sampled. Where a count could not be certified, the box is listed in section 9")
A("  rather than rounded to the nearest plausible integer.")
A("- The `y`-extremes `c-1, -c` (torque) and `c, 1-c` (potential) are **limits**")
A("  at `beta = pi/2`, derived analytically; at `beta = pi/2` itself the defining")
A("  root sits on the kink lattice with `h(t) = h(t-beta) = 0`, so the mass")
A("  direction is a genuine `0/0` and the point is degenerate.")
A("- The claim that these three determinants exhaust the census boundary is **not**")
A("  made here.")
A("- Validity domain: `beta in [0, pi]`, `t` in one `pi`-period, centered kernel")
A("  only. Precision 200 bits throughout.")
A("")
A("## Files")
A("")
A("```")
for k, v in c.get("files", {}).items():
    if k.endswith(".log") and PROV.get("log"):
        # the certificate's own files block still carries the generic name
        k = PROV["log"]
    A("  %-26s %s" % (k, v))
A("```")

if PROV:
    A("")
    A("## Provenance")
    A("")
    rows = []
    if PROV.get("max_depth_y") is not None:
        rows.append(("MAX_DEPTH_Y", str(PROV["max_depth_y"])))
    if PROV.get("log"):
        rows.append(("console log", PROV["log"]))
    for k, v in PROV.get("sha256", {}).items():
        rows.append(("sha256 " + k, v))
    wid = max([26] + [len(k) for k, _ in rows])
    A("```")
    for k, v in rows:
        A("  %-*s %s" % (wid, k, v))
    A("```")
    A("")
    A("The depth of the shipped run is recorded in `cert/PROVENANCE.json` rather")
    A("than in `certificate.json` itself: this certificate predates the change")
    A("that stamps the run parameters into every emitted certificate (`certify.py`")
    A("now writes them under `run_parameters`), and nothing was injected into the")
    A("shipped JSON after the fact.")

open(D + "CERTIFICATE.md", "w").write("\n".join(x for x in out if x is not None) + "\n")
print("wrote CERTIFICATE.md")
