"""
Compare the certified arrangement with the SHIPPED mosaic (website/census-map.js,
depth 8, and the depth 6/7 reruns), and explain every discrepancy.

Both discrepancies below have since been fixed in website/ -- the selector's
absolute cutoff was made relative, and the generator now takes the potential
curve as a separator -- so the shipped mosaic reports 4 census keys in 5
finite-grid flood-fill groups, and its representatives agree with the certified
face signatures where they land.  That renderer output is not a component
certificate.  What follows is the measurement that identified the defects; it
is kept because the two mechanisms are what the fixes had to address.  See
ARRANGEMENT.md, "SUPERSEDED -- the mosaic comparison".

Two things had to be measured, not asserted:

(1) the four 514-cell "fit:global" slivers at the corners of the upper half.
    They are NOT faces of the arrangement: no certified curve passes within
    0.35 of them in y.  They are the region where widgets.js
    `coincidenceTypeAt` returns "topological saddle" because the ABSOLUTE cutoff
    SEL_EPS = 1e-7 fires on the factor W = s0|sin t| + s1|sin(t-beta)| evaluated
    at the coincidence root.  W is proportional to the weight Wronskian Wwgt,
    which vanishes exactly ON the strata y = 0 and y = +-1, so the cutoff
    thickens those measure-zero strata into a band; near beta = 0 and beta = pi
    the band is wide enough to be resolved by the depth-8 grid, elsewhere it is
    not.  This script measures the band's width as a function of beta and checks
    it against the sliver outlines actually stored in census-map.js.

(2) the potential lens pinches to the single point (0, 1/2) at beta = 0 and to
    (pi, 1/2) at beta = pi.  In the arrangement the region BELOW the lower
    branch and the region ABOVE the upper branch are two different faces
    (F5/F8 and F7/F10); the mosaic coalesces samples from both bands into one
    finite-grid flood-fill group because at its finest column the lens is
    thinner than one cell.  This script computes the certified lens width at
    the mosaic's first column and compares.
"""
import json
import math
import subprocess
import mpmath
import cert_core as K

C = K.C
iv = K.iv

NX, NY = 6144, 3072
DBETA = math.pi / NX
DY = 2.0 / NY


def lens_width(beta):
    """certified y-width of the potential lens at beta."""
    bs = K.curve_branches(beta, "potential")
    ys = sorted(mpmath.mpf(b["y_lo"]) for b in bs if b["y_lo"] is not None)
    if len(ys) != 2:
        return None
    hi = sorted(mpmath.mpf(b["y_hi"]) for b in bs if b["y_hi"] is not None)
    return float(ys[1] - hi[0]), [mpmath.nstr(ys[0], 20), mpmath.nstr(ys[1], 20)]


def main():
    out = {}

    # ---- (2) the pinch ---------------------------------------------------
    print("potential lens width vs beta (mosaic cell height %.6e):" % DY)
    rows = []
    for k in [0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 64.5, 256.5, 1024.5, 3072.5]:
        beta = float(k * DBETA)
        w, ys = lens_width(beta)
        rows.append({"column": k, "beta": beta, "lens_y_width": w,
                     "cells": w / DY, "branch_y": ys})
        print("   column %8.1f  beta=%.8e  width=%.6e  = %.3f mosaic cells"
              % (k, beta, w, w / DY), flush=True)
    out["potential_lens_pinch"] = {
        "mosaic_cell_height_dy": DY, "mosaic_cell_width_dbeta": DBETA,
        "rows": rows,
        "conclusion":
            "the two potential-curve branches are distinct for every beta in "
            "(0, pi) and pinch only AT beta = 0 and beta = pi; the historical "
            "mosaic coalesced samples from the two outer POS face bands because "
            "at its first columns the lens was thinner than one cell.  This is "
            "a finite-grid observation, not a connected-component claim."}

    # ---- (1) the SEL_EPS band -------------------------------------------
    # W(t_root) = -kappa * Wwgt(beta, t_root); the sliver is {|W| < 1e-7}.
    js = r'''
const api = require("./harness.js").load();
const PI=Math.PI;
function mod(x,m){var r=x%m;return r<0?r+m:r;}
function phiC(t){var x=mod(t,PI);return (PI/2-x)*Math.cos(x)+Math.sin(x);}
function HC(t){var x=mod(t,PI);return (PI/2-x)*Math.sin(x);}
function P(T,t){return T.s0*phiC(t)+T.s1*phiC(t-T.beta);}
function Wg(T,t){return T.s0*Math.abs(Math.sin(t))+T.s1*Math.abs(Math.sin(t-T.beta));}
const SEL=1e-7;
// for a given beta, bisect in y for the edge of the band {min|W| over the
// coincidence roots that would otherwise be a trap} = SEL
function minWtrap(beta,y){
  const m=api.massesAt(y), T={beta:beta,s0:m.s0,s1:m.s1};
  let best=Infinity;
  for(const r of api.classificationRows("centered",T)){
    if(r.family!=="coincidence root") continue;
    const t=r.theta[0], p=P(T,t), w=Wg(T,t), tau=p-2*w;
    if(tau*p>0) best=Math.min(best,Math.abs(w));
  }
  return best;
}
const out=[];
for(const beta of [0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.3,0.8,1.4,2.4,3.0,3.10,3.13,3.139,3.1406,3.1411]){
  // upper edge of the band above y=0
  let lo=1e-12, hi=0.2, f=y=>minWtrap(beta,y)-SEL;
  let r=null;
  if(f(lo)<0 && f(hi)>0){ for(let i=0;i<200;i++){const m=0.5*(lo+hi); if(f(m)<0) lo=m; else hi=m;} r=0.5*(lo+hi); }
  // and below y=+1
  let lo2=1-1e-12, hi2=0.8, g=y=>minWtrap(beta,y)-SEL, r2=null;
  if(g(lo2)<0 && g(hi2)>0){ let a=lo2,b=hi2; for(let i=0;i<200;i++){const m=0.5*(a+b); if(g(m)<0) a=m; else b=m;} r2=0.5*(a+b); }
  out.push({beta:beta, band_above_y0:r, band_below_y1:(r2===null?null:1-r2)});
}
console.log(JSON.stringify(out));
'''
    with open("_band.js", "w") as fh:
        fh.write(js)
    band = json.loads(subprocess.check_output(["node", "_band.js"], text=True))
    print("\nSEL_EPS = 1e-7 band around the strata y = 0 and y = +-1:")
    print("   (mosaic cell height %.4e)" % DY)
    for r in band:
        print("   beta=%9.5f  band above y=0: %s cells   band below y=1: %s cells"
              % (r["beta"],
                 ("%.3f" % (r["band_above_y0"] / DY)) if r["band_above_y0"] else "  -  ",
                 ("%.3f" % (r["band_below_y1"] / DY)) if r["band_below_y1"] else "  -  "))
    out["sel_eps_band"] = {"SEL_EPS": 1e-7, "rows": band,
                           "mosaic_cell_height_dy": DY}

    # ---- the four slivers actually stored in census-map.js ---------------
    js2 = r'''
global.window={};
new Function(require("fs").readFileSync(
  "/home/jzimmermann/code/Skip-Connections-Avoid-Spurious-Local-Minima/website/census-map.js","utf8"))();
const M=window.CensusMap.models.centered, keys=M.keys;
const out=[];
for(const p of M.pieces){
  const xs=[].concat(...p.loops.map(l=>l.map(q=>q[0])));
  const ys=[].concat(...p.loops.map(l=>l.map(q=>q[1])));
  out.push({key:keys[p.id], cells:p.cells,
            beta_range:[Math.min(...xs)*Math.PI/M.nx, Math.max(...xs)*Math.PI/M.nx],
            y_range:[2*Math.min(...ys)/M.ny-1, 2*Math.max(...ys)/M.ny-1],
            rep_beta:p.at[0]*M.per, rep_y:2*p.at[1]-1, loops:p.loops.length});
}
console.log(JSON.stringify(out));
'''
    with open("_pieces.js", "w") as fh:
        fh.write(js2)
    pieces = json.loads(subprocess.check_output(["node", "_pieces.js"], text=True))
    print("\nshipped mosaic pieces (depth 8):")
    for p in pieces:
        print("   %-62s cells=%8d beta in [%.5f,%.5f] y in [%+.5f,%+.5f] loops=%d"
              % (p["key"][:62], p["cells"], p["beta_range"][0], p["beta_range"][1],
                 p["y_range"][0], p["y_range"][1], p["loops"]))
    out["mosaic_pieces"] = pieces

    with open("mosaic_comparison.json", "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
