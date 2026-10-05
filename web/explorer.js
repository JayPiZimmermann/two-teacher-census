// Two-teacher census explorer.
//
// Extracted from widgets.js of the project website
// (github.com/JayPiZimmermann/Skip-Connections-Avoid-Spurious-Local-Minima,
// website/widgets.js, function initLandscapeExplorer) so that the explorer
// runs on its own page.  The mathematics is unchanged: every table row is a
// constructor of the two complete Lean classifications, evaluated from its
// closed forms; type labels are theorem data, never inferred numerically.
// precompute/precompute_census_map.js drives THIS file through a DOM shim to
// generate census-map.js, so the map cannot drift from the table.
(function () {
  "use strict";

  var PI = Math.PI;
  var SVG_NS = "http://www.w3.org/2000/svg";

  function css(name, fallback) {
    var value = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
    return value || fallback;
  }

  function svgEl(tag, attrs, text) {
    var element = document.createElementNS(SVG_NS, tag);
    Object.keys(attrs || {}).forEach(function (key) { element.setAttribute(key, attrs[key]); });
    if (text !== undefined) element.textContent = text;
    return element;
  }

  function clearSvg(svg) {
    while (svg.firstChild) svg.removeChild(svg.firstChild);
  }

  // One thin width for every line that carries data, in every widget here;
  // the frame's own lines -- axes, guide circles, rules, ticks -- go thinner.
  var LINE = 1.5;
  var HAIR = 1;

  function line(svg, x1, y1, x2, y2, color, width, dash) {
    var attrs = {x1:x1, y1:y1, x2:x2, y2:y2, stroke:color, "stroke-width":width || LINE};
    if (dash) attrs["stroke-dasharray"] = dash;
    svg.appendChild(svgEl("line", attrs));
  }

  function textNode(svg, x, y, value, attrs) {
    var merged = {x:x, y:y, fill:css("--ink", "#222"), "font-size":13, "text-anchor":"middle"};
    Object.keys(attrs || {}).forEach(function (key) { merged[key] = attrs[key]; });
    svg.appendChild(svgEl("text", merged, value));
  }

  function circle(svg, x, y, r, fill, stroke) {
    svg.appendChild(svgEl("circle", {cx:x, cy:y, r:r, fill:fill, stroke:stroke || fill, "stroke-width":LINE}));
  }

  // The frame -- crossed hairlines and the outer circle -- may be drawn about a
  // given centre and at a given size, so a zooming panel can hold it still on
  // screen while the data inside magnifies.
  function drawAxes(svg, radius, cx, cy) {
    var rule = css("--rule", "#bbb");
    cx = cx || 0; cy = cy || 0;
    line(svg, cx - radius, cy, cx + radius, cy, rule, HAIR);
    line(svg, cx, cy - radius, cx, cy + radius, rule, HAIR);
    svg.appendChild(svgEl("circle", {cx:cx, cy:cy, r:radius, fill:"none",
      stroke:rule, "stroke-width":HAIR, "vector-effect":"non-scaling-stroke"}));
  }

  function drawArrow(svg, x0, y0, x1, y1, color, width) {
    line(svg, x0, y0, x1, y1, color, width || LINE);
    var angle = Math.atan2(y1-y0, x1-x0);
    var size = 8;
    var points = [
      [x1, y1],
      [x1-size*Math.cos(angle-0.45), y1-size*Math.sin(angle-0.45)],
      [x1-size*Math.cos(angle+0.45), y1-size*Math.sin(angle+0.45)]
    ];
    svg.appendChild(svgEl("polygon", {points:points.map(function (p) { return p.join(","); }).join(" "), fill:color}));
  }

  //
  // The table mirrors the constructors of the two complete Lean classifications.
  // Each row's angles and masses are evaluated from that constructor's
  // closed forms; scalar roots of displayed equations (torque roots, dead
  // residual roots, separated balance roots) are located by bisection or a
  // certified-bracket Newton polish.  Type labels are theorem data: where a
  // family's proved type depends on a displayed sign selector, the widget
  // evaluates that exact selector at the drawn representative.  No type is
  // ever inferred from sampled losses or a numerical Hessian.
  // ---------------------------------------------------------------------

  function initLandscapeExplorer() {
    var canvas = document.getElementById("ex-canvas");
    if (!canvas) return;

    var EPS = 1e-9;

    function mod(x, m) { var r = x % m; return r < 0 ? r + m : r; }

    // Kernels: closed branch forms of phiCos / phiCosJ and their couplings.
    function phiC(t) { var x = mod(t, PI); return (PI / 2 - x) * Math.cos(x) + Math.sin(x); }
    function HC(t)   { var x = mod(t, PI); return (PI / 2 - x) * Math.sin(x); }
    function dHC(t)  { var x = mod(t, PI); return (PI / 2 - x) * Math.cos(x) - Math.sin(x); }
    function phiJ(t) { var x = mod(t, 2 * PI); if (x > PI) x = 2 * PI - x; return (PI - x) * Math.cos(x) + Math.sin(x); }
    function HJ(t)   { var x = mod(t, 2 * PI); return x <= PI ? (PI - x) * Math.sin(x) : (x - PI) * Math.sin(x); }
    function dHJ(t)  { var x = mod(t, 2 * PI); return x <= PI ? (PI - x) * Math.cos(x) - Math.sin(x) : Math.sin(x) + (x - PI) * Math.cos(x); }

    function kernelOf(model) {
      return model === "centered"
        ? {phi: phiC, H: HC, dH: dHC, kap: PI / 2, per: PI}
        : {phi: phiJ, H: HJ, dH: dHJ, kap: PI, per: 2 * PI};
    }
    function Pot(K, T, t) { return T.s0 * K.phi(t) + T.s1 * K.phi(t - T.beta); }
    function Trq(K, T, t) { return -(T.s0 * K.H(t) + T.s1 * K.H(t - T.beta)); }
    function Wgt(T, t) { return T.s0 * Math.abs(Math.sin(t)) + T.s1 * Math.abs(Math.sin(t - T.beta)); }
    function Tau(K, T, t) { return Pot(K, T, t) - 2 * Wgt(T, t); }
    // The torque slope A' = 2W - P, i.e. -Tau: both kernels obey K'' + K = 2|sin|,
    // so a torque zero is a potential MAXIMUM exactly when Tau > 0 there.
    function TorqueSlope(K, T, t) { return -Tau(K, T, t); }

    function findRoots(f, lo, hi, n) {
      n = n || 1440;
      var roots = [], prevX = lo, prevV = f(lo);
      for (var i = 1; i <= n; i += 1) {
        var x = lo + (hi - lo) * i / n, v = f(x);
        if (isFinite(prevV) && isFinite(v) && prevV * v < 0) {
          var a = prevX, b = x;
          for (var k = 0; k < 80; k += 1) {
            var m = 0.5 * (a + b);
            if (f(a) * f(m) <= 0) b = m; else a = m;
          }
          var r = 0.5 * (a + b);
          if (Math.abs(f(r)) < 1e-7 &&
              !roots.some(function (q) { return Math.abs(q - r) < 1e-5; })) roots.push(r);
        }
        prevX = x; prevV = v;
      }
      return roots;
    }

    // -------------------------------------------------------------------
    // The census map.
    //
    // Horizontal: the teacher gap over the kernel's own period (pi centered,
    // 2pi noncentered) -- one fixed box for both, so the noncentered map is
    // the compressed one.  Vertical: the mass ratio through
    //     y = (2/pi)*atan(s0/s1),   (s0, s1) = (sin(y*pi/2), cos(y*pi/2)),
    // which is bounded, so ratio -> +-infinity is the TOP and BOTTOM EDGE of
    // the box rather than a limit off the plot.  That coordinate is also the
    // right quotient: it is unchanged by negating BOTH masses, which is an
    // exact isometry of the loss, so the (+,+) and (-,-) teachers are one
    // point of the map and only two sign sectors remain -- same sign above
    // y = 0, mixed sign below it.
    //
    // A cell's shade is its CENSUS: how many torque zeros the teacher
    // potential has and how many of them are maxima.  That is what the
    // classification's family list is indexed by, so a shade boundary is
    // exactly where the census changes.  The boundary curve itself is the
    // mass-free equation proved in
    // LeanFormalization/Planar/GeneralTeachers/PhaseBoundary.lean.
    var MAP_W = 168, MAP_H = 120;
    var curveCache = {};

    // The slope atom h'(y) = K(y) - 2|sin y| (hasDerivAt_couplingH), so the
    // Wronskian below is the mass-free boundary equation of
    // LeanFormalization/Planar/GeneralTeachers/PhaseBoundary.lean.
    function slopeAtom(K, y) { return K.phi(y) - 2 * Math.abs(Math.sin(y)); }
    function wronskian(K, beta, t) {
      return K.H(t) * slopeAtom(K, t - beta) - K.H(t - beta) * slopeAtom(K, t);
    }

    // The region boundaries, drawn from that equation rather than traced around
    // the shading grid: for each gap the degenerate torque roots are found in t,
    // and the mass ratio each one forces is read off as (-h(t-beta), h(t)).
    // The result is a smooth curve instead of a staircase of cell edges.
    function boundaryCurves(model) {
      if (curveCache[model]) return curveCache[model];
      var K = kernelOf(model), per = K.per, COLS = 260;
      var columns = [];
      for (var c = 0; c <= COLS; c += 1) {
        var beta = per * c / COLS;
        var ys = [];
        findRoots(function (t) { return wronskian(K, beta, t); }, 1e-6, per - 1e-6, 900)
          .forEach(function (t) {
            var a = K.H(t), b = K.H(t - beta);
            if (Math.abs(a) < 1e-9) return;
            var rho = -b / a;
            if (rho > 1e-9) ys.push(ratioCoord(-b, a));
          });
        ys.sort(function (p, q) { return p - q; });
        columns.push({beta: beta, ys: ys});
      }
      // Chain points across neighbouring columns into polylines by nearest y.
      var curves = [], open = [];
      columns.forEach(function (col) {
        var used = col.ys.map(function () { return false; });
        open = open.filter(function (poly) {
          var last = poly[poly.length - 1], best = -1, bestd = 0.09;
          col.ys.forEach(function (y, i) {
            var d = Math.abs(y - last.y);
            if (!used[i] && d < bestd) { bestd = d; best = i; }
          });
          if (best < 0) { curves.push(poly); return false; }
          used[best] = true; poly.push({beta: col.beta, y: col.ys[best]});
          return true;
        });
        col.ys.forEach(function (y, i) {
          if (!used[i]) open.push([{beta: col.beta, y: y}]);
        });
      });
      curves = curves.concat(open).filter(function (p) { return p.length > 2; });
      curveCache[model] = curves;
      return curves;
    }

    // The centered map's other certified separator, Wpot = 0
    // (precompute_census_map.js SEPARATOR_CURVES / certificates/
    // census_map_centered/ARRANGEMENT.md section 2): same shape of equation
    // as the torque lens's wronskian above, with phi in place of the slope
    // atom, so it is traced by the identical column-and-chain method.
    function potentialWronskian(K, beta, t) {
      return K.phi(t) * K.H(t - beta) - K.phi(t - beta) * K.H(t);
    }
    var curveCacheExt = {};
    // Generalisation of the column-and-chain tracer above to any mass-free
    // equation of this shape: locate its t-roots per column, read off the
    // mass ratio each one forces, and chain nearby y across columns into
    // polylines.  `keepAll` skips the opposite-sign filter that keeps only
    // the torque lens's y < 0 branch -- the potential curve needs both signs
    // (it lives in y in (0,1), same-sign masses; certificates/
    // census_map_centered/ARRANGEMENT.md section 2(b)-(c)).
    function tracedCurve(model, key, eqFn, keepAll) {
      var cacheKey = model + ":" + key;
      if (curveCacheExt[cacheKey]) return curveCacheExt[cacheKey];
      var K = kernelOf(model), per = K.per, COLS = 260;
      var columns = [];
      for (var c = 1; c <= COLS; c += 1) {   // c = 0 (beta = 0) is degenerate: both equations vanish identically
        var beta = per * c / COLS;
        var ys = [];
        findRoots(function (t) { return eqFn(K, beta, t); }, 1e-6, per - 1e-6, 900)
          .forEach(function (t) {
            var a = K.H(t), b = K.H(t - beta);
            if (Math.abs(a) < 1e-9) return;
            var rho = -b / a;
            if (keepAll || rho > 1e-9) ys.push(ratioCoord(-b, a));
          });
        ys.sort(function (p, q) { return p - q; });
        columns.push({beta: beta, ys: ys});
      }
      var curves = [], open = [];
      columns.forEach(function (col) {
        var used = col.ys.map(function () { return false; });
        open = open.filter(function (poly) {
          var last = poly[poly.length - 1], best = -1, bestd = 0.09;
          col.ys.forEach(function (y, i) {
            var d = Math.abs(y - last.y);
            if (!used[i] && d < bestd) { bestd = d; best = i; }
          });
          if (best < 0) { curves.push(poly); return false; }
          used[best] = true; poly.push({beta: col.beta, y: col.ys[best]});
          return true;
        });
        col.ys.forEach(function (y, i) {
          if (!used[i]) open.push([{beta: col.beta, y: y}]);
        });
      });
      curves = curves.concat(open).filter(function (p) { return p.length > 2; });
      curveCacheExt[cacheKey] = curves;
      return curves;
    }

    // The two certified curve strata beyond the straight lines of
    // lineStrata: the torque lens (both kernels; the same mass-free
    // boundary boundaryCurves already traces, section 6.1 of both
    // certificates) and, centered only, the potential curve (section 2 of
    // certificates/census_map_centered/ARRANGEMENT.md -- no such curve is
    // certified for the noncentered arrangement).  Clipped to the map's own
    // fundamental half-domain [0, per/2]: the reflection beta -> per - beta
    // fixes every y on both curves (the same symmetry that lets lineStrata's
    // "v" strata stop at F.half), so the far half only repeats them.
    function curveStrata(model) {
      var half = kernelOf(model).per / 2, out = [];
      function clipped(curves) {
        return curves
          .map(function (p) { return p.filter(function (q) { return q.beta <= half + 1e-6; }); })
          .filter(function (p) { return p.length > 2; });
      }
      clipped(tracedCurve(model, "lens", wronskian, false)).forEach(function (p) {
        out.push({name: "torque lens", curve: p});
      });
      if (model === "centered") {
        clipped(tracedCurve(model, "potential", potentialWronskian, true)).forEach(function (p) {
          out.push({name: "potential curve", curve: p});
        });
      }
      return out;
    }

    // The separated stratum's fold walls (certificates/census_map_
    // noncentered/CERTIFICATE.md 6.1 and 6.3): certified only as a BRACKET
    // per beta, not a closed form, so there is no equation to trace.  Drawn
    // instead as the certified data itself -- a shaded band between the
    // bracket's own lo/hi, per-beta -- with a dashed line through the
    // bracket midpoints as the click target.  `wa` is certificate.json's
    // `ywall` (the mixed-sector wall, hugging y = 0 from below); `collar` is
    // CERTIFICATE.md 6.3's table (the same-sign collar wall, hugging y = 0
    // from above).  The mirror rows are NOT separately measured: they follow
    // from the exact symmetries y -> -1-y (mixed sector, CERTIFICATE.md 2.4)
    // and y -> 1-y (same-sign collar, 6.3) proved for the whole boundary set.
    var WALL_DATA = {
      wa: [
        [0.30, -0.080, -0.060], [0.75, -0.080, -0.060], [1.50, -0.080, -0.060],
        [2.20, -0.100, -0.080], [2.60, -0.100, -0.080], [3.00, -0.120, -0.100],
        [3.10, -0.120, -0.100], [3.13, -0.150, -0.120]
      ],
      collar: [
        [0.15, 0.004, 0.200], [0.30, 0.000, 0.080], [0.50, 0.000, 0.030],
        [0.75, 0.000, 0.012], [1.00, 0.000, 0.020]
      ]
    };
    function foldWallStrata(model) {
      if (model !== "noncentered") return [];
      var out = [];
      function withMirror(rows, fn, name0, name1) {
        out.push({name: name0, rows: rows});
        out.push({name: name1, rows: rows.map(function (r) { return [r[0], fn(r[2]), fn(r[1])]; })});
      }
      withMirror(WALL_DATA.wa, function (y) { return -1 - y; },
        "separated fold wall (mixed sector, lower)", "separated fold wall (mixed sector, upper)");
      withMirror(WALL_DATA.collar, function (y) { return 1 - y; },
        "same-sign collar wall (near the mass-zero line)", "same-sign collar wall (near the seam)");
      return out;
    }

    // The lens's own closing points (its two branches meet at a double root
    // of Wtau -- certificates section "the lens closes at its two double
    // roots") and, centered only, the four crossings of the two curves above
    // with the existing straight strata, all exact certified constants
    // rather than located numerically, the same way lineStrata/pointStrata
    // hardcode PI/2 and 0.  Folds onto the SAME point as its mirror twin
    // (beta2* for the lens, the beta = pi endpoint for the potential curve),
    // so only one representative of each pair is listed.
    function extraPointStrata(model) {
      var out = [];
      function add(beta, y, name) {
        var m = massesAt(y);
        out.push({beta: beta, y: y, name: name, rep: [beta, m.s0, m.s1]});
      }
      if (model === "centered") {
        var c = (2 / PI) * Math.atan(2 / PI);
        var beta1Star = 1.420925475551033713494856535;   // CERTIFICATE.md/ARRANGEMENT.md section 2
        add(beta1Star, -0.5, "torque lens closing at its double root");
        add(PI / 2, -c, "torque lens meeting orthogonal teachers");
        add(PI / 2, c - 1, "torque lens meeting orthogonal teachers");
        add(PI / 2, c, "potential curve meeting orthogonal teachers");
        add(PI / 2, 1 - c, "potential curve meeting orthogonal teachers");
        add(0, 0.5, "potential curve closing at coincident teachers");
      } else {
        var betaStar = 2.22566963095871802977134451709;   // CERTIFICATE.md section 3.3
        add(betaStar, -0.5, "torque lens closing at its double root");
      }
      return out;
    }

    function censusAt(K, per, beta, s0, s1) {
      var T = {beta: beta, s0: s0, s1: s1};
      var N = 144, nroot = 0, nmax = 0, npos = 0, nneg = 0;
      var px = 0, pv = Trq(K, T, 0);
      for (var i = 1; i <= N; i += 1) {
        var x = per * i / N, v = Trq(K, T, x);
        if (pv * v < 0) {
          var a = px, b = x, fa = pv;
          for (var k = 0; k < 20; k += 1) {
            var m = 0.5 * (a + b), fm = Trq(K, T, m);
            if (fa * fm <= 0) { b = m; } else { a = m; fa = fm; }
          }
          var r = 0.5 * (a + b);
          nroot += 1;
          if (Tau(K, T, r) > 0) {              // A' < 0: a potential maximum
            nmax += 1;
            if (Pot(K, T, r) > 0) npos += 1; else if (Pot(K, T, r) < 0) nneg += 1;
          }
        }
        px = x; pv = v;
      }
      return nroot + "|" + nmax + "|" + npos + "|" + nneg;
    }

    function censusLabel(key) {
      var p = key.split("|");
      var roots = +p[0], maxima = +p[1], pos = +p[2], neg = +p[3];
      var text = roots + (roots === 1 ? " torque zero" : " torque zeros") + ", " +
        maxima + (maxima === 1 ? " maximum" : " maxima");
      if (maxima > 0 && neg === maxima) text += ", potential negative there";
      else if (maxima > 0 && pos === maxima) text += ", potential positive there";
      return text;
    }

    // Masses on the map's vertical coordinate; never both zero.
    // The mass pair on the projective line, with the SEAM AT s0 = 0.
    // psi = (y+1)*pi/2 runs over (0, pi], one full projective period, so
    //   y = -1, +1  ->  s0 = 0   (the seam, and the least interesting stratum)
    //   y =  0      ->  s1 = 0   (ratio -> +-infinity, ONE line, mid-plot)
    //   y <  0      ->  BOTH teacher masses positive -- the sector the centered
    //                   benignity theorem covers; above it a mass is negative.
    function massesAt(y) {
      var psi = (y + 1) * PI / 2;
      return {s0: Math.sin(psi), s1: Math.cos(psi)};
    }
    // Fold onto the projective line of mass pairs: negating BOTH teacher masses
    // (together with both student masses) is an exact isometry of the loss, so
    // (s0, s1) and (-s0, -s1) are ONE point of the map.  Folding also puts the
    // result in (-1, 1], which atan2 alone does not: atan2 ranges over (-pi, pi],
    // so every teacher with s1 < 0 would otherwise land outside the box.
    function ratioCoord(s0, s1) {
      var a = Math.atan2(s0, s1);          // (-pi, pi]
      a = mod(a, PI);                      // fold by the (-,-) isometry into (0, pi]
      if (a <= 0) a += PI;
      return 2 * a / PI - 1;               // (-1, 1], inverse of massesAt
    }

    // The map is PRECOMPUTED (precompute_census_map.js) because its regions are
    // level sets of the classification the table shows, and that classification
    // costs ~13 ms per teacher -- a grid of it is minutes, not milliseconds.
    // The asset ships beside the page and is fetched once, like the hero's
    // cached frames.
    // The teacher's own census, in exactly the form the generator hashes: the
    // multiset of (family kind, proved type) over the rows the table lists.
    function censusSignature(model, T) {
      var per = model === "centered" ? PI : 2 * PI;
      return buildRows(model, T).rows.map(function (r) {
        var kind = "separate";
        if (r.c && Math.abs(r.c[0]) < 1e-9 && Math.abs(r.c[1]) < 1e-9) kind = "dead";
        else if (/exact|representation|zero student/.test(r.family)) kind = "fit";
        else if (r.theta) {
          var d = mod(r.theta[0] - r.theta[1], per);
          if (d < 1e-4 || per - d < 1e-4) kind = "coincident";
        }
        // The SIGN of the student weights is part of the type: an all-positive
        // trap and an outer-split trap are different phenomena, and only the
        // first is visible to a positive-mass search.  A family whose split is
        // free is a whole line, so its two branches are classified separately.
        var out = [];
        if (r.split && r.typeAt) {
          var mu = r.split.mu;
          var same = Math.abs(mu) > 1e-9 ? r.typeAt(0.6 * mu, 0.4 * mu) : r.typeAt(1, -1);
          var opp = Math.abs(mu) > 1e-9 ? r.typeAt(1.6 * mu, -0.6 * mu) : r.typeAt(1, -1);
          if (/spurious/.test(same)) out.push(kind + ":trap@positive");
          if (/spurious/.test(opp)) out.push(kind + ":trap@mixed");
          if (/global/.test(same) || /global/.test(opp)) out.push(kind + ":global");
        } else {
          var t = r.type || "";
          if (/spurious/.test(t)) {
            var pos = r.c && r.c[0] > -1e-9 && r.c[1] > -1e-9;
            out.push(kind + (pos ? ":trap@positive" : ":trap@mixed"));
          } else if (/global/.test(t)) out.push(kind + ":global");
        }
        return out;
      }).reduce(function (a, b) { return a.concat(b); }, []).sort().join(" | ");
    }

    var censusRequest = null;
    function loadCensusMap(done) {
      var root = typeof window !== "undefined" ? window : globalThis;
      if (root.CensusMap) { done(root.CensusMap); return; }
      if (censusRequest) { censusRequest.push(done); return; }
      censusRequest = [done];
      var el = document.createElement("script");
      el.src = "census-map.js";
      function settle() {
        var waiting = censusRequest || [];
        censusRequest = null;
        waiting.forEach(function (cb) { cb(root.CensusMap || null); });
      }
      el.onload = settle;
      el.onerror = settle;
      document.head.appendChild(el);
    }
    function censusFor(model) {
      var root = typeof window !== "undefined" ? window : globalThis;
      var asset = root.CensusMap;
      return asset && asset.models ? asset.models[model] || null : null;
    }

    // Strata that are single LINES of the map rather than areas: the
    // classification treats each as its own case, and a line has no area to
    // shade, so each is drawn wide enough to see and to recognise as a region.
    // The map's ZERO-MEASURE pieces: strata that are single LINES rather than
    // areas.  The classification treats each as its own case, and a line has no
    // area to shade, so each is drawn wide enough to see and to click.
    //
    // Which mass vanishes on which line follows from the seam convention
    // psi = (y+1)*pi/2, (s0, s1) = (sin psi, cos psi):
    //   y = 0        ->  psi = pi/2  ->  (1, 0)   ->  the SECOND teacher is massless
    //   y = -1, +1   ->  psi = 0, pi ->  (0, ±1)  ->  the FIRST teacher is massless
    // and those two edges are ONE stratum: the y-axis is a circle of period pi,
    // so (0,1) and (0,-1) are the same teacher up to the global sign flip.
    //
    // `rep` is the configuration a click loads, chosen in the interior of the
    // stratum rather than at one of its crossings.
    function lineStrata(model) {
      var per = model === "centered" ? PI : 2 * PI;
      var g = per / 3, r2 = Math.sqrt(0.5);
      var out = [
        {kind: "v", at: 0, name: "coincident teachers", rep: [0, r2, r2]},
        {kind: "h", at: 0, name: "second teacher massless", rep: [g, 1, 0]},
        {kind: "h", at: 1, name: "first teacher massless", rep: [g, 0, -1]},
        {kind: "h", at: -1, name: "first teacher massless", rep: [g, 0, 1]}
      ];
      if (model === "centered") {
        out.push({kind: "v", at: PI / 2, name: "orthogonal teachers", rep: [PI / 2, r2, r2]});
      } else {
        out.push({kind: "v", at: PI, name: "antipodal teachers", rep: [PI, r2, r2]});
      }
      return out;
    }

    // Where two one-dimensional strata cross, the crossing is a piece of its
    // own: one single teacher, of lower dimension than either line through it.
    function pointStrata(model) {
      var per = model === "centered" ? PI : 2 * PI, half = per / 2, out = [];
      var vs = [{at: 0, name: "coincident teachers"},
                {at: half, name: model === "centered" ? "orthogonal teachers" : "antipodal teachers"}];
      var hs = [{at: 0, name: "the second teacher massless", s: [1, 0]},
                {at: -1, name: "the first teacher massless", s: [0, 1]},
                {at: 1, name: "the first teacher massless", s: [0, -1]}];
      vs.forEach(function (v) {
        hs.forEach(function (h) {
          out.push({beta: v.at, y: h.at, name: v.name + " with " + h.name,
                    rep: [v.at, h.s[0], h.s[1]]});
        });
      });
      // The certified curve strata (torque lens, centered's potential curve)
      // add crossings and closing points of their own; see extraPointStrata.
      return out.concat(extraPointStrata(model));
    }

    // The map's own geometry, shared by the renderer and the drag handler.
    var MAP_BOX = {W: 620, H: 330, L: 46, R: 14, T: 12, B: 30};
    var pickedCensus = null;      // a census chosen from the legend, or null
    var viewOf = {};              // live viewBox per panel, for frame placement
    // Zooming shrinks the viewBox, which would otherwise inflate every mark
    // with the geometry.  Strokes are held at screen size by the browser
    // (vector-effect); radii are held by multiplying in the current factor,
    // which is 1 at full view and shrinks as you zoom in.
    var markScale = {map: 1, canvas: 1};
    // The census is invariant under reflecting the gap, beta -> per - beta:
    // measured identical at 1416 of 1416 sampled teachers, in BOTH models.  The
    // second half of the axis therefore carries no information, so the map
    // draws the FUNDAMENTAL DOMAIN [0, per/2] at twice the resolution and shows
    // a teacher from the mirrored half at its reflection.
    function mapFold(model, beta) {
      var per = kernelOf(model).per, b = mod(beta, per);
      return b <= per / 2 ? {beta: b, mirrored: false} : {beta: per - b, mirrored: true};
    }
    var MAP_PAD = 7;                    // data inset, so edge strata are visible
    function mapFrame(model) {
      var per = kernelOf(model).per, half = per / 2;
      var pw = MAP_BOX.W - MAP_BOX.L - MAP_BOX.R, ph = MAP_BOX.H - MAP_BOX.T - MAP_BOX.B;
      var dx = MAP_BOX.L + MAP_PAD, dw = pw - 2 * MAP_PAD;
      var dy = MAP_BOX.T + MAP_PAD, dh = ph - 2 * MAP_PAD;
      return {
        per: per, half: half, pw: pw, ph: ph, dx: dx, dw: dw, dy: dy, dh: dh,
        X: function (beta) { return dx + dw * Math.min(1, mod(beta, per) / half); },
        Y: function (y) { return dy + dh * (1 - y) / 2; },
        betaAt: function (px) { return half * Math.min(1, Math.max(0, (px - dx) / dw)); },
        yAt: function (py) { return Math.min(1, Math.max(-1, 1 - 2 * (py - dy) / dh)); }
      };
    }

    // The gap the MAP shows.  `readTeacher` folds the gap into [0, pi] for both
    // models, because reflecting it is an isometry; the map draws the kernel's
    // whole period, so the marker uses the ORIENTED gap and can sit anywhere on
    // the axis -- otherwise the noncentered marker never leaves the left half.
    function orientedGap(model) {
      if (trapArmed) return TRAP.beta;
      var b0 = parseFloat(document.getElementById("t-beta0").value) * PI;
      var b1 = parseFloat(document.getElementById("t-beta1").value) * PI;
      return mod(b1 - b0, kernelOf(model).per);
    }

    // Colour carries the meaning, not the index: a region with an all-positive
    // trap is red (those are the ones a positive-mass search can meet), a region
    // whose traps need opposite-sign weights is slate, and a region where only
    // the exact fit is a minimum is pale green.  Shade separates one trap from
    // several.
    function typeCounts(key) {
      var n = {cp: 0, cm: 0, sp: 0, sm: 0};
      (key || "").split(" | ").forEach(function (t) {
        if (t === "coincident:trap@positive") n.cp += 1;
        else if (t === "coincident:trap@mixed") n.cm += 1;
        else if (t === "separate:trap@positive") n.sp += 1;
        else if (t === "separate:trap@mixed") n.sm += 1;
      });
      n.pos = n.cp + n.sp; n.mixed = n.cm + n.sm; n.all = n.pos + n.mixed;
      return n;
    }
    // Every census gets its OWN colour: two different censuses never share one.
    // The palette runs warm to cool and the assignment is ordered, so a census
    // carrying an all-positive trap -- the kind a positive-mass search can meet
    // -- still lands at the warm end, but no two of them look alike.
    var censusOrder = {};
    function censusRank(built) {
      var key = built.keys.join("\u0000");
      if (censusOrder[key]) return censusOrder[key];
      var idx = built.keys.map(function (k, i) { return i; });
      idx.sort(function (a, b) {
        var A = typeCounts(built.keys[a]), B = typeCounts(built.keys[b]);
        if ((B.pos > 0) - (A.pos > 0)) return (B.pos > 0) - (A.pos > 0);
        if (A.all !== B.all) return A.all - B.all;
        return built.keys[a] < built.keys[b] ? -1 : 1;
      });
      var rank = {};
      idx.forEach(function (t, k) { rank[t] = k; });
      censusOrder[key] = rank;
      return rank;
    }
    function typeColour(built, id, lifted) {
      if (!built) return "var(--lx-c0)";
      var r = censusRank(built)[id];
      var k = (r === undefined ? 0 : r % 12);
      return "var(--lx-" + (lifted ? "l" : "c") + k + ")";
    }
    // One line per census, so the legend reads as a list rather than a
    // paragraph.  "besides the exact fit" was on every entry -- constant text
    // carries no information, so it is stated once at the foot of the legend.
    function describeType(key) {
      // The paper's compact census label: collided (co) or separated (sep)
      // spurious-minimum families, same-sign (+) or opposite-sign (±) student
      // masses, a repeat carrying a multiplier; exact fits are not counted.
      var n = typeCounts(key), bits = [];
      [[n.cp, "co", "+"], [n.cm, "co", "±"], [n.sp, "sep", "+"], [n.sm, "sep", "±"]].forEach(function (t) {
        if (t[0]) bits.push((t[0] > 1 ? t[0] + "×" : "") + t[1] + "<sup>" + t[2] + "</sup>");
      });
      return bits.length ? bits.join(", ") : "none";
    }

    function renderMapLegend(model, built, selectedId) {
      var box = document.getElementById("ex-map-legend");
      if (!box) return;
      if (!built) { box.innerHTML = "<span>loading the region map…</span>"; return; }
      var order = built.keys.map(function (k, i) { return i; })
        .sort(function (a, b) { return typeCounts(built.keys[a]).all - typeCounts(built.keys[b]).all; });
      box.innerHTML = "<span class=\"lx-row\">" + order.map(function (i) {
        var key = built.keys[i], sel = i === selectedId;
        return "<span class=\"lx-legend-pick" + (sel ? " lx-here" : "") +
          (i === pickedCensus ? " lx-picked" : "") + "\" data-census=\"" + i +
          "\" tabindex=\"0\" role=\"button\" title=\"" + (sel ? "the current teacher's census; " : "") + "click to highlight every region with this census\">" +
          "<svg width=\"15\" height=\"11\" aria-hidden=\"true\"><rect x=\"0.5\" y=\"0.5\" width=\"14\" height=\"10\" fill=\"" +
          typeColour(built, i) + "\" stroke=\"var(--ink)\" stroke-opacity=\"0.35\"/></svg><span class=\"lx-lab\">" +
          describeType(key) + "</span></span>";
      }).join("") + "</span>" +
        "<span class=\"lx-row lx-note\">" +
        "<span><svg width=\"15\" height=\"11\" aria-hidden=\"true\"><circle cx=\"7.5\" cy=\"5.5\" r=\"2.6\" fill=\"var(--ink)\" fill-opacity=\"0.62\"/></svg>teacher network from each census</span>" +
        "<span><svg width=\"15\" height=\"11\" aria-hidden=\"true\"><circle cx=\"7.5\" cy=\"5.5\" r=\"4\" fill=\"none\" stroke=\"var(--accent)\" stroke-width=\"1.6\"/></svg>the current teacher network</span></span>";

      // Clicking a legend entry picks that census: every region carrying it
      // lifts at once, and one example teacher is loaded.  That is the only way
      // to reach a census whose regions are all slivers.
      Array.prototype.forEach.call(box.querySelectorAll(".lx-legend-pick"), function (el) {
        function pick() {
          var id = parseInt(el.getAttribute("data-census"), 10);
          pickedCensus = pickedCensus === id ? null : id;
          if (pickedCensus !== null) {
            var ex = built.pieces.filter(function (q) { return q.id === id; })
              .sort(function (a, b) { return b.cells - a.cells; })[0];
            if (ex) { snapTo(ex.rep, true); return; }
          }
          redrawOnly();
        }
        el.addEventListener("click", pick);
        el.addEventListener("keydown", function (e) {
          if (e.key === "Enter" || e.key === " ") { e.preventDefault(); pick(); }
        });
      });
    }

    function drawMap(model, T) {
      var svg = document.getElementById("ex-map");
      if (!svg) return;
      var F = mapFrame(model);
      var built = censusFor(model);
      var here = ratioCoord(T.s0, T.s1), gap = orientedGap(model);

      var ms = markScale.map;            // 1 at full view, smaller as you zoom in
      clearSvg(svg);
      svg.appendChild(svgEl("rect", {x: MAP_BOX.L, y: MAP_BOX.T, width: F.pw, height: F.ph,
        fill: "var(--paper)", stroke: "var(--rule)", "stroke-width": 1,
        "vector-effect": "non-scaling-stroke"}));

      var selectedId = -1;
      if (built) {
        selectedId = built.keys.indexOf(censusSignature(model, T));
        var nx = built.nx, ny = built.ny;
        // The grid spans a whole period, the plot only the fundamental domain,
        // so clamp at the midline: for rectilinear outlines that is an exact
        // clip, and the part clipped away is the mirror of what is drawn.
        var X = function (gx) { return F.dx + F.dw * Math.min(1, 2 * gx / nx); };
        var Y = function (gy) { return F.dy + F.dh * (1 - gy / ny); };
        // Suppress only a TWIN's dot -- same census, same folded spot.  Two
        // different censuses whose examples happen to land close together must
        // both keep their dot, or a region loses the only way to reach it.
        var drawnDots = [];
        // The example dot of the region holding a network type's default
        // teacher IS that teacher: the page opens with the ring on the dot.
        if (!built._pinned) {
          built._pinned = true;
          var dflt = DEFAULT_TEACHER[model];
          if (dflt) {
            var db = mod((parseFloat(dflt[1]) - parseFloat(dflt[0])) * PI, F.per);
            var dy = ratioCoord(parseFloat(dflt[2]), parseFloat(dflt[3]));
            var px = F.X(mapFold(model, db).beta), py = F.Y(dy);
            var hit = null;
            built.pieces.forEach(function (piece) {
              if (hit || !piece.big || !piece.loops) return;
              piece.loops.forEach(function (loop) {
                if (hit) return;
                var inside = false;
                for (var i = 0, j = loop.length - 1; i < loop.length; j = i++) {
                  var xi = X(loop[i][0]), yi = Y(loop[i][1]), xj = X(loop[j][0]), yj = Y(loop[j][1]);
                  if ((yi > py) !== (yj > py) && px < (xj - xi) * (py - yi) / (yj - yi) + xi) inside = !inside;
                }
                if (inside) hit = piece;
              });
            });
            if (hit) {
              hit.rep = [db, parseFloat(dflt[2]) / 2, parseFloat(dflt[3]) / 2];
              hit.at = [db / F.per, (1 + dy) / 2];
            }
          }
        }
        // The map is drawn folded onto half a period, so a region's two
        // mirror halves (and a piece's reflected twin) land on one another.
        // One dot per VISIBLE region: a dot is skipped when a dot of the same
        // census already sits at the same height on a folded x-range that
        // overlaps this piece's.
        function foldedRange(piece) {
          if (piece._fx) return piece._fx;
          var lo = Infinity, hi = -Infinity, nx = built.nx || 6144;
          (piece.loops || []).forEach(function (loop) {
            loop.forEach(function (q) {
              var u = Math.min(q[0], nx - q[0]) / (nx / 2);
              if (u < lo) lo = u;
              if (u > hi) hi = u;
            });
          });
          piece._fx = [lo, hi];
          return piece._fx;
        }
        function dotIsNew(x, y, census, piece) {
          var r = foldedRange(piece);
          for (var k = 0; k < drawnDots.length; k += 1) {
            var d = drawnDots[k];
            if (d[2] !== census || Math.abs(d[1] - y) >= 7) continue;
            var o = Math.min(r[1], d[3][1]) - Math.max(r[0], d[3][0]);
            if (o > 0.5 * Math.min(r[1] - r[0], d[3][1] - d[3][0])) return false;
          }
          drawnDots.push([x, y, census, r]);
          return true;
        }
        built.pieces.forEach(function (piece) {
          if (!piece.loops || !piece.loops.length) return;
          // entirely in the mirrored half: its twin draws the same ground
          var maxGx = 0;
          piece.loops.forEach(function (loop) {
            loop.forEach(function (q) { if (q[0] > maxGx) maxGx = q[0]; });
          });
          var minGx = Infinity;
          piece.loops.forEach(function (loop) {
            loop.forEach(function (q) { if (q[0] < minGx) minGx = q[0]; });
          });
          if (minGx >= nx / 2) return;
          // The outline path depends on the frame only, not on zoom or the
          // teacher, so it is built once per piece and cached on the asset.
          if (!piece._d) {
            piece._d = piece.loops.map(function (loop) {
              return loop.map(function (q, k) {
                return (k ? "L" : "M") + X(q[0]).toFixed(2) + "," + Y(q[1]).toFixed(2);
              }).join(" ") + " Z";
            }).join(" ");
          }
          var d = piece._d;
          // Highlight by LIFTING the picked census towards white with its own
          // solid token.  Outlining it and dimming the rest produced the
          // boundary artifacts: a piece overlaps its mirror twin, so a
          // translucent fill blends twice and an outline is drawn on both.
          var picked = pickedCensus !== null && piece.id === pickedCensus;
          var node = svgEl("path", {
            d: d, "fill-rule": "evenodd",
            fill: typeColour(built, piece.id, picked), stroke: "none"
          });
          if (piece.big) {
            node.setAttribute("class", "lx-piece");
            node.setAttribute("tabindex", "0");
            node.setAttribute("role", "button");
            node.setAttribute("aria-label", "example configuration in this region");
            node.addEventListener("click", function () {
              if (mapPanMoved) { mapPanMoved = false; return; }
              snapTo(piece.rep);
            });
            node.addEventListener("keydown", function (e) {
              if (e.key === "Enter" || e.key === " ") { e.preventDefault(); snapTo(piece.rep); }
            });
          }
          svg.appendChild(node);
        });
        // one representative per piece: the point a click lands on
        built.pieces.forEach(function (piece) {
          if (!piece.big) return;
          // The dot marks the region's example teacher.  Its stored position is
          // in whole-period coordinates, so fold it into the fundamental domain
          // the same way the outlines are; the teacher it LOADS is unchanged,
          // and a folded one shows the mirror bar on the ring.
          var dx = Math.min(piece.at[0], 1 - piece.at[0]);
          var ex = F.dx + F.dw * 2 * dx, ey = F.dy + F.dh * (1 - piece.at[1]);
          if (!dotIsNew(ex, ey, piece.id, piece)) return;   // the twin already marked it
          svg.appendChild(svgEl("circle", {cx: ex, cy: ey, r: 3 * ms,
            fill: "var(--ink)", "fill-opacity": 0.72, stroke: "var(--paper)",
            "stroke-width": 1, "vector-effect": "non-scaling-stroke",
            "pointer-events": "none"}));
          // The DOT is a click target in its own right.  Some regions are
          // slivers -- the two-mixed-trap census is a thread along the
          // orthogonal-teacher line -- and their area cannot be hit at all.
          hitTarget(svgEl("circle", {cx: ex, cy: ey, r: 11 * ms, fill: "transparent"}),
            piece.rep, "example teacher for this region");
        });
      }

      // THE ENCODING, stated in the legend: a DRAWN LINE is a piece in its own
      // right; two colours meeting with no line between them is a boundary
      // belonging to neither region.  Removing the per-piece outline above is
      // what makes that readable.  The visible stroke stays thin -- width is
      // never a channel -- and a separate transparent stroke over it supplies a
      // finger-sized hit target without being part of the drawing.
      function hitTarget(node, rep, label) {
        node.setAttribute("class", "lx-piece");
        node.setAttribute("tabindex", "0");
        node.setAttribute("role", "button");
        node.setAttribute("aria-label", label);
        node.addEventListener("click", function (e) {
          e.stopPropagation();
          if (mapPanMoved) { mapPanMoved = false; return; }
          snapTo(rep);
        });
        node.addEventListener("keydown", function (e) {
          if (e.key === "Enter" || e.key === " ") { e.preventDefault(); snapTo(rep); }
        });
        svg.appendChild(node);
      }
      // The degenerate teacher networks (coincident, orthogonal/antipodal, a
      // vanishing mass) and the plain fold walls are strata of the critical-point
      // classification, not census boundaries; for a readable map they are not
      // drawn.  DRAW_STRATA keeps the code for the research site.
      var DRAW_STRATA = false;
      if (DRAW_STRATA) lineStrata(model).forEach(function (st) {
        if (st.kind === "v" && (st.at < -1e-9 || st.at > F.half + 1e-9)) return;
        var a = st.kind === "v"
          ? {x1: F.X(st.at), y1: F.dy, x2: F.X(st.at), y2: F.dy + F.dh}
          : {x1: F.dx, y1: F.Y(st.at), x2: F.dx + F.dw, y2: F.Y(st.at)};
        svg.appendChild(svgEl("line", {x1: a.x1, y1: a.y1, x2: a.x2, y2: a.y2,
          stroke: "var(--ink)", "stroke-opacity": 0.6, "stroke-width": 1.6,
          "vector-effect": "non-scaling-stroke", "pointer-events": "none"}));
        // A one-dimensional piece carries an example teacher too, and it needs
        // to be visible and hittable -- the line itself is a thread.
        var sx = st.kind === "v" ? a.x1 : F.X(mapFold(model, st.rep[0]).beta);
        var sy = st.kind === "v" ? F.Y(ratioCoord(st.rep[1], st.rep[2])) : a.y1;
        svg.appendChild(svgEl("circle", {cx: sx, cy: sy, r: 3 * ms,
          fill: "var(--ink)", "fill-opacity": 0.72, stroke: "var(--paper)",
          "stroke-width": 1, "vector-effect": "non-scaling-stroke",
          "pointer-events": "none"}));
        hitTarget(svgEl("line", {x1: a.x1, y1: a.y1, x2: a.x2, y2: a.y2,
          stroke: "transparent", "stroke-width": 18 * ms, "pointer-events": "stroke"}),
          st.rep, st.name + ", a one-dimensional piece of the map");
        hitTarget(svgEl("circle", {cx: sx, cy: sy, r: 11 * ms, fill: "transparent"}),
          st.rep, "example teacher on " + st.name);
      });
      // The certified CURVE strata -- torque lens, centered's potential
      // curve -- drawn exactly like the straight lines above (same colour,
      // same weight: the encoding is "this is a one-dimensional piece", not
      // which one), just traced as a polyline instead of drawn straight.
      if (DRAW_STRATA) curveStrata(model).forEach(function (st) {
        var pts = st.curve;
        if (pts.length < 2) return;
        var d = pts.map(function (q, k) {
          return (k ? "L" : "M") + F.X(q.beta).toFixed(2) + "," + F.Y(q.y).toFixed(2);
        }).join(" ");
        svg.appendChild(svgEl("path", {d: d, fill: "none",
          stroke: "var(--ink)", "stroke-opacity": 0.6, "stroke-width": 1.6,
          "vector-effect": "non-scaling-stroke", "pointer-events": "none"}));
        var mid = pts[Math.floor(pts.length / 2)];
        var m = massesAt(mid.y), rep = [mid.beta, m.s0, m.s1];
        var sx2 = F.X(mid.beta), sy2 = F.Y(mid.y);
        svg.appendChild(svgEl("circle", {cx: sx2, cy: sy2, r: 3 * ms,
          fill: "var(--ink)", "fill-opacity": 0.72, stroke: "var(--paper)",
          "stroke-width": 1, "vector-effect": "non-scaling-stroke",
          "pointer-events": "none"}));
        hitTarget(svgEl("path", {d: d, fill: "none",
          stroke: "transparent", "stroke-width": 18 * ms, "pointer-events": "stroke"}),
          rep, st.name + ", a one-dimensional piece of the map");
        hitTarget(svgEl("circle", {cx: sx2, cy: sy2, r: 11 * ms, fill: "transparent"}),
          rep, "example teacher on " + st.name);
      });
      // The separated stratum's fold walls, drawn through the midpoints of
      // their certified per-beta brackets as one-dimensional strata like the
      // others; the bracket widths are documented in certificates/.
      if (DRAW_STRATA) foldWallStrata(model).forEach(function (st) {
        var rows = st.rows.filter(function (r) { return r[0] <= F.half + 1e-6; });
        if (rows.length < 2) return;
        rows = rows.slice().sort(function (a, b) { return a[0] - b[0]; });
        var mids = rows.map(function (r) { return {beta: r[0], y: (r[1] + r[2]) / 2}; });
        var dLine = mids.map(function (q, k) {
          return (k ? "L" : "M") + F.X(q.beta).toFixed(2) + "," + F.Y(q.y).toFixed(2);
        }).join(" ");
        svg.appendChild(svgEl("path", {d: dLine, fill: "none",
          stroke: "var(--ink)", "stroke-opacity": 0.6, "stroke-width": 1.6,
          "vector-effect": "non-scaling-stroke", "pointer-events": "none"}));
        var mp = mids[Math.floor(mids.length / 2)];
        var m2 = massesAt(mp.y), rep2 = [mp.beta, m2.s0, m2.s1];
        var sx3 = F.X(mp.beta), sy3 = F.Y(mp.y);
        svg.appendChild(svgEl("circle", {cx: sx3, cy: sy3, r: 3 * ms,
          fill: "var(--ink)", "fill-opacity": 0.72, stroke: "var(--paper)",
          "stroke-width": 1, "vector-effect": "non-scaling-stroke",
          "pointer-events": "none"}));
        hitTarget(svgEl("path", {d: dLine, fill: "none",
          stroke: "transparent", "stroke-width": 18 * ms, "pointer-events": "stroke"}),
          rep2, st.name + ", a one-dimensional census");
        hitTarget(svgEl("circle", {cx: sx3, cy: sy3, r: 11 * ms, fill: "transparent"}),
          rep2, "teacher on " + st.name);
      });
      if (DRAW_STRATA) pointStrata(model).forEach(function (pt) {
        var cx = F.X(pt.beta), cy = F.Y(pt.y);
        svg.appendChild(svgEl("circle", {cx: cx, cy: cy, r: 3.6 * ms, fill: "var(--paper)",
          stroke: "var(--ink)", "stroke-width": 1.6, "vector-effect": "non-scaling-stroke",
          "pointer-events": "none"}));
        hitTarget(svgEl("circle", {cx: cx, cy: cy, r: 12 * ms, fill: "transparent"}),
          pt.rep, pt.name + ", a single teacher and a zero-dimensional piece");
      });

      var xt = model === "centered" ? [0, 0.125, 0.25, 0.375, 0.5] : [0, 0.25, 0.5, 0.75, 1];
      xt.forEach(function (u) {
        var x = F.dx + F.dw * u / (F.half / PI);
        svg.appendChild(svgEl("line", {x1: x, y1: MAP_BOX.T + F.ph, x2: x, y2: MAP_BOX.T + F.ph + 4,
          stroke: "var(--rule)", "stroke-width": 1}));
        svg.appendChild(svgEl("text", {x: x, y: MAP_BOX.T + F.ph + 16, "text-anchor": "middle",
          "font-size": 11, fill: "var(--muted)"}, u === 0 ? "0" : u + "π"));
      });
      // the vertical axis is the mass ratio; s1 = 0 is the mid line
      [[-1, "0"], [-0.5, "1"], [0, "±∞"], [0.5, "−1"], [1, "0"]].forEach(function (q) {
        svg.appendChild(svgEl("text", {x: MAP_BOX.L - 6, y: F.Y(q[0]) + 4, "text-anchor": "end",
          "font-size": 11, fill: "var(--muted)"}, q[1]));
      });

      // The current teacher.  A gap in the mirrored half is drawn at its
      // reflection, and the marker says so with a bar rather than by pretending
      // the teacher sits somewhere it does not.
      var fold = mapFold(model, gap), mx = F.X(fold.beta), my = F.Y(here);
      svg.appendChild(svgEl("circle", {cx: mx, cy: my, r: 5 * ms,
        fill: "none", stroke: "var(--accent)", "stroke-width": 2,
        "vector-effect": "non-scaling-stroke", "pointer-events": "none"}));
      if (fold.mirrored) {
        svg.appendChild(svgEl("line", {x1: mx - 8 * ms, y1: my, x2: mx + 8 * ms, y2: my,
          stroke: "var(--accent)", "stroke-width": 1.2,
          "vector-effect": "non-scaling-stroke", "pointer-events": "none"}));
      }
      renderMapLegend(model, built, selectedId);
    }

    // Clicking a piece moves the teacher to that piece's representative.
    function snapTo(rep, keepPick) {
      if (!keepPick) pickedCensus = null;    // any other click clears the pick
      var model = modelName(), per = kernelOf(model).per;
      var b0 = parseFloat(document.getElementById("t-beta0").value);
      var b1 = mod(b0 * PI + rep[0], 2 * PI) / PI;
      var slider = document.getElementById("t-beta1");
      // The slider is discrete, so the representative is written on its own
      // grid rather than left for the browser to round; and the unit mass pair
      // is scaled up so the arrows are visible.  Scaling the teacher masses by a
      // positive factor is an exact symmetry of the census: the torque is linear
      // in them, so the roots and their types are the same.
      b1 = Math.max(0, Math.min(parseFloat(slider.max), b1));
      slider.value = b1.toFixed(6);       // the sliders accept any value, so the ring lands on the dot
      document.getElementById("t-s0").value = (rep[1] * 2).toFixed(5);
      document.getElementById("t-s1").value = (rep[2] * 2).toFixed(5);
      releaseTrap();
      render();
    }


    // Format helpers: angles in units of pi, masses to 3 decimals.
    function fmtA(t, per) { return (mod(t, per) / PI).toFixed(5) + "π"; }
    function fmtM(x) { return (x >= 0 ? "+" : "−") + Math.abs(x).toFixed(5); }
    function fmtPair(a, b) { return fmtM(a) + ", " + fmtM(b); }

    var ids = ["t-beta0", "t-beta1", "t-s0", "t-s1"];
    var selection = null;
    var rows = [];
    var trapArmed = false;

    var REGIMES = {
      centered: {
        model: "centered",
        theorem: "complete centered two-student classification",
        anchor: "thm-centered-two-student-enumeration",
        note: "The canonical projective gap lies in [0, π); β = 0 is the coincident stratum. Angles are read mod π."
      },
      noncentered: {
        model: "noncentered",
        theorem: "complete noncentered two-student classification",
        anchor: "thm-plain-critical-enumeration",
        note: "The oriented gap lies in [0, π]; β = 0 merges the rays and β = π is the antipodal stratum. Angles are read mod 2π."
      }
    };

    function currentRegime() {
      var sel = document.getElementById("ex-regime");
      return REGIMES[sel.value] || REGIMES.centered;
    }
    function modelName() { return currentRegime().model; }
    // The classification is stated in the frame where the first teacher sits at
    // angle 0 and the gap is canonical.  The two sliders carry the teachers'
    // ABSOLUTE angles; the frame is the rotation (and, if the second teacher
    // has gone past the first, the reflection) that maps one to the other.
    // Both are isometries of the loss, so every certified angle, mass and type
    // carries over unchanged, and only the drawing knows about them.
    function readFrame() {
      if (trapArmed) return {phi: 0, sigma: 1};
      var b0 = parseFloat(document.getElementById("t-beta0").value)*PI;
      var b1 = parseFloat(document.getElementById("t-beta1").value)*PI;
      var diff = (b1 - b0) % (2*PI); if (diff < 0) diff += 2*PI;
      return {phi: b0, sigma: diff > PI ? -1 : 1};
    }
    function readTeacher() {
      if (trapArmed) return {beta: TRAP.beta, s0: TRAP.s0, s1: TRAP.s1};
      var b0 = parseFloat(document.getElementById("t-beta0").value)*PI;
      var b1 = parseFloat(document.getElementById("t-beta1").value)*PI;
      var diff = (b1 - b0) % (2*PI); if (diff < 0) diff += 2*PI;
      var gap = diff > PI ? 2*PI - diff : diff;
      return {
        beta: Math.min(gap, 0.99*PI),
        s0: parseFloat(document.getElementById("t-s0").value),
        s1: parseFloat(document.getElementById("t-s1").value)
      };
    }
    // The angle sliders keep the site's grid in both regimes -- 0 to 2 in steps
    // of 0.01, printed as a multiple of pi -- so a distinguished angle is one
    // keypress away and the two regimes do not offer the reader two different
    // controls.  Nothing has to be narrowed for the centered model: a centered
    // direction is an unoriented line, so an angle past pi names the same line
    // as the angle pi below it, and readTeacher folds the gap either way.
    function applyRegimeToControls() {
      var trapRowEl = document.getElementById("ex-trap-row");
      if (trapRowEl) trapRowEl.hidden =
        document.getElementById("ex-regime").value !== "noncentered";
    }

    // The certified positive-weight L_R trap (thm-positive-trap-continuum).
    // Bisection only locates the theorem's bracketed root.
    var TRAP = (function () {
      var sin = Math.sin;
      var E = 1 / 2, p = 3 / 40;
      function Hh(x) { return (PI - x) * sin(x); }
      var Fp = (17 / 40) * (2 * PI - 1 / 2) * sin(3 / 40)
             + (3 / 40) * sin(17 / 40) * sin(1 / 2);
      var Gp = (3 / 40) * (2 * PI - 1 / 2) * sin(17 / 40)
             + (17 / 40) * sin(3 / 40) * sin(1 / 2);
      function Fr(q) {
        return -(1 / 2 * (PI + q - 1 / 2) * sin(q)
               + (PI - q) * sin(1 / 2) * sin(q - 1 / 2));
      }
      function Gr(q) {
        return -(1 / 2 * (PI - q) * sin(q - 1 / 2)
               + (PI - 1 / 2 + q) * sin(1 / 2) * sin(q));
      }
      function locus(q) { return Fp * Gr(q) - Fr(q) * Gp; }
      var a = 3 / 4, b = 755 / 1000, it, m;
      for (it = 0; it < 90; it += 1) {
        m = 0.5 * (a + b);
        if (locus(a) * locus(m) <= 0) b = m; else a = m;
      }
      var q = 0.5 * (a + b), D = -Fr(q);
      var HU = q * sin(q), HV = (q - 1 / 2) * sin(q - 1 / 2);
      var s0 = Fp / D, s1 = 1;
      var c0 = (((PI - p) * sin(p)) * D - Fp * HU) / (Hh(E) * D);
      var c1 = (((PI - 17 / 40) * sin(17 / 40)) * D + Fp * HV) / (Hh(E) * D);
      // Every mass -- teacher and the students the certificate pins -- is
      // scaled by one positive factor, which multiplies the loss by its square
      // and leaves the critical points and their types exactly where they are.
      // The factor puts the largest teacher atom just inside the rim, so the
      // configuration is drawn as large as the fixed mass scale allows.
      var lift = 2.3 / Math.max(Math.abs(s0), Math.abs(s1));
      return {
        q: q,
        beta: PI + 3 / 40 - q,
        s0: s0 * lift,
        s1: s1 * lift,
        theta: [PI + 1 / 2 - q, PI - q],
        c: [c0 * lift, c1 * lift]
      };
    }());

    // One row of the classification table.  angleText/massText may describe a
    // whole family; theta/c is the drawn representative.  When `split` is
    // present the family is a flat line of critical points sharing one total
    // mass, and `typeAt` evaluates that family's exact sign selector at a
    // representative of each sign of the split, so the row states the rule
    // along the whole line rather than a label at one point.
    function branch(o) {
      return {
        family: o.family, selector: o.selector, type: o.type, note: o.note,
        angleText: o.angleText, massText: o.massText,
        theta: o.theta || null, c: o.c || null,
        split: o.split || null, typeAt: o.typeAt || null,
        splitNote: o.splitNote || ""
      };
    }

    // The coincidence selector of the classification, evaluated exactly:
    // positive definite  <=>  tau*P > 0  and  W*a*b*tau < 0.
    // The classification's selector is a pair of SIGN tests, so the tolerance
    // belongs on each FACTOR and never on the product.  An absolute cutoff on
    // `W*c0*c1*tau` is a fourth power: all four factors shrink together near the
    // degenerate teacher strata, so the product crosses a fixed epsilon while its
    // sign never changes, and the map grows a thin boundary beside `s = 0` that
    // is an artifact of the cutoff rather than mathematics.  (Measured: at
    // β = 5.01210 the selector reads −3.99e−10 on one side of y = 0.00225 and
    // −2.15e−9 on the other — same sign, opposite verdicts.)
    // The tolerance must be RELATIVE to the teacher's own scale.  An absolute
    // cutoff was wrong twice over: first on the PRODUCT (a fourth power, fixed
    // earlier), and then still on each FACTOR -- because P, W and tau all carry
    // the teacher mass, so a fixed epsilon is a band of NONZERO width hugging
    // the strata where they vanish.  Certified measurement: that band is
    // 1e-6..1e-5 wide at every beta, and it is what turned the measure-zero
    // lines y = 0 and y = +-1 into four spurious `fit:global` slivers -- four
    // of the eight pieces the mosaic reported.  It is also why the piece count
    // agreed at depths 6, 7 and 8: the band is simply too thin to resolve below
    // about depth 14, so the agreement was never convergence.
    var SEL_REL = 1e-9;
    function coincidenceTypeAt(P, W, tau, a, b, scale) {
      var SEL_EPS = SEL_REL * Math.max(scale || 1, 1e-12);
      if (Math.abs(a) < SEL_EPS || Math.abs(b) < SEL_EPS) return "topological saddle";
      if (Math.abs(tau) < SEL_EPS) return "flat boundary (cell coefficient decides)";
      if (Math.abs(P) < SEL_EPS || Math.abs(W) < SEL_EPS) return "topological saddle";
      if (tau * P > 0 && W * a * b * tau < 0) return "spurious local minimum";
      return "topological saddle";
    }

    // The antipodal endpoint null cubic N(a,b) of the classification.
    function endpointNullCubic(S, a, b) {
      var mu = a + b;
      return a * b * Math.pow(Math.abs(mu), 3)
           - S * (a * Math.pow(Math.abs(b), 3) + b * Math.pow(Math.abs(a), 3));
    }

    // The unified coincidence spurious selector, both kernels, all teacher
    // signs (eq-noncentered-coincidence-selector / eq-centered-mixed-line-
    // selector; for a positive teacher it reduces to "potential max and
    // opposite-sign split").  Returns the required split pattern for a
    // nonglobal minimum, or null if every live split is a saddle.
    function coincidenceMinPattern(P, W, tau) {
      if (P > EPS && tau > EPS) {
        if (W > EPS) return "opposite";
        if (W < -EPS) return "positive";
      }
      if (P < -EPS && tau < -EPS) {
        if (W < -EPS) return "opposite";
        if (W > EPS) return "negative";
      }
      return null;
    }
    function splitFor(pattern, mu) {
      if (pattern === "opposite") {
        return Math.abs(mu) < EPS ? [0.5, -0.5] : [1.4 * mu, -0.4 * mu];
      }
      return [0.6 * mu, 0.4 * mu]; // same-sign split matching mu's sign
    }

    // Torque roots of the plain-ReLU teacher.  Two PAIRS of them collide, at
    // t = 0 and at t = pi, as the teacher gap approaches pi, and a uniform scan
    // steps over both pairs; the root the collision at t = 0 leaves behind
    // sits at t of order (pi - beta)^2, which is inside the margin a scan of
    // the open period never enters.  So the scan runs over the CLOSED period
    // and grades its nodes towards the kink lattice {0, beta, pi, beta+pi},
    // which is where the collisions happen.
    var TRQ_UNIFORM = 1440, TRQ_RING = 60, TRQ_R_MIN = 1e-13, TRQ_R_MAX = 0.05;
    function torqueRootsJ(K, T) {
      var per = K.per, nodes = [], i, j;
      for (i = 0; i <= TRQ_UNIFORM; i += 1) nodes.push(per * i / TRQ_UNIFORM);
      var q = Math.pow(TRQ_R_MAX / TRQ_R_MIN, 1 / (TRQ_RING - 1));
      [0, PI, mod(T.beta, per), mod(T.beta + PI, per)].forEach(function (c) {
        var r = TRQ_R_MIN;
        for (j = 0; j < TRQ_RING; j += 1) {
          if (c + r > 0 && c + r < per) nodes.push(c + r);
          if (c - r > 0 && c - r < per) nodes.push(c - r);
          r *= q;
        }
      });
      nodes.sort(function (a, b) { return a - b; });
      var f = function (t) { return Trq(K, T, t); };
      var sref = Math.abs(T.s0) + Math.abs(T.s1), roots = [];
      var prevX = nodes[0], prevV = f(prevX);
      for (i = 1; i < nodes.length; i += 1) {
        var x = nodes[i];
        if (!(x > prevX)) continue;
        var v = f(x);
        // A sign change between two values that are both at the rounding level
        // of the torque is rounding, not a root -- the graded nodes go down to
        // 1e-13 of a kink, where a degenerate teacher has nothing else to give.
        if (isFinite(prevV) && isFinite(v) && prevV * v < 0 &&
            Math.max(Math.abs(prevV), Math.abs(v)) > 1e-13 * sref) {
          var a = prevX, b = x;
          for (j = 0; j < 80; j += 1) {
            var m = 0.5 * (a + b);
            if (f(a) * f(m) <= 0) b = m; else a = m;
          }
          var r = 0.5 * (a + b);
          if (Math.abs(f(r)) < 1e-7 * sref && !roots.some(function (z) {
            var d = Math.abs(z - r);
            return Math.min(d, per - d) < 1e-5;
          })) roots.push(r);
        }
        prevX = x; prevV = v;
      }
      return roots;
    }
    function torqueRootsOf(model, K, T) {
      if (model === "noncentered") return torqueRootsJ(K, T);
      // At the orthogonal centered teacher beta = pi/2 the two torque atoms
      // have the COMMON kink-lattice zeros t = 0 and t = pi/2, for every
      // choice of teacher masses.  The generic open-period sign-change scan
      // deliberately excludes t = 0 and cannot see a common zero there; the
      // missing row used to turn the exact beta = pi/2 census into a false
      // fit-only column.  Handle the exact stratum algebraically.  This is an
      // equality test on purpose: a nearby nondegenerate teacher must still be
      // classified by its displaced roots, never snapped to this column.
      if (T.beta === PI / 2) return [0, PI / 2];
      return findRoots(function (t) { return Trq(K, T, t); },
                       1e-6, K.per - 1e-6);
    }

    function coincidenceRows(model, K, T, opts) {
      opts = opts || {};
      var out = [];
      var roots = torqueRootsOf(model, K, T);
      var sref = Math.max(Math.abs(T.s0), Math.abs(T.s1), 1);
      var tscale = Math.abs(T.s0) + Math.abs(T.s1);   // the teacher's own scale
      roots.forEach(function (t) {
        var P = Pot(K, T, t), W = Wgt(T, t), tau = Tau(K, T, t);
        var mu = P / K.kap;
        var selName = opts.selector ||
          (model === "centered" ? "CenteredLineExplicitPoint" : "GeneralJSameDirectionSolvedPoint");
        out.push(branch({
          family: "coincidence root",
          selector: selName,
          type: coincidenceTypeAt(P, W, tau, 0.6 * mu, 0.4 * mu, tscale),
          note: "Both students at a root of the torque A(t) = 0, with the total mass pinned by the radial rows and the split entirely free. "
            + "Here P = " + P.toFixed(3) + ", W = " + W.toFixed(3) + ", tau = P - 2W = " + tau.toFixed(3)
            + ". The classification selector is tau*P > 0 and W*c0*c1*tau < 0.",
          angleText: "θ₀ = θ₁ = " + fmtA(t, K.per),
          massText: "c₀ + c₁ = " + mu.toFixed(5),
          theta: [t, t], c: [0.6 * mu, 0.4 * mu],
          split: {mu: mu, sref: sref},
          splitNote: "sweeping the flat line c₀ + c₁ = " + mu.toFixed(3),
          typeAt: function (a, b) { return coincidenceTypeAt(P, W, tau, a, b, tscale); }
        }));
      });
      return out;
    }

    function deadRows(model, K, T) {
      var out = [];
      var roots = torqueRootsOf(model, K, T);
      roots.forEach(function (y) {
        var mu = Pot(K, T, y) / K.kap;
        if (Math.abs(mu) < 1e-8) return;
        var g = function (x) { return mu * K.phi(x - y) - Pot(K, T, x); };
        findRoots(g, 1e-6, K.per - 1e-6).forEach(function (x) {
          var d = mod(x - y, K.per);
          if (d < 1e-4 || K.per - d < 1e-4) return;
          out.push(branch({
            family: "one student dead",
            selector: model === "centered" ? "CenteredMixedFirstStudentOffRoot" : "SignedTeacherJFirstStudentOffExplicit",
            type: "topological saddle",
            note: "The live student sits at a torque root with pinned mass; the dead direction solves the displayed residual-root equation. Reviving the dead student descends.",
            angleText: "dead θ₀ = " + fmtA(x, K.per) + ", live θ₁ = " + fmtA(y, K.per),
            massText: "c₀ = 0, c₁ = " + mu.toFixed(5),
            theta: [x, y], c: [0, mu]
          }));
        });
      });
      return out;
    }

    // With both masses zero the angle derivative vanishes identically, so
    // criticality only asks each parked direction to be a zero of the teacher
    // potential.  That is one family, not one row per pair of zeros.
    function nullRows(model, K, T) {
      var zeros = findRoots(function (t) { return Pot(K, T, t); }, 1e-6, K.per - 1e-6);
      if (!zeros.length) return [];
      return [branch({
        family: "both students dead",
        selector: model === "centered" ? "CenteredMixedNullRoot" : "SignedTeacherJBothStudentsOffExplicit",
        type: "topological saddle",
        note: "With both masses zero the angular derivative vanishes identically, so criticality constrains the directions only through the mass rows: each parked direction must be a zero of the teacher potential, and any zero will do for either student, independently. The row is that whole family; the picture draws one member of it. Reviving either student descends.",
        // No angles: the directions are NOT determined here.  Printing two of
        // them would name a configuration, when what is critical is the family
        // of all choices from the zero set.
        angleText: "each direction at a zero of the teacher potential",
        massText: "c₀ = c₁ = 0",
        theta: [zeros[0], zeros[zeros.length > 1 ? 1 : 0]], c: [0, 0]
      })];
    }

    // Separated critical points of the CENTERED kernel: masses solved from the
    // torque rows, c0 = A(th1)/H(D), c1 = -A(th0)/H(D); the two radial rows are
    // the remaining balance.  A grid scan plus Newton polish locates roots of
    // that balance; each located point is a representative of the proved
    // separated branch (beam family).  The plain-ReLU kernel needs the other
    // elimination and has its own locator below.
    function separatedScanCentered(K, T) {
      function residual(t0, t1) {
        var D = t0 - t1, H = K.H(D);
        if (Math.abs(H) < 5e-3) return null;
        var c0 = Trq(K, T, t1) / H, c1 = -Trq(K, T, t0) / H;
        return [
          K.kap * c0 + K.phi(D) * c1 - Pot(K, T, t0),
          K.phi(D) * c0 + K.kap * c1 - Pot(K, T, t1),
          c0, c1
        ];
      }
      var found = [];
      var N = 110;
      for (var i = 0; i < N; i += 1) {
        for (var j = 0; j < N; j += 1) {
          var t0 = K.per * (i + 0.5) / N, t1 = K.per * (j + 0.5) / N;
          var r = residual(t0, t1);
          if (!r) continue;
          if (Math.abs(r[0]) + Math.abs(r[1]) > 0.35) continue;
          // Newton polish on (t0,t1)
          var x = [t0, t1], ok = false;
          for (var it = 0; it < 60; it += 1) {
            var f = residual(x[0], x[1]);
            if (!f) break;
            var h = 1e-6;
            var fa = residual(x[0] + h, x[1]), fb = residual(x[0], x[1] + h);
            if (!fa || !fb) break;
            var j00 = (fa[0] - f[0]) / h, j01 = (fb[0] - f[0]) / h;
            var j10 = (fa[1] - f[1]) / h, j11 = (fb[1] - f[1]) / h;
            var det = j00 * j11 - j01 * j10;
            if (Math.abs(det) < 1e-12) break;
            var dx0 = (-f[0] * j11 + f[1] * j01) / det;
            var dx1 = (-f[1] * j00 + f[0] * j10) / det;
            x = [x[0] + 0.8 * dx0, x[1] + 0.8 * dx1];
            if (Math.abs(dx0) + Math.abs(dx1) < 1e-12) { ok = true; break; }
          }
          var fr = ok && residual(x[0], x[1]);
          if (!fr || Math.abs(fr[0]) + Math.abs(fr[1]) > 1e-7) continue;
          var p = {t0: mod(x[0], K.per), t1: mod(x[1], K.per), c0: fr[2], c1: fr[3]};
          if (Math.abs(p.c0) < 1e-6 || Math.abs(p.c1) < 1e-6) continue;
          // drop exact fits (listed separately): both students on teacher lattice
          var lat = function (t) {
            var d0 = Math.min(mod(t, K.per), K.per - mod(t, K.per));
            var db = Math.min(mod(t - T.beta, K.per), K.per - mod(t - T.beta, K.per));
            return Math.min(d0, db) < 1e-4;
          };
          if (lat(p.t0) && lat(p.t1)) continue;
          // The two students are UNLABELLED, so orient the pair before storing:
          // comparing the gap against half a period does NOT quotient the swap
          // (at a gap of exactly half a period BOTH orderings pass it, which is
          // how antipodal noncentered pairs came to be listed twice).
          if (p.t1 < p.t0 - 1e-12 || (Math.abs(p.t1 - p.t0) <= 1e-12 && p.c1 < p.c0)) {
            p = {t0: p.t1, t1: p.t0, c0: p.c1, c1: p.c0};
          }
          if (!found.some(function (q) {
            return Math.abs(q.t0 - p.t0) < 1e-4 && Math.abs(q.t1 - p.t1) < 1e-4;
          })) found.push(p);
        }
      }
      return found;
    }

    // -------------------------------------------------------------------
    // SEPARATED CRITICAL POINTS OF THE PLAIN-RELU KERNEL: LOCATE, THEN CHECK
    //
    // Solve the mass pair from the RADIAL rows rather than the torque rows.
    // Their determinant is
    //     massDet = kappa^2 - phi(D)^2 = pi^2 - phi_R(D)^2,     D = t0 - t1,
    // and phi_R has range [0, pi] with phi_R = pi only at D = 0, so massDet
    // vanishes ONLY at student coincidence -- in particular not at the
    // antipodal gap D = pi, where H(D) = 0 and the torque elimination above
    // has no solution at all.  Multiplying the two remaining (torque) rows by
    // that matrix clears the division altogether, leaving the pair
    //     Ea(t0,t1) = h(D) P(t0) + phi(D) A(t0) - kappa A(t1)
    //     Eb(t0,t1) = h(D) P(t1) + kappa A(t0)  - phi(D) A(t1)
    // whose zero set is the separated stratum together with the coincidence
    // stratum D = 0, on which Ea and Eb vanish identically.  Newton runs on
    // (Ea, Eb), which is smooth and division-free, and a candidate is ACCEPTED
    // only against the true balance g below -- a near-degenerate ray on which
    // Ea and Eb are small is not a critical point, and reading criticality off
    // the surrogate alone manufactures whole ladders of them near D = 0.
    var SEP_REL = 1e-9;                     // relative to the teacher's scale
    var SEP_JAC_REL = 1e-7;                 // relative to the Jacobian's size
    // phi_R and H_R together, out of ONE sine and ONE cosine.  Evaluating them
    // apart doubles the trigonometry, and the seed lattice below evaluates both
    // at every one of its nodes -- this pair is the classification's inner loop.
    // The third entry is |sin t|, which the folded angle already has in hand:
    // it costs nothing here and it is what the derivatives below need.
    var ATOM_A = [0, 0, 0], ATOM_B = [0, 0, 0];
    var LOAD_0 = [0, 0, 0], LOAD_1 = [0, 0, 0];
    function atomJ(t, out) {
      var x = mod(t, 2 * PI), sg = 1;
      if (x > PI) { x = 2 * PI - x; sg = -1; }
      var s = Math.sin(x), u = PI - x;
      out[0] = u * Math.cos(x) + s;
      out[1] = sg * u * s;
      out[2] = s;
    }
    // the teacher's potential P(t), torque A(t) and weight W(t) at one angle
    function loadJ(T, t, out) {
      atomJ(t, ATOM_A); atomJ(t - T.beta, ATOM_B);
      out[0] = T.s0 * ATOM_A[0] + T.s1 * ATOM_B[0];
      out[1] = -(T.s0 * ATOM_A[1] + T.s1 * ATOM_B[1]);
      out[2] = T.s0 * ATOM_A[2] + T.s1 * ATOM_B[2];
    }
    function sepMassesJ(K, T, t0, t1) {
      var phiD = K.phi(t0 - t1), det = K.kap * K.kap - phiD * phiD;
      loadJ(T, t0, LOAD_0); loadJ(T, t1, LOAD_1);
      var P0 = LOAD_0[0], P1 = LOAD_1[0];
      return [(K.kap * P0 - phiD * P1) / det, (K.kap * P1 - phiD * P0) / det];
    }
    // The residual pair and its exact Jacobian, out of the same evaluation.
    // Differentiating uses only phi' = -h, h' = phi - 2|sin| (the slope atom),
    // P' = A and A' = 2W - P, and the h(D) A terms cancel in both diagonal
    // entries, which is why the Jacobian is this short:
    //     dEa/dt0 = 2 (phi(D) W(t0) - |sin D| P(t0))
    //     dEa/dt1 = -(phi(D) - 2|sin D|) P(t0) + h(D) A(t0)
    //               - kappa (2 W(t1) - P(t1))
    //     dEb/dt0 = (phi(D) - 2|sin D|) P(t1) + h(D) A(t1)
    //               + kappa (2 W(t0) - P(t0))
    //     dEb/dt1 = 2 (|sin D| P(t1) - phi(D) W(t1))
    function sepSurrogateJ(K, T, t0, t1, out, jac) {
      var D = t0 - t1;
      atomJ(D, ATOM_A);
      var phiD = ATOM_A[0], hD = ATOM_A[1], sD = ATOM_A[2];
      loadJ(T, t0, LOAD_0); loadJ(T, t1, LOAD_1);
      var P0 = LOAD_0[0], A0 = LOAD_0[1], W0 = LOAD_0[2];
      var P1 = LOAD_1[0], A1 = LOAD_1[1], W1 = LOAD_1[2];
      out[0] = hD * P0 + phiD * A0 - K.kap * A1;
      out[1] = hD * P1 + K.kap * A0 - phiD * A1;
      if (jac) {
        var slope = phiD - 2 * sD;
        jac[0] = 2 * (phiD * W0 - sD * P0);
        jac[1] = -slope * P0 + hD * A0 - K.kap * (2 * W1 - P1);
        jac[2] = slope * P1 + hD * A1 + K.kap * (2 * W0 - P0);
        jac[3] = 2 * (sD * P1 - phiD * W1);
      }
      return out;
    }
    function sepBalanceJ(K, T, t0, t1) {
      var c = sepMassesJ(K, T, t0, t1);   // leaves the two loads in LOAD_0/1
      var A0 = LOAD_0[1], A1 = LOAD_1[1], hD = K.H(t0 - t1);
      return [c[1] * hD + A0, c[0] * hD - A1, c[0], c[1]];
    }

    // The seed lattice in the student angle t1, graded towards the small values
    // of the teacher weight W(t) = s0|sin t| + s1|sin(t - beta)|.  Both
    // residuals are 2 pi D W(t1) + O(D^2), so a separated family whose gap D
    // tends to zero forces W(t1) to zero: the near-coincident families -- the
    // ones that carry the census where the teacher is nearly antipodal, at gaps
    // of a few hundredths -- live exactly at the small values of W, and a
    // uniform lattice resolves them only at a spacing no browser can afford.
    // A ring is spent only where W really is small on the teacher's own scale:
    // the four kinks 0, pi, beta, beta+pi carry the values s1|sin beta| and
    // s0|sin beta|, which are small exactly at the degenerate columns beta -> 0
    // and beta -> pi and at the strata where a teacher mass vanishes -- the
    // places the near-coincident families live.  Away from those the kink rings
    // buy nothing and cost most of the scan.
    var SEP_S_UNIFORM = 120, SEP_S_RING = 26, SEP_R_MIN = 1e-4, SEP_R_MAX = 0.15;
    var SEP_W_REL = 0.1;
    function sepLatticeJ(K, T) {
      var per = K.per, pts = [], i, k;
      for (i = 0; i < SEP_S_UNIFORM; i += 1) pts.push(per * i / SEP_S_UNIFORM);
      var W = function (t) { return Wgt(T, t); };
      var wref = SEP_W_REL * (Math.abs(T.s0) + Math.abs(T.s1));
      var centres = [0, PI, mod(T.beta, per), mod(T.beta + PI, per)]
        .filter(function (c) { return Math.abs(W(c)) < wref; });
      var M = 4 * SEP_S_UNIFORM, pv = W(0);
      for (i = 1; i <= M; i += 1) {
        var t = per * i / M, v = W(t);
        if (pv * v < 0) {
          var lo = per * (i - 1) / M, hi = t, flo = pv;
          for (k = 0; k < 60; k += 1) {
            var m = 0.5 * (lo + hi), fm = W(m);
            if (flo * fm <= 0) hi = m; else { lo = m; flo = fm; }
          }
          centres.push(0.5 * (lo + hi));
        }
        pv = v;
      }
      var q = Math.pow(SEP_R_MAX / SEP_R_MIN, 1 / (SEP_S_RING - 1));
      centres.forEach(function (c) {
        var r = SEP_R_MIN;
        for (var j = 0; j < SEP_S_RING; j += 1) {
          pts.push(mod(c + r, per), mod(c - r, per));
          r *= q;
        }
      });
      pts.sort(function (a, b) { return a - b; });
      var out = [pts[0]];
      for (i = 1; i < pts.length; i += 1) {
        if (pts[i] - out[out.length - 1] > 1e-13) out.push(pts[i]);
      }
      out.push(out[0] + per);
      return out;
    }
    // The gap lattice: GEOMETRIC from SEP_D_MIN up to SEP_D_SPLIT, then uniform
    // to the antipodal gap pi.  A uniform gap lattice misses the near-coincident
    // families outright.
    // SEP_D_MIN is a float64 limit, not a modelling choice: the mass solve
    // divides by massDet ~ pi^2 D^2, so the balance below carries the rounding
    // of the residual amplified by 1e-16 / D^2, and a family at a gap of 1e-4
    // cannot be told from a near-solution at all.  Below this gap the locator
    // reports nothing rather than reporting jitter.
    var SEP_D_MIN = 1e-3, SEP_D_SPLIT = 0.35, SEP_D_FINE = 26, SEP_D_COARSE = 38;
    function sepGapsJ() {
      var out = [], q = Math.pow(SEP_D_SPLIT / SEP_D_MIN, 1 / SEP_D_FINE), d = SEP_D_MIN, k;
      for (k = 0; k <= SEP_D_FINE; k += 1) { out.push(d); d *= q; }
      for (k = 1; k <= SEP_D_COARSE; k += 1) {
        out.push(SEP_D_SPLIT + (PI - SEP_D_SPLIT) * k / SEP_D_COARSE);
      }
      return out;
    }

    function separatedScanJ(K, T) {
      var per = K.per, S = sepLatticeJ(K, T), DV = sepGapsJ(), N = S.length;
      var sref = Math.abs(T.s0) + Math.abs(T.s1);
      var P1 = new Float64Array(N), A1 = new Float64Array(N), i, k;
      for (i = 0; i < N; i += 1) {
        loadJ(T, S[i], LOAD_0);
        P1[i] = LOAD_0[0]; A1[i] = LOAD_0[1];
      }
      // A cell of the lattice is kept when the value box of BOTH residuals over
      // its four corners contains zero: sign is scale free, which is what the
      // residuals are not -- they carry the teacher's mass and vanish to first
      // order as the gap closes.
      var cells = [], ga = new Float64Array(N), gb = new Float64Array(N);
      var pa = new Float64Array(N), pb = new Float64Array(N), swap;
      for (k = 0; k < DV.length; k += 1) {
        atomJ(DV[k], ATOM_A);
        var D = DV[k], phiD = ATOM_A[0], hD = ATOM_A[1];
        for (i = 0; i < N; i += 1) {
          loadJ(T, S[i] + D, LOAD_0);
          var p0 = LOAD_0[0], a0 = LOAD_0[1];
          ga[i] = hD * p0 + phiD * a0 - K.kap * A1[i];
          gb[i] = hD * P1[i] + K.kap * a0 - phiD * A1[i];
        }
        if (k > 0) {
          for (i = 0; i + 1 < N; i += 1) {
            if (!(Math.min(pa[i], pa[i + 1], ga[i], ga[i + 1]) <= 0 &&
                  Math.max(pa[i], pa[i + 1], ga[i], ga[i + 1]) >= 0)) continue;
            if (!(Math.min(pb[i], pb[i + 1], gb[i], gb[i + 1]) <= 0 &&
                  Math.max(pb[i], pb[i + 1], gb[i], gb[i + 1]) >= 0)) continue;
            cells.push(0.5 * (S[i] + S[i + 1]), 0.5 * (DV[k - 1] + D));
          }
        }
        swap = pa; pa = ga; ga = swap;
        swap = pb; pb = gb; gb = swap;
      }
      var found = [], escale = sref * K.kap * K.kap;
      for (var c = 0; c < cells.length; c += 2) {
        // Newton is run to a fixed point rather than to a residual level, and
        // the BEST iterate is what leaves the loop: near the folds where two
        // families are about to merge the Jacobian is nearly singular, and a
        // residual cutoff there stops at a whole cloud of points that are one
        // family.  The verdict is the balance test below, not the surrogate.
        var x1 = cells[c], x0 = cells[c] + cells[c + 1];
        var bx0 = x0, bx1 = x1, best = Infinity;
        var f = [0, 0], jj = [0, 0, 0, 0], stall = 0;
        for (var it = 0; it < 30; it += 1) {
          sepSurrogateJ(K, T, x0, x1, f, jj);
          var r = Math.abs(f[0]) + Math.abs(f[1]);
          if (r < best) { best = r; bx0 = x0; bx1 = x1; stall = 0; }
          else if ((stall += 1) > 2) break;    // nothing more to gain
          if (r < 1e-15 * escale) break;
          var det = jj[0] * jj[3] - jj[1] * jj[2];
          if (!(Math.abs(det) > 1e-16 * escale)) break;
          var dx0 = (-f[0] * jj[3] + f[1] * jj[1]) / det;
          var dx1 = (-f[1] * jj[0] + f[0] * jj[2]) / det;
          x0 += dx0; x1 += dx1;
          if (Math.abs(dx0) + Math.abs(dx1) < 1e-14) break;
        }
        x0 = bx0; x1 = bx1;
        // TRANSVERSALITY, against the Jacobian's own size.  A separated family
        // is an isolated point of the balance, and the degenerate strata --
        // beta -> 0, where one teacher atom is left and a student may be parked
        // anywhere with a mass tending to zero -- carry instead a near-curve of
        // near-solutions.  Sampling one of those returns one "family" per row of
        // the gap lattice, and the count then grows with the lattice: measured
        // 41 at this resolution against 104 at a much finer one, at
        // beta = 0.001, all with the same normalised determinant 4e-10.  The
        // separation from a genuine family is three orders wide in either
        // direction (1e-4 to 1e-1 for every certified family measured), which
        // is why a single relative cut does the whole job.  A point this test
        // rejects is one the certificate's Krawczyk step would also decline.
        sepSurrogateJ(K, T, x0, x1, f, jj);
        var jn = Math.max(Math.abs(jj[0]), Math.abs(jj[1]),
                          Math.abs(jj[2]), Math.abs(jj[3]));
        if (!(Math.abs(jj[0] * jj[3] - jj[1] * jj[2]) > SEP_JAC_REL * jn * jn)) continue;
        var g = sepBalanceJ(K, T, x0, x1);
        if (!(Math.abs(g[0]) + Math.abs(g[1]) < SEP_REL * sref)) continue;
        // one live student is the "one student dead" family, listed separately
        if (Math.abs(g[2]) < 1e-7 * sref || Math.abs(g[3]) < 1e-7 * sref) continue;
        var b = x1, gap = x0 - x1;
        if (gap < 0) { b = x0; gap = -gap; }
        gap = mod(gap, per); b = mod(b, per);
        if (gap > PI) { gap = per - gap; b = mod(b + per - gap, per); }
        if (gap < SEP_D_MIN) continue;
        var mtol = 1e-6 + 1e-3 * gap;
        // Two candidates are ONE object either when they agree to the precision
        // the gap allows -- the located position carries the residual divided by
        // the smallest singular value of the Jacobian, and that value falls off
        // with the gap -- or when the balance is satisfied all the way ALONG the
        // segment between them, which is what happens on the degenerate strata:
        // there the solution set is a curve and every seed lands on a different
        // point of it, so that a fixed tolerance reports one family per seed.
        if (found.some(function (z) {
          if (Math.abs(z.s - b) < mtol && Math.abs(z.gap - gap) < mtol) return true;
          if (Math.abs(z.s - b) > 0.05 || Math.abs(z.gap - gap) > 0.05) return false;
          var ms = 0.5 * (z.s + b), mg = 0.5 * (z.gap + gap);
          var mb = sepBalanceJ(K, T, ms + mg, ms);
          return Math.abs(mb[0]) + Math.abs(mb[1]) < SEP_REL * sref;
        })) continue;
        // drop the exact fit (listed separately): both students on the teacher
        // lattice {0, beta}
        var lat = function (t) {
          var d0 = Math.min(mod(t, per), per - mod(t, per));
          var db = Math.min(mod(t - T.beta, per), per - mod(t - T.beta, per));
          return Math.min(d0, db) < 1e-5;
        };
        if (lat(mod(b, per)) && lat(mod(b + gap, per))) continue;
        var cc = sepMassesJ(K, T, b + gap, b);
        var p = {t0: mod(b + gap, per), t1: b, c0: cc[0], c1: cc[1], s: b, gap: gap};
        // The two students are UNLABELLED, so orient the pair before storing.
        if (p.t1 < p.t0 - 1e-12 || (Math.abs(p.t1 - p.t0) <= 1e-12 && p.c1 < p.c0)) {
          p = {t0: p.t1, t1: p.t0, c0: p.c1, c1: p.c0, s: b, gap: gap};
        }
        found.push(p);
      }
      return found;
    }

    // Schur label selector of the noncentered separated stratum
    // (eq-separated-schur): evaluated as the exact closed-form sign test the
    // classification uses.  Returns a type string.
    function schurLabel(K, T, p) {
      var D = p.t0 - p.t1;
      var Hp = K.dH(D), Hd = K.H(D), Fd = K.phi(D);
      var A00 = -p.c0 * p.c1 * Hp + p.c0 * Tau(K, T, p.t0);
      var A11 = -p.c0 * p.c1 * Hp + p.c1 * Tau(K, T, p.t1);
      var A01 = p.c0 * p.c1 * Hp;
      var massDet = K.kap * K.kap - Fd * Fd;
      var T00 = massDet * A00 - K.kap * Hd * Hd * p.c0 * p.c0;
      var T11 = massDet * A11 - K.kap * Hd * Hd * p.c1 * p.c1;
      var T01 = massDet * A01 - Fd * Hd * Hd * p.c0 * p.c1;
      var detT = T00 * T11 - T01 * T01;
      // The pivots are quadratic in the student masses and in the teacher's,
      // so their scale runs over orders of magnitude across the map: the sign
      // test is taken against the form's OWN size, never against a fixed
      // epsilon, which would be a band of nonzero width around every stratum
      // where the form degenerates.
      var Tref = Math.max(Math.abs(T00), Math.abs(T11), Math.abs(T01), 1e-300);
      if (T00 > SEP_REL * Tref && detT > SEP_REL * Tref * Tref) return "spurious local minimum";
      if (T00 < -SEP_REL * Tref || T11 < -SEP_REL * Tref ||
          detT < -SEP_REL * Tref * Tref) return "topological saddle";
      // Neither definite nor indefinite: the second variation is positive
      // SEMI-definite with a null direction, so quadratic order cannot decide
      // the type at all -- the verdict needs the cubic and is an OPEN residual
      // of the classification (Lean: generalJ_separatedPair_degeneratePivot_true_null,
      // "quadratic order is constitutionally unable to decide the boundary").
      // Naming it as though it were a proved type would be a false citation.
      return "undecided at second order";
    }

    function exactFitRow(model, K, T) {
      return branch({
        family: "exact fit",
        selector: model === "centered" ? "SignedCenteredExactFit" : "GeneralJNoncenteredExactFit (first-moment matched)",
        type: "global minimum",
        note: model === "centered"
          ? "The two student atoms reproduce the teacher's signed measure mod π (either assignment order)."
          : "The centered exact fit in the orientation with u = v; the other orientation lifts are not critical for the plain loss.",
        angleText: "0.000π, " + fmtA(T.beta, K.per),
        massText: fmtPair(T.s0, T.s1),
        theta: [0, T.beta], c: [T.s0, T.s1]
      });
    }

    function zeroTeacherRows(model, K, T) {
      return [
        branch({
          family: "zero student measure",
          selector: model === "centered" ? "CenteredZeroTeacherExact" : "ZeroTeacherJExact",
          type: "global minimum",
          note: "The zero teacher field is matched exactly: both students dead (angles free).",
          angleText: "free, free",
          massText: "c₀ = c₁ = 0",
          theta: [0, T.beta], c: [0, 0]
        }),
        branch({
          family: "cancelling coincident pair",
          selector: model === "centered" ? "CenteredZeroTeacherExact" : "ZeroTeacherJExact",
          type: "global minimum",
          note: "The two students share one " + (model === "centered" ? "projective line" : "oriented ray") + " and their masses cancel, so the network represents the zero function; direction and common mass are free.",
          angleText: "θ₀ = θ₁ (free)",
          massText: "c₁ = −c₀ (free)",
          theta: [0, 0], c: [1, -1],
          split: {mu: 0, sref: Math.max(Math.abs(T.s0), Math.abs(T.s1), 1)},
          splitNote: "sweeping the cancelling line c₁ = −c₀",
          typeAt: function () { return "global minimum"; }
        })
      ];
    }

    function oneAtomRows(model, K, T, angle, mass) {
      var out = [
        branch({
          family: "one-atom exact representation",
          selector: model === "centered" ? "CenteredOneTeacherExact" : "OneTeacherJExact",
          type: "global minimum",
          note: "Every live student lies on the teacher " + (model === "centered" ? "line" : "ray") + " and the masses sum to the effective teacher mass, so the loss is zero for every split.",
          angleText: fmtA(angle, K.per) + " (dead slots free)",
          massText: "c₀ + c₁ = " + mass.toFixed(5),
          theta: [angle, angle], c: [mass / 2, mass / 2],
          split: {mu: mass, sref: Math.max(Math.abs(mass), 1)},
          splitNote: "sweeping the exact line c₀ + c₁ = " + mass.toFixed(3),
          typeAt: function () { return "global minimum"; }
        })
      ];
      if (model === "centered") {
        out.push(branch({
          family: "perpendicular pair",
          selector: "CenteredOneTeacherPerpendicular",
          type: "topological saddle",
          note: "Both students on the perpendicular line with total mass 2σ/π. Co-rotating the pair changes the loss at second order by exactly −(2/π)σ² < 0, so every split of this line is a saddle.",
          angleText: "θ₀ = θ₁ = " + fmtA(angle + PI / 2, K.per),
          massText: "c₀ + c₁ = " + (2 * mass / PI).toFixed(5),
          theta: [angle + PI / 2, angle + PI / 2], c: [mass / PI, mass / PI],
          split: {mu: 2 * mass / PI, sref: Math.max(Math.abs(mass), 1)},
          splitNote: "sweeping the perpendicular line c₀ + c₁ = " + (2 * mass / PI).toFixed(3),
          typeAt: function () { return "topological saddle"; }
        }));
      } else {
        out.push(branch({
          family: "antipodal cancellation",
          selector: "OneTeacherJAntipodalCancel",
          type: "topological saddle",
          note: "Both students on the opposite ray, where Φ_R(π) = 0 pins the total mass to zero. The rotation line is exactly flat and the quadratic form degenerates, so the corner cubic supplies the descent.",
          angleText: "θ₀ = θ₁ = " + fmtA(angle + PI, K.per),
          massText: "c₁ = −c₀ (free)",
          theta: [angle + PI, angle + PI], c: [0.5, -0.5],
          split: {mu: 0, sref: Math.max(Math.abs(mass), 1)},
          splitNote: "sweeping the cancelling line c₁ = −c₀",
          typeAt: function () { return "topological saddle"; }
        }));
      }
      return out;
    }

    // ----- centered strata ------------------------------------------------

    function centeredQuarterRows(K, T) {
      var out = [];
      if (T.s0 !== T.s1) return out;
      var u, v;
      if (T.beta === PI / 2) {
        u = PI / 4; v = -PI / 4;
      } else {
        u = T.beta / 2; v = T.beta / 2 + PI / 2;
      }
      var D = u - v, det = K.kap * K.kap - K.phi(D) * K.phi(D);
      var c0 = (K.kap * Pot(K, T, u) - K.phi(D) * Pot(K, T, v)) / det;
      var c1 = (K.kap * Pot(K, T, v) - K.phi(D) * Pot(K, T, u)) / det;
      out.push(branch({
        family: T.beta === PI / 2 ? "orthogonal diagonal" : "quarter pair",
        selector: T.beta === PI / 2 ? "CenteredOrthoQuarterDiagonalPoint" : "CenteredPositiveNonorthQuarterPoint",
        type: "topological saddle",
        note: "Equal teacher masses admit the bisector-and-perpendicular pair (a quarter turn apart); masses are the Cramér solutions of the radial rows. Swapping the two students gives the mirror point.",
        angleText: fmtA(u, K.per) + ", " + fmtA(v, K.per),
        massText: fmtPair(c0, c1),
        theta: [u, v], c: [c0, c1]
      }));
      return out;
    }

    function centeredBeamRows(K, T) {
      if (T.s0 === T.s1) return [];
      return separatedScanCentered(K, T).map(function (p) {
        return branch({
          family: "separated beam root",
          selector: "CenteredPositiveBeamSaddle",
          type: "topological saddle",
          note: "The interlaced separated family: one scalar mass-ratio equation in beam coordinates. The drawn point is the bisection representative of the certified root.",
          angleText: fmtA(p.t0, K.per) + ", " + fmtA(p.t1, K.per),
          massText: fmtPair(p.c0, p.c1),
          theta: [p.t0, p.t1], c: [p.c0, p.c1]
        });
      });
    }

    function centeredRows(T) {
      var K = kernelOf("centered");
      // Centered directions are projective: beta = pi is the same teacher
      // stratum as beta = 0.  The interactive controls already return the
      // canonical half-open representative, while the exported classifier is
      // also used by certificate checks on the displayed closed square.  Fold
      // that duplicate endpoint here so both callers see the same census.
      var beta = mod(T.beta, K.per);
      if (beta !== T.beta) T = {beta: beta, s0: T.s0, s1: T.s1};
      var coincident = T.beta === 0;
      if (coincident) {
        var merged = T.s0 + T.s1;
        return merged === 0
          ? {stratum: "coincident teacher, cancelling total mass", rows: zeroTeacherRows("centered", K, T)}
          : {stratum: "coincident teacher, one effective atom", rows: oneAtomRows("centered", K, T, 0, merged)};
      }
      if (T.s0 === 0 && T.s1 === 0) {
        return {stratum: "two distinct zero-mass teacher slots", rows: zeroTeacherRows("centered", K, T)};
      }
      if (T.s0 === 0 || T.s1 === 0) {
        var angle = T.s0 !== 0 ? 0 : T.beta;
        var mass = T.s0 !== 0 ? T.s0 : T.s1;
        return {stratum: "one effective centered teacher atom", rows: oneAtomRows("centered", K, T, angle, mass)};
      }
      var rows;
      if (T.s0 * T.s1 > 0) {
        rows = [exactFitRow("centered", K, T)]
          .concat(coincidenceRows("centered", K, T))
          .concat(centeredQuarterRows(K, T))
          .concat(centeredBeamRows(K, T));
        return {stratum: "genuine same-sign centered teacher", rows: rows};
      }
      rows = [exactFitRow("centered", K, T)]
        .concat(coincidenceRows("centered", K, T, {selector: "CenteredMixedLineRoot"}))
        .concat(deadRows("centered", K, T))
        .concat(nullRows("centered", K, T));
      return {stratum: "genuine mixed-sign centered teacher", rows: rows};
    }

    // ----- noncentered strata --------------------------------------------

    function cuspHolds(mu, S, a, b) {
      if (a * b > 0) return mu * (mu * mu * mu - S * (a * a + b * b)) > 0;
      if (a * b < 0) return mu !== 0 && mu * S * (Math.abs(a) + Math.abs(b)) - Math.abs(mu) * Math.abs(mu) * Math.abs(mu) > 0;
      return false;
    }

    function antipodalEndpointRows(T, atZero) {
      var mu = atZero ? T.s0 : T.s1;
      var S = T.s0 + T.s1;
      var other = atZero ? T.s1 : T.s0;
      var angle = atZero ? 0 : PI;
      var name = atZero ? "endpoint 0" : "endpoint π";
      var sel = atZero ? "AntipodalSameDirectionZeroJ" : "AntipodalSameDirectionPiJ";
      var sref = Math.max(Math.abs(T.s0), Math.abs(T.s1), 1);
      function typeAt(a, b) {
        if (Math.abs(other) < 1e-9) return "global minimum";
        if (Math.abs(a) < 1e-7 || Math.abs(b) < 1e-7) return "topological saddle";
        return endpointNullCubic(S, a, b) > 1e-9
          ? "spurious local minimum" : "topological saddle";
      }
      return [branch({
        family: "coincident " + name,
        selector: sel,
        type: typeAt(0.6 * mu, 0.4 * mu),
        note: "Both students on a teacher ray, total mass pinned to " + mu.toFixed(3)
          + " with the split free. The second variation is the perfect square (π/2)(c₀x + c₁y)², "
          + "always degenerate, so the null cubic N(c₀,c₁) decides; here S = " + S.toFixed(3) + ".",
        angleText: "θ₀ = θ₁ = " + fmtA(angle, 2 * PI),
        massText: "c₀ + c₁ = " + mu.toFixed(5),
        theta: [angle, angle], c: [0.6 * mu, 0.4 * mu],
        split: {mu: mu, sref: sref},
        splitNote: "sweeping the flat line c₀ + c₁ = " + mu.toFixed(3),
        typeAt: typeAt
      })];
    }

    function antipodalRows(T) {
      var K = kernelOf("noncentered");
      var S = T.s0 + T.s1;
      var out = [
        branch({
          family: "opposite exact fit",
          selector: "AntipodalOppositeExactJ",
          type: "global minimum",
          note: "The oriented student rays and signed masses copy the antipodal teacher.",
          angleText: "0.000π, 1.000π",
          massText: fmtPair(T.s0, T.s1),
          theta: [0, PI], c: [T.s0, T.s1]
        })
      ];
      // interior torque roots t+ = pi s0 / S, t- = pi + pi s1 / S
      if (Math.abs(S) > EPS) {
        [PI * T.s0 / S, PI + PI * T.s1 / S].forEach(function (t, idx) {
          var interior = idx === 0 ? (t > 1e-6 && t < PI - 1e-6) : (t > PI + 1e-6 && t < 2 * PI - 1e-6);
          if (!interior) return;
          var mu = Pot(K, T, t) / PI;
          out.push(branch({
            family: "coincident interior root",
            selector: "AntipodalSameDirectionInteriorSolvedJ",
            type: "topological saddle",
            note: "Interior torque root " + (idx === 0 ? "t₊ = πs₀/S" : "t₋ = π + πs₁/S")
              + ". At an antipodal interior root tau = −P, so tau·P = −P² ≤ 0 and the rigid rotation descends for every split.",
            angleText: "θ₀ = θ₁ = " + fmtA(t, 2 * PI),
            massText: "c₀ + c₁ = " + mu.toFixed(5),
            theta: [t, t], c: [mu / 2, mu / 2],
            split: {mu: mu, sref: Math.max(Math.abs(T.s0), Math.abs(T.s1), 1)},
            splitNote: "sweeping the flat line c₀ + c₁ = " + mu.toFixed(3),
            typeAt: function () { return "topological saddle"; }
          }));
        });
      }
      out = out.concat(antipodalEndpointRows(T, true), antipodalEndpointRows(T, false));
      if (T.s0 === T.s1) {
        out.push(branch({
          family: "opposite perpendicular pair",
          selector: "AntipodalOppositePerpendicularJ",
          type: "topological saddle",
          note: "Equal antipodal masses admit the perpendicular pair with both masses 2s/π.",
          angleText: "0.500π, 1.500π",
          massText: fmtPair(2 * T.s0 / PI, 2 * T.s0 / PI),
          theta: [PI / 2, 3 * PI / 2], c: [2 * T.s0 / PI, 2 * T.s0 / PI]
        }));
      }
      // separated midpoint families (symmetric pairs about a teacher ray)
      [[T.s0, 0, "bisector on the 0 ray"], [T.s1, PI, "bisector on the π ray"]].forEach(function (fam) {
        var muT = fam[0], axis = fam[1];
        var f = function (E) {
          var a = S * Math.sin(E / 2) / (E + Math.sin(E));
          return PI * muT - S * E / 2 - 2 * a * (PI - E) * Math.cos(E / 2);
        };
        findRoots(f, 1e-4, PI - 1e-4, 900).forEach(function (E) {
          var a = S * Math.sin(E / 2) / (E + Math.sin(E));
          out.push(branch({
            family: "separated midpoint pair",
            selector: "AntipodalSeparatedMidpointTransportJ",
            type: "topological saddle",
            note: "Equal-mass pair symmetric about the teacher axis (" + fam[2] + "), opening E from the two displayed scalar equations.",
            angleText: fmtA(axis + E / 2, 2 * PI) + ", " + fmtA(axis - E / 2, 2 * PI),
            massText: fmtPair(a, a),
            theta: [axis + E / 2, axis - E / 2], c: [a, a]
          }));
        });
      });
      return out.concat(deadRows("noncentered", K, T), nullRows("noncentered", K, T));
    }

    function noncenteredGenuineRows(T) {
      var K = kernelOf("noncentered");
      var out = [exactFitRow("noncentered", K, T)];
      if (trapArmed) {
        out.push(branch({
          family: "certified positive-weight trap",
          selector: "L_R_positiveWeight_trap_at_explicit_teacher",
          type: "spurious local minimum",
          note: "Both teacher and student masses are positive. Lean proves strict local minimality and positive loss; bisection only draws the certified root.",
          angleText: fmtA(TRAP.theta[0], 2 * PI) + ", " + fmtA(TRAP.theta[1], 2 * PI),
          massText: fmtPair(TRAP.c[0], TRAP.c[1]),
          theta: TRAP.theta, c: TRAP.c
        }));
      }
      out = out.concat(coincidenceRows("noncentered", K, T));
      if (T.s0 === T.s1) {
        [T.beta / 2, T.beta / 2 + PI].forEach(function (t0) {
          var t1 = t0 + PI;
          out.push(branch({
            family: "opposite bisector pair",
            selector: "GeneralJOppositeNonseparatedSolvedGeometry",
            type: "topological saddle",
            note: "Opposite students exist only at equal teacher masses, on the bisector line; each mass is individually pinned.",
            angleText: fmtA(t0, 2 * PI) + ", " + fmtA(t1, 2 * PI),
            massText: fmtPair(Pot(K, T, t0) / PI, Pot(K, T, t1) / PI),
            theta: [t0, t1], c: [Pot(K, T, t0) / PI, Pot(K, T, t1) / PI]
          }));
        });
      }
      out = out.concat(deadRows("noncentered", K, T), nullRows("noncentered", K, T));
      separatedScanJ(K, T).forEach(function (p) {
        var D = mod(p.t0 - p.t1, 2 * PI);
        // Opposite students are listed above only when the teacher masses are
        // equal, which is the only gap at which that family exists; at every
        // other teacher an antipodal separated pair is a row of its own, and
        // dropping it wholesale is what made the old locator blind at D = pi.
        if (T.s0 === T.s1 && Math.abs(D - PI) < 1e-4) return;
        var label = schurLabel(K, T, p);
        out.push(branch({
          family: "separated pair",
          selector: label === "spurious local minimum"
            ? "GeneralJSolvedSeparatedTransport ∧ positive Schur pivots"
            : label === "topological saddle"
              ? "GeneralJSolvedSeparatedTransport ∧ Schur descent"
              : "GeneralJSolvedSeparatedTransport ∧ rank-one germ",
          type: label,
          note: label === "undecided at second order"
            ? "The second variation is positive semi-definite with a null direction here, so quadratic order cannot decide the type; the verdict needs the cubic term and is an open residual of the classification. This is the separated type-change locus det T = 0."
            : "Masses are the Cramér solutions of the radial rows, whose determinant π² − Φ_R(E)² vanishes only at coincident students; the torque rows are the remaining balance and fix the orientation. The type column evaluates the classification's exact Schur sign selector at this root.",
          angleText: fmtA(p.t0, 2 * PI) + ", " + fmtA(p.t1, 2 * PI),
          massText: fmtPair(p.c0, p.c1),
          theta: [p.t0, p.t1], c: [p.c0, p.c1]
        }));
      });
      return out;
    }

    function noncenteredRows(T) {
      var K = kernelOf("noncentered");
      // The gap lives on the circle: a gap of a whole period is the COINCIDENT
      // teacher, and taking it for a generic one leaves the torque identically
      // zero, so that every angle is a root and the scan below reads rounding.
      if (T.beta !== 0 && mod(T.beta, K.per) === 0) {
        T = {beta: 0, s0: T.s0, s1: T.s1};
      }
      if (T.beta === 0) {
        var merged = T.s0 + T.s1;
        return merged === 0
          ? {stratum: "coincident teacher, cancelling total mass", rows: zeroTeacherRows("noncentered", K, T)}
          : {stratum: "coincident teacher, one effective atom", rows: oneAtomRows("noncentered", K, T, 0, merged)};
      }
      if (T.s0 === 0 && T.s1 === 0) {
        return {stratum: "zero teacher", rows: zeroTeacherRows("noncentered", K, T)};
      }
      if (T.s0 === 0 || T.s1 === 0) {
        var angle = T.s0 !== 0 ? 0 : T.beta;
        var mass = T.s0 !== 0 ? T.s0 : T.s1;
        return {stratum: "one effective oriented teacher atom", rows: oneAtomRows("noncentered", K, T, angle, mass)};
      }
      if (T.beta === PI) {
        return {stratum: "genuine antipodal teacher", rows: antipodalRows(T)};
      }
      return {stratum: "genuine open-gap teacher", rows: noncenteredGenuineRows(T)};
    }

    // A two-student critical point does not know which student is which, so a
    // row and the row with the students exchanged are ONE point.  Every family
    // is expected to emit a canonical representative; this is the guard that no
    // corner case slips past, since a family that folds the swap by an
    // inequality on the gap fails exactly at the self-paired configurations
    // (equal angles, or a gap of half a period).
    function studentPairKey(r, per) {
      if (!r.theta || r.theta.length !== 2) return null;
      var q = function (x) { return Math.round(mod(x, per) / (1e-6)) * 1e-6; };
      var m = function (x) { return Math.round((x || 0) / 1e-6) * 1e-6; };
      var a = [q(r.theta[0]), m(r.c && r.c[0])], b = [q(r.theta[1]), m(r.c && r.c[1])];
      if (b[0] < a[0] || (b[0] === a[0] && b[1] < a[1])) { var t = a; a = b; b = t; }
      return a[0] + "|" + a[1] + "|" + b[0] + "|" + b[1];
    }

    function buildRows(model, T) {
      var classification = model === "centered" ? centeredRows(T) : noncenteredRows(T);
      var per = kernelOf(model).per, seen = {};
      classification.rows = classification.rows.filter(function (r) {
        if (r.split) return true;              // a whole line, not a single point
        var k = studentPairKey(r, per);
        if (k === null) return true;
        var tag = r.family + "#" + k;
        if (seen[tag]) return false;
        seen[tag] = true;
        return true;
      });
      // The separated root finder re-finds the one-student-dead critical points
      // as the c -> 0 limit of the separated balance system (one mass at the
      // level of round-off, the other at the dead row's pinned mass).  Those are
      // the dead rows already listed, not further critical points, and at
      // c = 0 the second-order test is degenerate by construction; drop the
      // duplicates rather than show them "undecided".
      var dead = classification.rows.filter(function (r) { return r.family === "one student dead" && r.theta; });
      function sameDir(a, b) { var d = mod(a - b, per); return d < 2e-3 || per - d < 2e-3; }
      classification.rows = classification.rows.filter(function (r) {
        if (r.family !== "separated pair" || !r.c || !r.theta) return true;
        var small = Math.min(Math.abs(r.c[0]), Math.abs(r.c[1]));
        var large = Math.max(Math.abs(r.c[0]), Math.abs(r.c[1]));
        if (small > 1e-4 * Math.max(large, 1e-9)) return true;
        return !dead.some(function (d) {
          return (sameDir(r.theta[0], d.theta[0]) && sameDir(r.theta[1], d.theta[1])) ||
                 (sameDir(r.theta[0], d.theta[1]) && sameDir(r.theta[1], d.theta[0]));
        });
      });
      classification.rows.forEach(function (r, i) {
        r.key = r.family + "#" + r.selector + "#" + (r.angleText || "") + "#" + i;
        r.row = i;
      });
      return classification;
    }

    function pickSelection(list) {
      if (!list.length) return null;
      var byKey = selection && list.filter(function (r) { return r.key === selection.key; });
      if (byKey && byKey.length) return byKey[0];
      var byFamily = selection && list.filter(function (r) { return r.family === selection.family; });
      if (byFamily && byFamily.length) return byFamily[0];
      var byType = selection && list.filter(function (r) { return r.type === selection.type; });
      if (byType && byType.length) return byType[0];
      // First view: a separated saddle shows the most -- two students off
      // the teacher directions -- so it is preferred before the exact fit.
      if (!selection) {
        var sep = list.filter(function (r) {
          return shortFamily(r) === "separate" && !r.split && r.type && /saddle/.test(r.type);
        });
        if (sep.length) return sep[0];
      }
      var global = list.filter(function (r) { return r.type === "global minimum"; });
      return global.length ? global[0] : list[0];
    }

    // --- drawing ----------------------------------------------------------
    // A unit is a spoke: direction = angle, length = |mass| on a FIXED scale
    // that never moves with the controls, so the picture never rescales.  The
    // colour says teacher or student, as everywhere else on the site; the sign
    // of the mass is the line style, solid for positive and dashed with a
    // hollow head for negative.  In the centered model a direction is an
    // unoriented line, so its spoke is drawn as the whole diameter.
    var RAY = 230, MIN_MARK = 14;
    // The scale covers the whole REACHABLE range instead of clipping into it:
    // the teacher sliders stop at 2.4 each, and the classification's own mass
    // rows can pin a total of twice that, so 4.8 is the largest mass any family
    // can display.  Nothing is ever clamped away, and the scale never moves --
    // the circle of mass one is the reference mark, named in the legend.
    var MASS_FULL = 4.8;                 // the mass that reaches the rim
    var UNIT_RADIUS = RAY/MASS_FULL;     // where mass one sits

    function firstMoment(cs, ths) {
      var x = 0, y = 0;
      cs.forEach(function (c, i) { x += c * Math.cos(ths[i]); y += c * Math.sin(ths[i]); });
      return {x: x, y: y};
    }

    // One atom: an arrow from the origin, its length the size of the mass on the
    // fixed scale.  Colour says teacher or student and nothing else; a negative
    // mass is a dashed shaft with a hollow head.
    function drawUnit(model, angle, mass, kind, frame) {
      angle = frame ? frame.phi + frame.sigma*angle : angle;
      var colour = kind === "teacher" ? css("--negative", "#2166ac")
                                      : css("--positive", "#c82d43");
      var paper = css("--paper", "#fff"), muted = css("--muted", "#666");
      var dead = Math.abs(mass) < 1e-9;
      var cf = markScale.canvas;           // marks keep their screen size under zoom
      var r = Math.min(RAY, RAY*Math.abs(mass)/MASS_FULL);
      var ends = model === "centered" ? [1, -1] : [1];   // an unoriented line is drawn both ways
      ends.forEach(function (sign) {
        var cx = sign*Math.cos(angle), cy = -sign*Math.sin(angle);
        if (dead) {
          // A unit of mass zero has no direction to speak of: its angle is free
          // (the angular force carries a factor c), so it sits at the origin.
          circle(canvas, 0, 0, 5*cf, paper, muted);
          return;
        }
        var head = (kind === "teacher" ? 14 : 10)*cf, tipX = r*cx, tipY = r*cy;
        if (r > head) {
          line(canvas, 0, 0, tipX - head*0.8*cx, tipY - head*0.8*cy, colour, LINE,
               mass < 0 ? "6 4" : null);
        }
        var a = Math.atan2(cy, cx);
        canvas.appendChild(svgEl("polygon", {
          points: [[tipX, tipY],
                   [tipX - head*Math.cos(a - 0.42), tipY - head*Math.sin(a - 0.42)],
                   [tipX - head*Math.cos(a + 0.42), tipY - head*Math.sin(a + 0.42)]]
            .map(function (q) { return q[0].toFixed(2) + "," + q[1].toFixed(2); }).join(" "),
          fill: mass < 0 ? paper : colour, stroke: colour, "stroke-width": LINE,
          "stroke-linejoin": "round"
        }));
      });
    }

    // --- dragging the teachers --------------------------------------------
    // The plot is the teacher's second control surface: dragging an atom sets
    // its mass from the radius and, for the second teacher, the gap from the
    // angle.  The first teacher's direction is the gauge, so only its mass
    // moves.  Everything is written back through the sliders, which stay the
    // single source of truth.
    var grabbed = -1;

    function canvasPoint(event) {
      var rect = canvas.getBoundingClientRect(), box = canvas.viewBox.baseVal;
      var x = box.x + (event.clientX - rect.left)/rect.width*box.width;
      var y = box.y + (event.clientY - rect.top)/rect.height*box.height;
      return {x: x, y: -y};
    }

    function setSlider(id, value) {
      var input = document.getElementById(id);
      input.value = String(value);
      return parseFloat(input.value);
    }

    function nearestTeacher(point) {
      var T = readTeacher(), model = modelName(), frame = readFrame();
      var best = -1, bestGap = Infinity;
      [[0, T.s0], [T.beta, T.s1]].forEach(function (pair, index) {
        var angle = frame.phi + frame.sigma*pair[0];
        var r = Math.min(RAY, RAY*Math.abs(pair[1])/MASS_FULL);
        (model === "centered" ? [1, -1] : [1]).forEach(function (sign) {
          var gap = Math.hypot(point.x - sign*r*Math.cos(angle),
                               point.y - sign*r*Math.sin(angle));
          if (gap < bestGap) { bestGap = gap; best = index; }
        });
      });
      return bestGap < 70 ? best : -1;
    }

    function dragTo(point) {
      if (grabbed < 0) return;
      var massId = grabbed ? "t-s1" : "t-s0", angleId = grabbed ? "t-beta1" : "t-beta0";
      var mass = Math.min(MASS_FULL, Math.hypot(point.x, point.y)*MASS_FULL/RAY);
      var sign = parseFloat(document.getElementById(massId).value) < 0 ? -1 : 1;
      setSlider(massId, (sign*mass).toFixed(2));
      // The angle wraps rather than clamping: dragging past the end of the
      // range comes back at zero, over pi in the centered model and over two pi
      // in the plain one.  The slider reaches 2 in both, so nothing is pinned at
      // the end of the track, and the drag quantizes to the slider's own 0.01
      // grid so the two controls agree on where a teacher can sit.  In the
      // centered model the drag lands on the representative in [0, pi) of the
      // line it points at, which is the same line as its antipode.
      var units = modelName() === "centered" ? 1 : 2;
      var angle = Math.atan2(point.y, point.x)/PI % units;
      if (angle < 0) angle += units;
      angle = Math.round(angle*100)/100;
      if (angle >= units) angle = 0;
      setSlider(angleId, angle.toFixed(2));
      releaseTrap();
      render();
    }

    // The whole line of a split-free family: the shared direction, with one
    // marker where the pinned sum of the two masses sits.
    function drawLine(model, angle, mu, frame) {
      var colour = css("--positive", "#c82d43"), paper = css("--paper", "#fff");
      var a = frame ? frame.phi + frame.sigma*angle : angle, cf = markScale.canvas;
      var ends = model === "centered" ? [1, -1] : [1];
      ends.forEach(function (sign) {
        var cx = sign*Math.cos(a), cy = -sign*Math.sin(a);
        line(canvas, 0, 0, RAY*cx, RAY*cy, colour, LINE);
        if (Math.abs(mu) < 1e-9) return;
        var r = Math.min(RAY, RAY*Math.abs(mu)/MASS_FULL);
        circle(canvas, r*cx, r*cy, 6*cf, mu < 0 ? paper : colour, colour);
      });
      if (Math.abs(mu) < 1e-9) circle(canvas, 0, 0, 6*cf, paper, colour);
    }

    function drawFigure(model, T, chosen) {
      clearSvg(canvas);
      var rule = css("--rule", "#bbb");
      // The frame does NOT zoom: the axes and the outer circle are drawn at the
      // current view's centre and counter-scaled, so they stay put on screen
      // while the mass-one circle and the direction arrows magnify with the
      // zoom.  Zooming is for reading small masses, not for enlarging a border.
      var cv = viewOf["ex-canvas"], cf = markScale.canvas;
      drawAxes(canvas, (RAY + 8) * cf,
        cv ? cv.x + cv.w / 2 : 0, cv ? cv.y + cv.h / 2 : 0);
      canvas.appendChild(svgEl("circle", {cx: 0, cy: 0, r: UNIT_RADIUS, fill: "none",
        stroke: rule, "stroke-width": HAIR, "stroke-dasharray": "3 5",
        "vector-effect": "non-scaling-stroke"}));
      var frame = readFrame();
      [[0, T.s0], [T.beta, T.s1]].forEach(function (pair) {
        drawUnit(model, pair[0], pair[1], "teacher", frame);
      });
      var v = firstMoment([T.s0, T.s1], [0, T.beta]);
      if (!chosen || !chosen.theta || !chosen.c) return;
      if (chosen.split) {
        drawLine(model, chosen.theta[0], chosen.split.mu, frame);
      } else {
        chosen.theta.forEach(function (angle, index) {
          drawUnit(model, angle, chosen.c[index], "student", frame);
        });
      }
      return {u: firstMoment(chosen.c, chosen.theta), v: v};
    }

    // Which kind of configuration a family is, read off the family itself: all
    // masses zero is dead, an exact representation is a fit, one shared
    // direction is coincident, and two directions is separate.
    // The table says the verdict in as few words as it takes: the long names
    // live in the theorems, not in a cell.
    var SHORT_TYPE = {
      "global minimum": "global min",
      "topological saddle": "saddle",
      "spurious local minimum": "spurious min",
      "undecided at second order": "not decided (2nd order flat)"
    };
    function shortType(type) { return SHORT_TYPE[type] || type; }

    // The classification is computed in the canonical frame, with the first
    // teacher at angle zero and the second one turning positively.  The table
    // prints the angles the plot actually draws, so every angle in a row's
    // text is mapped back through that frame.
    function inFrame(text, frame, per) {
      if (!text) return text;
      return text.replace(/(-?[0-9]+\.[0-9]+)π/g, function (all, number) {
        return fmtA(frame.phi + frame.sigma*parseFloat(number)*PI, per);
      });
    }

    function shortFamily(r) {
      var per = currentRegime().model === "centered" ? PI : 2*PI;
      if (r.c && Math.abs(r.c[0]) < 1e-9 && Math.abs(r.c[1]) < 1e-9) return "dead";
      if (/exact|representation|zero student/.test(r.family)) return "fit";
      if (r.theta) {
        var d = mod(r.theta[0] - r.theta[1], per);
        if (d < 1e-4 || per - d < 1e-4) return "coincident";
      }
      return "separate";
    }

    function render() {
      var regime = currentRegime(), model = modelName(), T = readTeacher();
      // The readouts are number fields (typing into one moves its slider) or,
      // in the headless shim, plain elements; angles are shown in units of pi.
      function showVal(id, text) {
        var el = document.getElementById(id);
        if (!el) return;
        if ("value" in el && el.tagName === "INPUT") { if (document.activeElement !== el) el.value = text; }
        else el.textContent = text;
      }
      ["t-beta0", "t-beta1"].forEach(function (id) {
        showVal(id + "-val", parseFloat(document.getElementById(id).value).toFixed(3));
      });
      showVal("t-s0-val", T.s0.toFixed(4));
      showVal("t-s1-val", T.s1.toFixed(4));

      var classification = buildRows(model, T);
      rows = classification.rows;
      var chosen = pickSelection(rows);
      if (chosen) selection = {key: chosen.key, family: chosen.family, type: chosen.type};

      // A family with a free mass split is a whole line of critical points: only
      // the sum of the two masses is pinned, so the picture shows the shared
      // direction as a bare line with one marker at the required sum, and the
      // label is the rule on the line rather than a value at one split.
      // A family whose mass split is free is a whole LINE of critical points, and
      // the type can differ along it.  That is a property of the family, not of
      // the reader's click, so every such row states the rule for both signs of
      // the split -- not only the selected one.  (see shortType below)
      rows.forEach(function (r) {
        if (!r.split || !r.typeAt) return;
        var mu = r.split.mu;
        r.massText = Math.abs(mu) > EPS ? "c₀ + c₁ = " + mu.toFixed(5) : "c₁ = −c₀";
        var same, opposite;
        if (Math.abs(mu) > EPS) {
          same = r.typeAt(0.6*mu, 0.4*mu);
          opposite = r.typeAt(1.6*mu, -0.6*mu);
        } else {
          same = opposite = r.typeAt(1, -1);
        }
        r.type = same;
        r.typeText = same === opposite ? shortType(same)
          : shortType(same) + " (c₀c₁>0), " + shortType(opposite) + " (c₀c₁<0)";
        r.typeStyle = same === opposite ? same
          : (same === "spurious local minimum" || opposite === "spurious local minimum"
             ? "spurious local minimum" : same);
      });

      var frame = readFrame();
      var period = model === "centered" ? PI : 2*PI;
      var body = document.getElementById("ex-crit-body");
      while (body.firstChild) body.removeChild(body.firstChild);
      rows.forEach(function (r) {
        var tr = document.createElement("tr");
        var selected = chosen && r.key === chosen.key;
        if (selected) tr.className = "is-selected";
        tr.setAttribute("aria-selected", selected ? "true" : "false");
        // The branch column says what kind of configuration it is, and nothing
        // else: fit, coincident, separate, or dead.  The Lean constructor that
        // certifies it belongs in the theorem, not in a table cell.
        var famTd = document.createElement("td");
        famTd.textContent = shortFamily(r);
        tr.appendChild(famTd);
        function cell(text, cls) {
          var td = document.createElement("td");
          td.textContent = text;
          if (cls) td.className = cls;
          tr.appendChild(td);
        }
        cell(inFrame(r.angleText, frame, period) || "—", "classification-num");
        cell(r.massText || "—", "classification-num");
        var style = r.typeStyle || r.type;
        cell(r.typeText || shortType(r.type),
             style === "global minimum" ? "label-globalmin"
           : style === "spurious local minimum" ? "label-min"
           : style === "topological saddle" ? "label-saddle"
           : style === "undecided at second order" ? "label-undecided" : "label-notcrit");
        tr.title = r.note;
        tr.tabIndex = 0;
        tr.addEventListener("click", function () {
          selection = {key: r.key, family: r.family, type: r.type};
          render();
        });
        tr.addEventListener("keydown", function (event) {
          if (event.key === "Enter" || event.key === " ") { event.preventDefault(); tr.click(); }
        });
        body.appendChild(tr);
      });

      lastDraw = {model: model, T: T, chosen: chosen};
      var uv = drawFigure(model, T, chosen);
      drawMap(model, T);
      if (!censusFor(model)) {
        loadCensusMap(function () {
          if (modelName() === model) drawMap(model, readTeacher());
        });
      }
      // The one quantity the drawing cannot carry: the first moment of the
      // students against the first moment of the teacher.  The two losses differ
      // by exactly one eighth of the squared mismatch (thm-general-harmonic-
      // split), and where the mismatch vanishes a critical point of one loss is
      // a critical point of the other (thm-general-matched-equivalence).
      // Only in the plain-ReLU regime, where a student is an oriented ray and
      // the configuration therefore fixes its first moment: a centered student
      // is an unoriented LINE, so flipping the representative of any one of them
      // changes u and the reported number would be about the drawing rather than
      // about the family.  Every free-split family here has its two students on
      // one direction with the total mass pinned, so u is the same all along the
      // line and the number is a property of the family, not of the drawn split.
      var matchedText = "";
      if (uv && chosen && chosen.c && model === "noncentered") {
        var dx = uv.u.x - uv.v.x, dy = uv.u.y - uv.v.y;
        var mis = dx * dx + dy * dy;
        matchedText = mis < 1e-9
          ? "First moments matched: the students reproduce the teacher's, so this is a critical point of the centered loss too."
          : "First-moment mismatch ‖u − v‖² = " + mis.toFixed(3) +
            ", so the plain-ReLU loss stands ⅛‖u − v‖² = " + (mis / 8).toFixed(4) +
            " above the centered loss here.";
      }
      // check_classification.js drives this widget headless against a stub DOM
      // that carries only the elements the classification itself needs, so the
      // readout is written when it is there, exactly as every init* here
      // returns early when its own canvas is not.
      var momentNote = document.getElementById("ex-moment-note");
      if (momentNote) momentNote.textContent = matchedText;
      // The left legend: two lines, teacher and student generators.  The marks
      // for a zero mass and for a free split appear only while the picture
      // shows one.
      var sw = function (col, dash) {
        return "<svg width=\"22\" height=\"8\" aria-hidden=\"true\"><line x1=\"1\" y1=\"4\" x2=\"21\" y2=\"4\" stroke=\"" + col + "\" stroke-width=\"1.5\"" + (dash ? " stroke-dasharray=\"6 4\"" : "") + "/></svg>";
      };
      var ring = "<svg width=\"10\" height=\"10\" aria-hidden=\"true\"><circle cx=\"5\" cy=\"5\" r=\"4\" fill=\"none\" stroke=\"var(--muted)\" stroke-width=\"1.5\"/></svg>";
      var teacherDead = Math.abs(T.s0) < 1e-9 || Math.abs(T.s1) < 1e-9;
      var studentDead = !!(chosen && chosen.c && chosen.c.some(function (c) { return Math.abs(c) < 1e-9; }));
      var studentSplit = !!(chosen && chosen.split);
      var extra = (teacherDead || studentDead ? "<span>" + ring + "zero mass generator</span>" : "") +
        (studentSplit ? "<span><svg width=\"34\" height=\"10\" aria-hidden=\"true\"><line x1=\"1\" y1=\"5\" x2=\"33\" y2=\"5\" stroke=\"var(--positive)\" stroke-width=\"1.5\"/><circle cx=\"14\" cy=\"5\" r=\"4\" fill=\"var(--positive)\"/></svg>free split: masses sum to the dot</span>" : "");
      document.getElementById("ex-legend").innerHTML =
        "<span class=\"lx-row\"><b>teacher generators</b>" +
        "<span>" + sw("var(--negative)") + "positive mass</span>" +
        "<span>" + sw("var(--negative)", true) + "negative mass</span>" +
        "<span>" + sw("var(--rule)", true) + "unit circle</span></span>" +
        "<span class=\"lx-row\"><b>student generators</b>" +
        "<span>" + sw("var(--positive)") + "positive mass</span>" +
        "<span>" + sw("var(--positive)", true) + "negative mass</span></span>" +
        (extra ? "<span class=\"lx-row\">" + extra + "</span>" : "");

      var trapNote = document.getElementById("ex-trap-note");
      if (trapNote) trapNote.textContent = trapArmed
        ? "Loaded the example trap: β₂ = π + 3/40 − q with q = " + TRAP.q.toFixed(6)
          + ", the theorem-certified root. Every mass is scaled by one positive factor, which "
          + "multiplies the loss by its square and leaves the critical points and their types "
          + "untouched. Moving a control releases it."
        : "";
    }

    function releaseTrap() { if (trapArmed) { trapArmed = false; selection = null; } }

    ids.forEach(function (id) {
      document.getElementById(id).addEventListener("input", function () {
        releaseTrap();
        render();
      });
    });
    // Each network type opens on its own teacher network: a plain one with
    // eight critical families (two collided and two separated spurious minima,
    // three saddles, the fit), a skip one inside the two-positive-trap lens.  Values are (beta0, beta1 in units of pi, s0, s1).
    var DEFAULT_TEACHER = {
      noncentered: ["0.500", "0.100", "1.0000", "-0.7000"],
      centered: ["0.000", "0.480", "1.6000", "-1.2000"]
    };
    function applyDefaultTeacher(model) {
      var d = DEFAULT_TEACHER[model];
      if (!d) return;
      ["t-beta0", "t-beta1", "t-s0", "t-s1"].forEach(function (id, k) {
        var el = document.getElementById(id);
        if (el) el.value = d[k];
      });
    }
    document.getElementById("ex-regime").addEventListener("change", function () {
      releaseTrap();
      selection = null;
      applyRegimeToControls();
      applyDefaultTeacher(modelName());
      render();
    });
    // Opt-in hook for precompute_census_map.js, which needs the very same
    // classification the table is built from.  A browser never defines this, so
    // it costs the page nothing; it exists so the generated map can never drift
    // from the census it claims to show.
    if (typeof globalThis !== "undefined" && globalThis.__censusExport) {
      globalThis.__censusExport({censusSignature: censusSignature,
                                 kernelOf: kernelOf, massesAt: massesAt,
                                 ratioCoord: ratioCoord,
                                 classificationRows: function (model, T) {
                                   return buildRows(model, T).rows;
                                 }});
    }
    canvas.addEventListener("pointerdown", function (event) {
      var point = canvasPoint(event);
      grabbed = nearestTeacher(point);
      if (grabbed < 0) return;
      canvas.classList.add("dragging");
      if (canvas.setPointerCapture) canvas.setPointerCapture(event.pointerId);
      dragTo(point);
      event.preventDefault();
    });
    canvas.addEventListener("pointermove", function (event) {
      if (grabbed < 0) return;
      dragTo(canvasPoint(event));
    });
    ["pointerup", "pointercancel"].forEach(function (name) {
      canvas.addEventListener(name, function () {
        grabbed = -1;
        canvas.classList.remove("dragging");
      });
    });
    var trapButton = document.getElementById("ex-trap");
    if (trapButton) trapButton.addEventListener("click", function () {
      // Prefer a teacher that is ALREADY a region's example, so the button and
      // the map agree: the loaded teacher lands on a dot the reader can see and
      // click again, instead of at a hand-picked point sitting somewhere in a
      // region with no marker on it.
      var built = censusFor(modelName());
      if (built) {
        var ids = built.keys.map(function (k, i) { return i; })
          .filter(function (i) { return /trap@positive/.test(built.keys[i]); });
        if (ids.length) {
          var best = null;
          built.pieces.forEach(function (q) {
            if (ids.indexOf(q.id) < 0 || !q.big) return;
            if (!best || q.cells > best.cells) best = q;
          });
          if (best) {
            pickedCensus = best.id;
            snapTo(best.rep, true);
            return;
          }
        }
      }
      // Fall back to the certified witness when the asset has not arrived yet.
      trapArmed = true;
      document.getElementById("t-beta0").value = "0";
      document.getElementById("t-beta1").value = String(TRAP.beta / PI);
      document.getElementById("t-s0").value = String(TRAP.s0);
      document.getElementById("t-s1").value = String(TRAP.s1);
      selection = {family: "certified positive-weight trap",
                   type: "spurious local minimum"};
      render();
    });
    // Redraw both panels from the classification `render` already computed.
    // Zooming needs this so marks can keep their screen size, and it must NOT
    // reclassify: a census costs about 13 ms and a wheel emits many events.
    var lastDraw = null;
    function redrawOnly() {
      if (!lastDraw) return;
      drawFigure(lastDraw.model, lastDraw.T, lastDraw.chosen);
      drawMap(lastDraw.model, lastDraw.T);
    }
    var pendingRedraw = false;
    function redrawPanel(svgId) {
      if (!lastDraw) return;
      if (svgId === "ex-map") drawMap(lastDraw.model, lastDraw.T);
      else drawFigure(lastDraw.model, lastDraw.T, lastDraw.chosen);
    }

    // ZOOM AND PAN.  Both panels zoom by shrinking the SVG's own viewBox about
    // the pointer, so the drawing keeps its pixel size and the page never
    // reflows.  The previous version grew the map's width inside a scrolling
    // box; it fought the container and did not work at all.  Nothing here is
    // rebuilt by `render`, which only clears children, so a zoom survives every
    // slider move.
    var mapPanMoved = false;      // read by the piece click handlers below
    function makeView(svgId, boxW, boxH, x0, y0, mayPan, lockCentre) {
      var node = document.getElementById(svgId);
      if (!node) return;
      var v = {x: x0, y: y0, w: boxW, h: boxH};
      viewOf[svgId] = v;
      function clamp() {
        // The window is always FULLY covered by the drawing: the viewBox is
        // kept inside it, so no pan or zoom can expose a gap at an edge.
        v.w = Math.min(v.w, boxW);
        v.h = Math.min(v.h, boxH);
        // A radial plot has one mathematical centre.  Keep coordinate zero at
        // the centre of the visible box so its origin, axes, outer frame and
        // mass-reference circle can never drift apart while zooming.
        if (lockCentre) {
          v.x = x0 + (boxW - v.w) / 2;
          v.y = y0 + (boxH - v.h) / 2;
          return;
        }
        v.x = Math.min(Math.max(v.x, x0), x0 + boxW - v.w);
        v.y = Math.min(Math.max(v.y, y0), y0 + boxH - v.h);
      }
      function apply() {
        clamp();
        node.setAttribute("viewBox",
          v.x.toFixed(3) + " " + v.y.toFixed(3) + " " + v.w.toFixed(3) + " " + v.h.toFixed(3));
        // Marks keep their SCREEN size: record the factor and redraw from the
        // cached classification.  This never re-solves anything -- `redrawOnly`
        // reuses the rows `render` already computed.
        markScale[svgId === "ex-map" ? "map" : "canvas"] = v.w / boxW;
        // The viewBox change is instant; the marker redraw is coalesced to one
        // per animation frame and limited to the zoomed panel, so a burst of
        // wheel events does not rebuild the whole map for each tick.
        if (!pendingRedraw) {
          pendingRedraw = true;
          (window.requestAnimationFrame || function (f) { setTimeout(f, 16); })(function () {
            pendingRedraw = false;
            redrawPanel(svgId);
          });
        }
      }
      function at(event) {                       // pointer, in the svg's own units
        var r = node.getBoundingClientRect();
        if (!r.width || !r.height) return {x: v.x + v.w / 2, y: v.y + v.h / 2};
        return {x: v.x + v.w * (event.clientX - r.left) / r.width,
                y: v.y + v.h * (event.clientY - r.top) / r.height};
      }
      node.addEventListener("wheel", function (event) {
        event.preventDefault();
        var wNew = Math.max(boxW / 12, Math.min(boxW, v.w * Math.exp(event.deltaY * 0.0016)));
        if (wNew === v.w) return;
        if (!lockCentre) {
          var q = at(event), f = wNew / v.w;
          v.x = q.x - (q.x - v.x) * f;
          v.y = q.y - (q.y - v.y) * f;
        }
        v.w = wNew;
        v.h = boxH * (wNew / boxW);
        apply();
      }, {passive: false});
      node.addEventListener("dblclick", function () {
        v.x = x0; v.y = y0; v.w = boxW; v.h = boxH; apply();
      });
      var pan = null;
      node.addEventListener("pointerdown", function (event) {
        if (lockCentre) return;                  // radial origin stays centred
        if (v.w >= boxW - 1e-6) return;          // nothing to pan at full view
        if (mayPan && !mayPan()) return;         // the circle's teacher drag wins
        mapPanMoved = false;                     // a new gesture, not the last pan
        pan = {sx: event.clientX, sy: event.clientY, x: v.x, y: v.y, id: event.pointerId,
               captured: false};
      });
      node.addEventListener("pointermove", function (event) {
        if (!pan) return;
        var r = node.getBoundingClientRect();
        if (!r.width || !r.height) return;
        var dx = event.clientX - pan.sx, dy = event.clientY - pan.sy;
        if (Math.abs(dx) + Math.abs(dy) <= 3) return;
        // a pan that ends over a region must not also count as a click on it
        mapPanMoved = true;
        // Capture only once the pointer has actually moved.  Capturing on
        // pointerdown retargets the eventual click to the svg itself, so no
        // region, line or dot could be clicked while zoomed in.  Without
        // capture the pan would die as soon as the pointer crosses one of the
        // region paths or hit targets, which is most of the surface.
        if (!pan.captured && node.setPointerCapture) {
          pan.captured = true;
          try { node.setPointerCapture(pan.id); } catch (e) {}
        }
        v.x = pan.x - v.w * dx / r.width;
        v.y = pan.y - v.h * dy / r.height;
        apply();
      });
      ["pointerup", "pointercancel"].forEach(function (name) {
        node.addEventListener(name, function (event) {
          if (pan && pan.captured && node.releasePointerCapture) {
            try { node.releasePointerCapture(pan.id); } catch (e) {}
          }
          pan = null;
        });
      });
      node.classList.add("lx-zoomable");
    }
    // The polar plot zooms only about its mathematical origin.  Panning would
    // move that origin and its mass-one circle away from the centre of the
    // panel; teacher atoms remain draggable inside the centred view.
    makeView("ex-canvas", 510, 510, -255, -255, null, true);
    makeView("ex-map", MAP_BOX.W, MAP_BOX.H, 0, 0);

    // Clicking empty map space clears a census picked from the legend, the same
    // way clicking a region does.  Without this the pick could only be undone
    // from the legend, which is not where the reader is looking.
    (function () {
      var mapSvg = document.getElementById("ex-map");
      if (!mapSvg) return;
      mapSvg.addEventListener("click", function (event) {
        if (mapPanMoved) { mapPanMoved = false; return; }
        if (event.target !== mapSvg) return;      // a piece handled it already
        if (pickedCensus === null) return;
        pickedCensus = null;
        redrawOnly();
      });
    }());

    applyRegimeToControls();
    if (typeof globalThis === "undefined" || !globalThis.__censusExport) applyDefaultTeacher(modelName());
    render();
  }

  initLandscapeExplorer();
})();
