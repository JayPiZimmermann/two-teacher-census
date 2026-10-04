// libm attribution harness for validate.py section (a2): prints, for each
// sample, node's OWN sin/cos values of the folded arguments together with the
// seven kernel quantities computed by the literal widgets.js source, so that
// the Python side can (i) measure the node-vs-glibc sin/cos difference and
// (ii) recompute the formulas from node's values and check bit-identity.
//
// Row format (17 numbers):
//   t beta  [x_t sin(x_t) cos(x_t) sin(t)]  [x_u sin(x_u) cos(x_u) sin(u)]
//   phiC(t) HC(t) dHC(t) slopeAtom(t) Wtau Wpot Wwgt
// where x_t = mod(t, PI) and u = t - beta.
"use strict";
const fs = require("fs");
const path = require("path");

const SRC = fs.readFileSync(
  path.join(__dirname, "..", "..", "..", "website", "widgets.js"), "utf8");
function extract(name, mustContain) {
  const ix = SRC.indexOf("function " + name + "(");
  if (ix < 0) throw new Error("function " + name + " not found");
  let depth = 0, end = -1;
  for (let i = SRC.indexOf("{", ix); i < SRC.length; i++) {
    if (SRC[i] === "{") depth++;
    else if (SRC[i] === "}") { depth--; if (depth === 0) { end = i + 1; break; } }
  }
  const text = SRC.slice(ix, end);
  for (const frag of mustContain) {
    if (!text.includes(frag)) throw new Error("guard failed: " + name);
  }
  return text;
}
const PI = Math.PI;
const body = [
  extract("mod", ["x % m"]),
  extract("phiC", ["(PI / 2 - x) * Math.cos(x) + Math.sin(x)"]),
  extract("HC", ["(PI / 2 - x) * Math.sin(x)"]),
  extract("dHC", ["(PI / 2 - x) * Math.cos(x) - Math.sin(x)"]),
  extract("slopeAtom", ["K.phi(y) - 2 * Math.abs(Math.sin(y))"]),
].join("\n");
const make = new Function(
  "PI",
  body + "\n" +
  "var K = {phi: phiC, H: HC, dH: dHC, kap: PI / 2, per: PI};\n" +
  "function sA(y) { return slopeAtom(K, y); }\n" +
  "return {mod: mod, phiC: phiC, HC: HC, dHC: dHC, sA: sA};");
const W = make(PI);

const N = parseInt(process.argv[2] || "6000", 10);
const G1 = 0.7548776662466927, G2 = 0.5698402909980532;
let u1 = 0.5, u2 = 0.5;
const out = [];
for (let k = 0; k < N; k++) {
  u1 = (u1 + G1) % 1;
  u2 = (u2 + G2) % 1;
  const t = -2 * PI + u1 * 5 * PI;
  const beta = u2 * PI;
  const uu = t - beta;
  const xt = W.mod(t, PI), xu = W.mod(uu, PI);
  const Wtau = W.HC(t) * W.sA(uu) - W.HC(uu) * W.sA(t);
  const Wpot = W.phiC(t) * W.HC(uu) - W.phiC(uu) * W.HC(t);
  const Wwgt = Math.abs(Math.sin(t)) * W.HC(uu) - Math.abs(Math.sin(uu)) * W.HC(t);
  out.push([t, beta,
            xt, Math.sin(xt), Math.cos(xt), Math.sin(t),
            xu, Math.sin(xu), Math.cos(xu), Math.sin(uu),
            W.phiC(t), W.HC(t), W.dHC(t), W.sA(t), Wtau, Wpot, Wwgt]
           .map(x => x.toPrecision(17)).join(" "));
}
process.stdout.write(out.join("\n") + "\n");
