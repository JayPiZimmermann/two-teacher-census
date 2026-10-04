// libm attribution harness for validate.py section (a2): prints, for each
// sample, node's OWN sin/cos values of the folded arguments together with the
// seven kernel quantities computed by the literal widgets.js source, so the
// Python side can (i) measure the node-vs-glibc sin/cos difference and
// (ii) recompute the formulas from node's values and check bit-identity.
//
// Row format (23 numbers):
//   t beta
//   [xh xf sin(xh) cos(xh) sin(xf) cos(xf) sin(t)]        for the argument t
//   [xh xf sin(xh) cos(xh) sin(xf) cos(xf) sin(u)]        for u = t - beta
//   phiJ(t) HJ(t) dHJ(t) slopeAtom(t) Wtau Wpot Wwgt
// where xh = mod(., 2 PI) and xf is xh folded by x > PI -> 2 PI - x.
"use strict";
const fs = require("fs");
const path = require("path");

const SRC = fs.readFileSync(
  path.join(__dirname, "..", "..", "website", "widgets.js"), "utf8");
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
  extract("phiJ", ["(PI - x) * Math.cos(x) + Math.sin(x)"]),
  extract("HJ", ["(PI - x) * Math.sin(x)", "(x - PI) * Math.sin(x)"]),
  extract("dHJ", ["(PI - x) * Math.cos(x) - Math.sin(x)"]),
  extract("slopeAtom", ["K.phi(y) - 2 * Math.abs(Math.sin(y))"]),
].join("\n");
const make = new Function(
  "PI",
  body + "\n" +
  "var K = {phi: phiJ, H: HJ, dH: dHJ, kap: PI, per: 2 * PI};\n" +
  "function sA(y) { return slopeAtom(K, y); }\n" +
  "return {mod: mod, phiJ: phiJ, HJ: HJ, dHJ: dHJ, sA: sA};");
const W = make(PI);

const N = parseInt(process.argv[2] || "6000", 10);
const G1 = 0.8191725133961645, G2 = 0.6710436067037893;
let u1 = 0.5, u2 = 0.5;
const out = [];
function pack(arg) {
  const xh = W.mod(arg, 2 * PI);
  const xf = xh > PI ? 2 * PI - xh : xh;
  return [xh, xf, Math.sin(xh), Math.cos(xh), Math.sin(xf), Math.cos(xf),
          Math.sin(arg)];
}
for (let k = 0; k < N; k++) {
  u1 = (u1 + G1) % 1;
  u2 = (u2 + G2) % 1;
  const t = -2 * PI + u1 * 6 * PI;
  const beta = u2 * 2 * PI;
  const uu = t - beta;
  const row = [t, beta].concat(pack(t), pack(uu), [
    W.phiJ(t), W.HJ(t), W.dHJ(t), W.sA(t),
    W.HJ(t) * W.sA(uu) - W.HJ(uu) * W.sA(t),
    W.phiJ(t) * W.HJ(uu) - W.phiJ(uu) * W.HJ(t),
    Math.abs(Math.sin(t)) * W.HJ(uu) - Math.abs(Math.sin(uu)) * W.HJ(t)]);
  out.push(row.map(x => Number(x).toPrecision(17)).join(" "));
}
process.stdout.write(out.join("\n") + "\n");
