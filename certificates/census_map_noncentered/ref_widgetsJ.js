// Reference generator for validate.py section (a): evaluates the LITERAL
// source of website/widgets.js (the noncentered kernel atoms, Pot/Trq, the
// separatedScan residual arithmetic and ratioCoord) on 25000 quasi-random
// (t, beta, t1) samples and writes ref_samplesJ.txt.
//
// Extraction is by FUNCTION NAME with brace counting (the file's line numbers
// move; its identifiers do not), with a guard asserting each extracted
// function still contains the expected expression fragment.  The extracted
// text is eval-ed verbatim, so every reference number is produced by the
// widget's own arithmetic order.
//
// Row format (12 numbers):
//   t beta t1  phiJ(t) HJ(t) dHJ(t) slopeAtom(t) Wtau(beta,t) Wpot(beta,t)
//   Wwgt(beta,t) Dsep(beta,t,t1) y(beta,t)
// with Dsep = H(D)^2 * (F0(e0) F1(e1) - F0(e1) F1(e0)), D = t - t1, the
// division form of the separatedScan balance at the two unit teachers, and
// y = ratioCoord(-HJ(t-beta), HJ(t)).
"use strict";
const fs = require("fs");
const path = require("path");

const SRC = fs.readFileSync(
  path.join(__dirname, "..", "..", "website", "widgets.js"), "utf8");

function extract(name, mustContain) {
  const ix = SRC.indexOf("function " + name + "(");
  if (ix < 0) throw new Error("function " + name + " not found in widgets.js");
  let depth = 0, end = -1;
  for (let i = SRC.indexOf("{", ix); i < SRC.length; i++) {
    if (SRC[i] === "{") depth++;
    else if (SRC[i] === "}") { depth--; if (depth === 0) { end = i + 1; break; } }
  }
  const text = SRC.slice(ix, end);
  for (const frag of mustContain) {
    if (!text.includes(frag)) {
      throw new Error("guard failed: " + name + " no longer contains " + JSON.stringify(frag));
    }
  }
  return text;
}

const PI = Math.PI;
if (!SRC.includes("var PI = Math.PI;")) throw new Error("guard failed: PI");

const body = [
  extract("mod", ["x % m"]),
  extract("phiJ", ["mod(t, 2 * PI)", "if (x > PI) x = 2 * PI - x", "(PI - x) * Math.cos(x) + Math.sin(x)"]),
  extract("HJ", ["mod(t, 2 * PI)", "(PI - x) * Math.sin(x)", "(x - PI) * Math.sin(x)"]),
  extract("dHJ", ["(PI - x) * Math.cos(x) - Math.sin(x)", "Math.sin(x) + (x - PI) * Math.cos(x)"]),
  extract("slopeAtom", ["K.phi(y) - 2 * Math.abs(Math.sin(y))"]),
  extract("Pot", ["T.s0 * K.phi(t) + T.s1 * K.phi(t - T.beta)"]),
  extract("Trq", ["-(T.s0 * K.H(t) + T.s1 * K.H(t - T.beta))"]),
  extract("ratioCoord", ["Math.atan2(s0, s1)", "mod(a, PI)", "2 * a / PI - 1"]),
].join("\n");
const make = new Function(
  "PI",
  body + "\n" +
  "var K = {phi: phiJ, H: HJ, dH: dHJ, kap: PI, per: 2 * PI};\n" +
  "function sA(y) { return slopeAtom(K, y); }\n" +
  "function Wtau(beta, t) { return HJ(t) * sA(t - beta) - HJ(t - beta) * sA(t); }\n" +
  "function Wpot(beta, t) { return phiJ(t) * HJ(t - beta) - phiJ(t - beta) * HJ(t); }\n" +
  "function Wwgt(beta, t) { return Math.abs(Math.sin(t)) * HJ(t - beta) - Math.abs(Math.sin(t - beta)) * HJ(t); }\n" +
  "function balances(beta, s0, s1, t0, t1) {\n" +
  "  var T = {beta: beta, s0: s0, s1: s1};\n" +
  "  var D = t0 - t1, H = HJ(D);\n" +
  "  var c0 = Trq(K, T, t1) / H, c1 = -Trq(K, T, t0) / H;\n" +
  "  return [K.kap * c0 + K.phi(D) * c1 - Pot(K, T, t0),\n" +
  "          K.phi(D) * c0 + K.kap * c1 - Pot(K, T, t1)];\n" +
  "}\n" +
  "function Dsep(beta, t0, t1) {\n" +
  "  var A = balances(beta, 1, 0, t0, t1), B = balances(beta, 0, 1, t0, t1);\n" +
  "  var HD = HJ(t0 - t1);\n" +
  "  return (A[0] * B[1] - B[0] * A[1]) * HD * HD;\n" +
  "}\n" +
  "return {phiJ: phiJ, HJ: HJ, dHJ: dHJ, sA: sA, Wtau: Wtau, Wpot: Wpot,\n" +
  "        Wwgt: Wwgt, Dsep: Dsep, ratioCoord: ratioCoord};");
const W = make(PI);

// deterministic quasi-random samples (additive Kronecker sequence)
const N = 25000;
const G1 = 0.8191725133961645, G2 = 0.6710436067037893, G3 = 0.5497004779019703;
let u1 = 0.5, u2 = 0.5, u3 = 0.5;
const rows = [];
for (let k = 0; k < N; k++) {
  u1 = (u1 + G1) % 1;
  u2 = (u2 + G2) % 1;
  u3 = (u3 + G3) % 1;
  const t = -2 * PI + u1 * 6 * PI;        // t in [-2pi, 4pi]
  const beta = u2 * 2 * PI;               // beta in [0, 2pi]
  const t1 = u3 * 2 * PI;                 // t1 in [0, 2pi]
  rows.push([t, beta, t1,
             W.phiJ(t), W.HJ(t), W.dHJ(t), W.sA(t),
             W.Wtau(beta, t), W.Wpot(beta, t), W.Wwgt(beta, t),
             W.Dsep(beta, t, t1),
             W.ratioCoord(-W.HJ(t - beta), W.HJ(t))]
            .map(x => Number(x).toPrecision(17)).join(" "));
}
fs.writeFileSync(path.join(__dirname, "ref_samplesJ.txt"), rows.join("\n") + "\n");
console.log("wrote ref_samplesJ.txt (" + N + " samples)");
