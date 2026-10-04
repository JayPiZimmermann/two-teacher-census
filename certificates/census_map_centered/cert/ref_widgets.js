// Reference generator for validate.py section (a): evaluates the LITERAL
// source of website/widgets.js (the centered kernel atoms) on 25000
// quasi-random (t, beta) samples and writes ref_samples.txt.
//
// Extraction is by FUNCTION NAME with brace counting (the file's line numbers
// move; its identifiers do not), with a guard asserting each extracted
// function still contains the expected expression fragment.  The extracted
// text is then eval-ed verbatim, so the numbers below are produced by the
// widget's own arithmetic order.
"use strict";
const fs = require("fs");
const path = require("path");

const SRC = fs.readFileSync(
  path.join(__dirname, "..", "..", "..", "website", "widgets.js"), "utf8");

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

const PI = Math.PI;                       // widgets.js line `var PI = Math.PI;`
if (!SRC.includes("var PI = Math.PI;")) throw new Error("guard failed: PI");

const body = [
  extract("mod", ["x % m"]),
  extract("phiC", ["mod(t, PI)", "(PI / 2 - x) * Math.cos(x) + Math.sin(x)"]),
  extract("HC", ["mod(t, PI)", "(PI / 2 - x) * Math.sin(x)"]),
  extract("dHC", ["mod(t, PI)", "(PI / 2 - x) * Math.cos(x) - Math.sin(x)"]),
  extract("slopeAtom", ["K.phi(y) - 2 * Math.abs(Math.sin(y))"]),
].join("\n");
const make = new Function(
  "PI",
  body + "\n" +
  "var K = {phi: phiC, H: HC, dH: dHC, kap: PI / 2, per: PI};\n" +
  "function sA(y) { return slopeAtom(K, y); }\n" +
  "function Wtau(beta, t) { return HC(t) * sA(t - beta) - HC(t - beta) * sA(t); }\n" +
  "function Wpot(beta, t) { return phiC(t) * HC(t - beta) - phiC(t - beta) * HC(t); }\n" +
  "function Wwgt(beta, t) { return Math.abs(Math.sin(t)) * HC(t - beta) - Math.abs(Math.sin(t - beta)) * HC(t); }\n" +
  "return {phiC: phiC, HC: HC, dHC: dHC, sA: sA, Wtau: Wtau, Wpot: Wpot, Wwgt: Wwgt};");
const W = make(PI);

// deterministic quasi-random samples (additive golden-ratio sequence)
const N = 25000;
const G1 = 0.7548776662466927, G2 = 0.5698402909980532;  // 2-d Kronecker
let rows = [];
let u = 0.5, v = 0.5;
for (let k = 0; k < N; k++) {
  u = (u + G1) % 1;
  v = (v + G2) % 1;
  const t = -2 * PI + u * 5 * PI;         // t in [-2pi, 3pi]
  const beta = v * PI;                    // beta in [0, pi]
  rows.push([t, beta,
             W.phiC(t), W.HC(t), W.dHC(t), W.sA(t),
             W.Wtau(beta, t), W.Wpot(beta, t), W.Wwgt(beta, t)]
            .map(x => x.toPrecision(17)).join(" "));
}
fs.writeFileSync(path.join(__dirname, "ref_samples.txt"), rows.join("\n") + "\n");
console.log("wrote ref_samples.txt (" + N + " samples)");
