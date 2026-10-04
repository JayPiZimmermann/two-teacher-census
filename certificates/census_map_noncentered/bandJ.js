// Shipped-classifier probe over a list of (beta, y) points.
//   node bandJ.js '[[beta,y],...]'   ->  [{beta,y,census}, ...]
const api = require("./harness.js").loadClassification();

const pts = JSON.parse(process.argv[2]);
const out = pts.map(function (p) {
  const m = api.massesAt(p[1]);
  return {beta: p[0], y: p[1],
          census: api.censusSignature("noncentered", {beta: p[0], s0: m.s0, s1: m.s1})};
});
process.stdout.write(JSON.stringify(out));
