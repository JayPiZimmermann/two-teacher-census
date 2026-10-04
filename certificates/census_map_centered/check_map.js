// Agreement of the SHIPPED centered classifier and mosaic with the certified
// arrangement, at the current website/ state:
//   (a) censusSignature("centered", .) at the ten certified face
//       representatives, compared with the certified per-face censuses of
//       arrangement_faces.json;
//   (b) the finite-grid flood-fill output of the shipped census-map.js
//       (centered block): group count, key count, and the shipped census at
//       every renderer representative.
// Writes face_widget.json.  Read-only on website/.
"use strict";
const crypto = require("crypto");
const fs = require("fs");
const path = require("path");
const api = require(path.join(__dirname, "..", "census_map_noncentered", "harness.js"))
  .loadClassification();

const ROOT = path.join(__dirname, "..", "..");
const PI = Math.PI;
const KAPPA = PI / 2;
const C = (2 / PI) * Math.atan(2 / PI);
const SPURIOUS = "spurious local minimum";
const FLAT = "flat boundary (cell coefficient decides)";

function sha256(file) {
  return crypto.createHash("sha256").update(fs.readFileSync(file)).digest("hex");
}

function near(a, b) { return Math.abs(a - b) <= 2e-15; }

function rootName(t) {
  if (t === 0) return "0";
  if (t === PI / 2) return "pi/2";
  return String(t / PI) + " pi";
}

function orthogonalCase(spec) {
  const T = {beta: PI / 2, s0: spec.s0, s1: spec.s1};
  const rows = api.classificationRows("centered", T)
    .filter(r => r.family === "coincidence root")
    .sort((a, b) => a.theta[0] - b.theta[0]);
  const signature = api.censusSignature("centered", T);
  const roots = rows.map(function (row) {
    const mu = row.split.mu;
    const same = row.typeAt(0.6 * mu, 0.4 * mu);
    const opposite = row.typeAt(1.6 * mu, -0.6 * mu);
    return {t: rootName(row.theta[0]), total_student_mass: mu,
            same_sign_split: same, opposite_sign_split: opposite};
  });
  const trapped = [];
  roots.forEach(function (r) {
    if (r.same_sign_split === SPURIOUS) trapped.push(r.t + ":same-sign");
    if (r.opposite_sign_split === SPURIOUS) trapped.push(r.t + ":opposite-sign");
  });
  const rootSetCorrect = roots.length === 2 && roots[0].t === "0" && roots[1].t === "pi/2";
  const y = api.ratioCoord(spec.s0, spec.s1);
  const nonStrictCorrect = !spec.non_strict || rows.some(function (row) {
    const n = spec.non_strict;
    return rootName(row.theta[0]) === n.t
      && (n.total_mass === undefined || row.split.mu === n.total_mass)
      && row.typeAt(n.probe_weights[0], n.probe_weights[1]) === n.type;
  });
  const agrees = rootSetCorrect && signature === spec.signature
    && JSON.stringify(trapped) === JSON.stringify(spec.trapped)
    && nonStrictCorrect && (spec.y === undefined || near(y, spec.y));
  return {id: spec.id, teacher_masses: [spec.s0, spec.s1], y: y,
          expected_signature: spec.signature, shipped_signature: signature,
          expected_trapped_roots: spec.trapped,
          expected_non_strict_candidate: spec.non_strict || null,
          shipped_roots: roots, root_set_correct: rootSetCorrect,
          non_strict_candidate_correct: nonStrictCorrect, agrees: agrees};
}

const faces = JSON.parse(fs.readFileSync(
  path.join(__dirname, "arrangement_faces.json"), "utf8"));

const out = {object: "shipped-classifier census at the certified face "
                     + "representatives, and the shipped mosaic's centered "
                     + "finite-grid flood-fill output",
             sources_sha256: {
               "website/widgets.js": sha256(path.join(ROOT, "website", "widgets.js")),
               "website/precompute_census_map.js": sha256(path.join(ROOT, "website", "precompute_census_map.js")),
               "website/census-map.js": sha256(path.join(ROOT, "website", "census-map.js")),
               "certificates/census_map_centered/build_arrangement.py":
                 sha256(path.join(__dirname, "build_arrangement.py")),
               "certificates/census_map_centered/check_map.js": sha256(__filename),
               "certificates/census_map_centered/arrangement_faces.json":
                 sha256(path.join(__dirname, "arrangement_faces.json"))
             },
             faces: {}, orthogonal_column: {}, beta_seam: {}, mosaic: {}};
let disagree = 0;
for (const f of faces.faces) {
  const beta = Number(f.interior_point.beta), y = Number(f.interior_point.y);
  const m = api.massesAt(y);
  const shipped = api.censusSignature("centered", {beta: beta, s0: m.s0, s1: m.s1});
  const same = shipped === f.census;
  if (!same) disagree += 1;
  out.faces[f.id] = {beta: beta, y: y, certified: f.census, shipped: shipped,
                     agree: same};
}
out.faces_disagreeing = disagree;

// The sign-chart ratio (-h(t-beta),h(t)) is 0/0 at beta = pi/2.  The true
// column is instead obtained by evaluating the two common torque roots t=0
// and t=pi/2 directly.  These representatives cover the one/two/one open
// intervals in both teacher-sign sectors and every transition point.  The
// transition masses are exact algebraic multiples of kappa, so no numerical
// root locator carries the verdict.
const MIX1 = "coincident:trap@mixed | fit:global";
const MIX2 = "coincident:trap@mixed | coincident:trap@mixed | fit:global";
const POS1 = "coincident:trap@positive | fit:global";
const POS2 = "coincident:trap@positive | coincident:trap@positive | fit:global";
const columnSpecs = [
  {id: "lower_open_below_c_minus_1", s0: 1, s1: 2,
   signature: MIX1, trapped: ["pi/2:opposite-sign"]},
  {id: "lower_threshold_c_minus_1", s0: 1, s1: KAPPA, y: C - 1,
   signature: MIX1, trapped: ["pi/2:opposite-sign"],
   non_strict: {t: "0", probe_weights: [1, -1], type: FLAT}},
  {id: "lower_open_middle", s0: 1, s1: 1,
   signature: MIX2, trapped: ["0:opposite-sign", "pi/2:opposite-sign"]},
  {id: "lower_threshold_minus_c", s0: KAPPA, s1: 1, y: -C,
   signature: MIX1, trapped: ["0:opposite-sign"],
   non_strict: {t: "pi/2", probe_weights: [1, -1], type: FLAT}},
  {id: "lower_open_above_minus_c", s0: 2, s1: 1,
   signature: MIX1, trapped: ["0:opposite-sign"]},
  {id: "upper_open_below_c", s0: 2, s1: -1,
   signature: POS1, trapped: ["0:same-sign"]},
  {id: "upper_threshold_c", s0: KAPPA, s1: -1, y: C,
   signature: POS1, trapped: ["0:same-sign"],
   non_strict: {t: "pi/2", total_mass: 0, probe_weights: [1, 1],
                type: "topological saddle"}},
  {id: "upper_open_middle", s0: 1, s1: -1,
   signature: POS2, trapped: ["0:same-sign", "pi/2:same-sign"]},
  {id: "upper_threshold_1_minus_c", s0: 1, s1: -KAPPA, y: 1 - C,
   signature: POS1, trapped: ["pi/2:same-sign"],
   non_strict: {t: "0", total_mass: 0, probe_weights: [1, 1],
                type: "topological saddle"}},
  {id: "upper_open_above_1_minus_c", s0: 1, s1: -2,
   signature: POS1, trapped: ["pi/2:same-sign"]}
];
const columnCases = columnSpecs.map(orthogonalCase);
const columnDisagree = columnCases.filter(r => !r.agrees).length;
disagree += columnDisagree;
const faceById = Object.fromEntries(faces.faces.map(f => [f.id, f]));
const caseById = Object.fromEntries(columnCases.map(c => [c.id, c]));

// A visual face pair is deliberately NOT an adjacency edge or a component
// certificate.  The check below verifies only (i) the five exact presentation
// labels emitted by build_arrangement.py, (ii) equality of the two open-face
// signatures, and (iii) agreement with the listed local exact-column witnesses.  It
// neither constructs paths nor excludes paths elsewhere in the parameter map.
const EXPECTED_VISUAL_FACE_PAIRS = [
  {id: "VP1", faces: ["F1", "F4"],
   column_witness_cases: ["lower_open_below_c_minus_1", "lower_open_above_minus_c"],
   column_witness_intervals: ["(-1,c-1)", "(-c,0)"]},
  {id: "VP2", faces: ["F2", "F3"], column_witness_cases: ["lower_open_middle"],
   column_witness_intervals: ["(c-1,-c)"]},
  {id: "VP3", faces: ["F5", "F8"], column_witness_cases: ["upper_open_below_c"],
   column_witness_intervals: ["(0,c)"]},
  {id: "VP4", faces: ["F6", "F9"], column_witness_cases: ["upper_open_middle"],
   column_witness_intervals: ["(c,1-c)"]},
  {id: "VP5", faces: ["F7", "F10"], column_witness_cases: ["upper_open_above_1_minus_c"],
   column_witness_intervals: ["(1-c,1)"]}
];
const emittedPairSpecs = ((faces.visual_face_pairs || {}).pairs || []);
const emittedById = Object.fromEntries(emittedPairSpecs.map(p => [p.id, p]));
const exactPairLabelSet = emittedPairSpecs.length === EXPECTED_VISUAL_FACE_PAIRS.length
  && EXPECTED_VISUAL_FACE_PAIRS.every(function (expected) {
    const emitted = emittedById[expected.id];
    return emitted
      && JSON.stringify(emitted.faces) === JSON.stringify(expected.faces)
      && JSON.stringify(emitted.column_witness_cases)
        === JSON.stringify(expected.column_witness_cases)
      && JSON.stringify(emitted.column_witness_intervals)
        === JSON.stringify(expected.column_witness_intervals);
  });
if (!exactPairLabelSet) disagree += 1;
const visualFacePairs = EXPECTED_VISUAL_FACE_PAIRS.map(function (pair) {
  const columns = pair.column_witness_cases.map(id => caseById[id]);
  const signatures = pair.faces.map(id => faceById[id] && faceById[id].census);
  const faceSignatureAgreement = signatures[0] !== undefined
    && signatures.every(s => s === signatures[0]);
  const localColumnWitnessAgreement = faceSignatureAgreement
    && columns.every(column => Boolean(column && column.agrees
      && signatures[0] === column.shipped_signature));
  const checkPassed = exactPairLabelSet && localColumnWitnessAgreement;
  if (!checkPassed) disagree += 1;
  return {id: pair.id, faces: pair.faces,
          column_witness_cases: pair.column_witness_cases,
          column_witness_intervals: pair.column_witness_intervals,
          face_signatures: signatures,
          column_witness_signatures: columns.map(column =>
            column ? column.shipped_signature : null),
          face_signature_agreement: faceSignatureAgreement,
          local_column_witness_agreement: localColumnWitnessAgreement,
          check_passed: checkPassed};
});
out.orthogonal_column = {
  beta: "pi/2",
  exact_common_torque_roots: ["0", "pi/2"],
  threshold_constant: {symbol: "c", formula: "(2/pi) atan(2/pi)", value: C},
  exact_root_data: {
    "0": {P: "(pi/2) s0 + s1", W: "s1", tau: "(pi/2) s0 - s1"},
    "pi/2": {P: "s0 + (pi/2) s1", W: "s0", tau: "-s0 + (pi/2) s1"}
  },
  interval_census: [
    {y: "(-1,c-1)", traps: 1, split: "opposite-sign", roots: ["pi/2"]},
    {y: "c-1", traps: 1, note: "the t=0 candidate has tau=0"},
    {y: "(c-1,-c)", traps: 2, split: "opposite-sign", roots: ["0", "pi/2"]},
    {y: "-c", traps: 1, note: "the t=pi/2 candidate has tau=0"},
    {y: "(-c,0)", traps: 1, split: "opposite-sign", roots: ["0"]},
    {y: "(0,c)", traps: 1, split: "same-sign", roots: ["0"]},
    {y: "c", traps: 1, note: "the t=pi/2 candidate has zero pinned total mass"},
    {y: "(c,1-c)", traps: 2, split: "same-sign", roots: ["0", "pi/2"]},
    {y: "1-c", traps: 1, note: "the t=0 candidate has zero pinned total mass"},
    {y: "(1-c,1)", traps: 1, split: "same-sign", roots: ["pi/2"]}
  ],
  cases: columnCases,
  cases_disagreeing: columnDisagree,
  visual_face_pairs: {
    count: visualFacePairs.length,
    claim_scope: "explicit non-topological presentation pairs only; no "
      + "connectivity, adjacency, or path-exclusion claim",
    connectivity_checked: false,
    exact_pair_labels_agree: exactPairLabelSet,
    pairs: visualFacePairs
  }
};

// The displayed square duplicates the projective beta seam.  Check that the
// source classifier identifies beta = pi with beta = 0, including the
// cancelling teacher, and that both mass-coordinate seams are trap-free.
const seamMasses = [[1, 2], [1, -1], [0, 1]];
const betaCases = seamMasses.map(function (s) {
  const zero = api.censusSignature("centered", {beta: 0, s0: s[0], s1: s[1]});
  const pi = api.censusSignature("centered", {beta: PI, s0: s[0], s1: s[1]});
  const agrees = zero === pi && !/:trap@/.test(zero);
  if (!agrees) disagree += 1;
  return {teacher_masses: s, beta_0: zero, beta_pi: pi, agrees: agrees};
});
const massCases = [[1, 0], [0, 1], [0, -1]].map(function (s) {
  const signature = api.censusSignature("centered", {beta: PI / 2, s0: s[0], s1: s[1]});
  const agrees = !/:trap@/.test(signature) && /global/.test(signature);
  if (!agrees) disagree += 1;
  return {teacher_masses: s, signature: signature, trap_free: agrees};
});
out.beta_seam = {identification: "beta = 0 ~ beta = pi", cases: betaCases,
                 mass_coordinate_boundaries_at_beta_pi_over_2: massCases};

// the shipped asset
const src = fs.readFileSync(path.join(__dirname, "..", "..", "website", "census-map.js"), "utf8");
const w = {};
new Function("window", src)(w);
const M = w.CensusMap.models.centered;
out.mosaic.depth = w.CensusMap.depth;
out.mosaic.keys = M.keys.length;
out.mosaic.algorithm_scope = "finite-grid 4-neighbour flood fill with the "
  + "source generator's presentation separator; not a topology certificate";
out.mosaic.flood_fill_groups = M.pieces.length;
out.mosaic.large_flood_fill_groups = M.pieces.filter(p => p.big).length;
out.mosaic.flood_fill_group_representatives = M.pieces.map(function (p) {
  const T = {beta: p.rep[0], s0: p.rep[1], s1: p.rep[2]};
  return {cells: p.cells, beta: p.rep[0],
          shipped: api.censusSignature("centered", T)};
});
fs.writeFileSync(path.join(__dirname, "face_widget.json"),
                 JSON.stringify(out, null, 1) + "\n");
console.log(`census checks disagreeing: ${disagree} | mosaic: ${out.mosaic.keys} keys, `
  + `${out.mosaic.flood_fill_groups} finite-grid groups `
  + `(${out.mosaic.large_flood_fill_groups} large)`);
for (const r of out.mosaic.flood_fill_group_representatives) {
  console.log("  cells", r.cells, r.shipped);
}
if (disagree) {
  console.error("check_map.js: shipped source disagrees with the centered census contract");
  process.exitCode = 1;
}
