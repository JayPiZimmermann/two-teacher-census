// Precompute the two-student CENSUS MAP shown in the landscape explorer.
//
//   node precompute_census_map.js            regenerate census-map.js
//   node precompute_census_map.js --check    read-only validation of census-map.js
//   node precompute_census_map.js --depth N  refinement depth (default 7)
//
// WHAT THE GRID COLORS AND GROUPS ARE.  A sampled cell is colored by the census
// the TABLE shows, with the SIGN of the student weights kept: for a teacher it
// is the multiset of
//   <family kind> : global | trap@positive | trap@mixed
// over the classification's rows.  An all-positive trap and an outer-split trap
// are different phenomena and only the first is visible to a positive-mass
// search, so they must not be merged.  Saddles are dropped: measured, every
// numerically unstable distinction in the full census was a saddle distinction
// (40 sign flips down one column became 0), while the minima are stable.
//
// WHY A QUADTREE.  Sampling the census on a uniform grid does not converge --
// the finite-grid group count roughly doubled with each refinement, because a
// sign-change root finder cannot see a tangential root and reports one family on and off at
// the sampling scale.  Refining only where neighbours disagree spends the work
// on boundaries, so edges are located to the depth tolerance while interiors
// stay cheap, and the finite-grid output can be checked across depths.  None of
// these group counts is asserted to be an exact connected-component count.
//
// The generator drives the SHIPPED ../explorer.js through a DOM shim, so the map
// cannot drift from the census it claims to show.
const fs = require("fs");
const path = require("path");
const {spawn} = require("child_process");
let api = null;

const ROOT = __dirname;
const WIDGETS = path.join(ROOT, "..", "explorer.js");

const argDepth = process.argv.indexOf("--depth");
const argWk = process.argv.indexOf("--worker");
const DEPTH = argWk > 0 ? parseInt(process.argv[argWk + 2], 10)
  : (argDepth > 0 ? parseInt(process.argv[argDepth + 1], 10) : 7);
const argOut = process.argv.indexOf("--out");
const ASSET_PATH = argOut > 0 ? process.argv[argOut + 1] : path.join(ROOT, "..", "census-map.js");
const argModel = process.argv.indexOf("--model");
const MODELS = argModel > 0 ? [process.argv[argModel + 1]] : ["centered", "noncentered"];
const argW = process.argv.indexOf("--workers");
const WORKERS = argW > 0 ? parseInt(process.argv[argW + 1], 10)
  : Math.max(1, Math.min(24, require("os").cpus().length - 2));
const BASE_X = 24, BASE_Y = 12;      // coarse grid before refinement

class El {
  constructor(t){this.tagName=t;this.attrs={};this.children=[];this.style={};this._text="";
    this.innerHTML="";this.dataset={};this.classList={add(){},remove(){},contains(){return false}};}
  setAttribute(k,v){this.attrs[k]=v;} getAttribute(k){return this.attrs[k];} removeAttribute(){}
  appendChild(c){this.children.push(c);return c;} removeChild(c){return c;}
  addEventListener(){} removeEventListener(){}
  getBoundingClientRect(){return{width:620,height:330,left:0,top:0};}
  querySelector(){return null;} querySelectorAll(){return[];}
  get firstChild(){return this.children.length?this.children[0]:null;}
  set textContent(v){this._text=v;this.children=[];} get textContent(){return this._text;}
}

function loadClassification() {
  const byId={}, mk=(id,x)=>(byId[id]=Object.assign(new El("div"),x||{}));
  ["ex-canvas","ex-map","ex-crit-body","ex-legend","ex-map-legend","ex-trap-note",
   "ex-crit-table","ex-trap-row","ex-trap"].forEach(i=>mk(i));
  mk("ex-regime",{value:"centered"});
  mk("t-beta0",{value:"0.00"}); mk("t-beta1",{value:"0.31"});
  mk("t-s0",{value:"1"});       mk("t-s1",{value:"1"});
  ["t-beta0","t-beta1","t-s0","t-s1"].forEach(i=>mk(i+"-val"));
  global.document={createElementNS:(n,t)=>new El(t),createElement:t=>new El(t),
    getElementById:i=>byId[i]||null,querySelectorAll:()=>[],querySelector:()=>null,
    addEventListener:()=>{},documentElement:new El("html"),body:new El("body"),head:new El("head")};
  global.window={addEventListener:()=>{},devicePixelRatio:1,
    matchMedia:()=>({matches:false,addEventListener(){}})};
  global.requestAnimationFrame=f=>setTimeout(f,0);
  global.getComputedStyle=()=>({getPropertyValue:()=>"#333"});
  globalThis.__censusExport=x=>{api=x;};
  new Function(fs.readFileSync(WIDGETS,"utf8"))();
  if(!api) throw new Error("explorer.js did not export the classification hook");
  return api;
}

// ---------------------------------------------------------------------------
// A CERTIFIED CURVE AS A SEPARATOR
//
// One census boundary of the centered map is invisible to any grid.  The
// potential curve  Wpot(beta, t) = phi(t) h(t-beta) - phi(t-beta) h(t) = 0
// runs the full width of 0 < beta < pi.  The two certified open-face bands
// 0 < y < y_lower and y_upper < y < 1 carry the SAME census -- one same-sign
// collapsed trap.  An unlabelled finite-grid flood fill coalesces their samples
// around the unresolved endpoint pinch.  This is a renderer statement, not a
// claim about connected components of an exact census level set.  Their
// two branches pinch together LINEARLY as beta -> 0,
//     y_upper - y_lower  ->  0.28113816169 * beta,
// while the first column of the grid sits at beta = per/(2 NX) and the cell
// height is 2/NY -- both proportional to 2^-depth.  The pinch is therefore
// 0.110403 cells wide at EVERY depth, and no refinement resolves it: the flood
// fill walks around the tip at any resolution.  So the curve is injected as a
// presentation label rather than sampled.
//
// It enters as a per-cell SIDE LABEL -- below the lower branch, between the
// branches, above the upper branch -- and two neighbouring cells are joined
// only when their census AND their side agree.  Blocking only the VERTICAL
// edges the curve crosses would leak: where the curve descends by a row between
// two columns, the horizontal step along that row crosses it.
//
// Only a curve that spans the whole width has this presentation contract.  The other certified
// curve of the centered map, the torque lens Wtau = 0, lives over
// beta in [1.42092547..., 1.72066717...], strictly inside (0, pi), so the
// generator has no certified full-width labelling convention for it.  Injecting
// a partial separator would split a finite-grid flood-fill group without such a
// convention.  `separatorSides` enforces this: a curve with no branch in some
// column is an error, not a partial separator.
const PI = Math.PI;
function modPI(x) { const r = x % PI; return r < 0 ? r + PI : r; }
function couplingH(t) { const x = modPI(t); return (PI / 2 - x) * Math.sin(x); }

// On the smooth piece that carries its roots, Wpot reduces to
//     Psi(b, w) = (b^2/4 - w^2) sin b + (b/2) (cos 2w + cos b),
// even in w, with  Psi(b, 0) > 0,  Psi(b, W) = (pi/4)(2b - pi) sin b < 0  for
// b < pi/2, and  Psi_w = -(2w sin b + b sin 2w) < 0  on (0, W], W = (pi - b)/2.
// So Psi has exactly one zero in (0, W) and bisection on that bracket converges
// to it; the two zeros +-w give the two branches.  The piece with b < pi/2 is
//     beta < pi/2:  b = beta,       t = (beta + pi)/2 + w
//     beta > pi/2:  b = pi - beta,  t = beta/2 + w
// and the root is carried to the map's vertical coordinate through the SHIPPED
// ratioCoord, the same one the sampler inverts.  Checked against the 200-bit
// interval enclosures of certificates/census_map_centered/cert_core.py
// (`curve_branches(beta, "potential")`): the two agree to 15 significant digits
// at every beta tested.
function reducedPsi(b, w) {
  return (b * b / 4 - w * w) * Math.sin(b) + (b / 2) * (Math.cos(2 * w) + Math.cos(b));
}
function potentialBranches(beta) {
  const lower = beta < PI / 2;
  const b = lower ? beta : PI - beta, W = (PI - b) / 2;
  // no certified bracket: beta = 0, pi/2, pi, the degenerate columns
  if (!(b > 0) || !(reducedPsi(b, 0) > 0) || !(reducedPsi(b, W) < 0)) return null;
  let lo = 0, hi = W;
  for (let k = 0; k < 200; k += 1) {
    const m = 0.5 * (lo + hi);
    if (m === lo || m === hi) break;
    if (reducedPsi(b, m) > 0) lo = m; else hi = m;
  }
  const w = 0.5 * (lo + hi);
  const ys = [w, -w].map(function (s) {
    const t = lower ? (beta + PI) / 2 + s : beta / 2 + s;
    return api.ratioCoord(-couplingH(t - beta), couplingH(t));
  });
  return ys[0] < ys[1] ? ys : [ys[1], ys[0]];
}

// Wired for the CENTERED model only.  The mechanism is model-general, but the
// noncentered map has a different boundary set, and no curve of it has been
// certified to span its width, so it keeps the plain census flood fill.
const SEPARATOR_CURVES = {
  centered: [{name: "potential curve (Wpot = 0)", branchesAt: potentialBranches}]
};

// The y-hull of each branch over a range of gaps.  Both branches are monotone
// in beta on either side of per/2 and turn there, so the two ends plus that
// turning point give the hull exactly; the nudge keeps the samples off the
// three degenerate gaps 0, per/2, per, where the reduction has no bracket.
function branchSpan(curve, b0, b1, per) {
  const E = 1e-12, half = per / 2, bs = [b0 + E, b1 - E];
  if (b0 < half && half < b1) bs.push(half - E, half + E);
  const lo = [Infinity, Infinity], hi = [-Infinity, -Infinity];
  for (const b of bs) {
    const br = curve.branchesAt(b);
    if (!br) continue;
    for (let k = 0; k < 2; k += 1) {
      if (br[k] < lo[k]) lo[k] = br[k];
      if (br[k] > hi[k]) hi[k] = br[k];
    }
  }
  return [lo, hi];
}

function separatorSides(model, per, NX, NY) {
  const curves = SEPARATOR_CURVES[model];
  if (!curves || !curves.length) return null;
  const side = new Int8Array(NX * NY);   // one base-3 digit per curve
  curves.forEach(function (curve, c) {
    const scale = Math.pow(3, c);
    for (let x = 0; x < NX; x += 1) {
      const beta = per * (x + 0.5) / NX, br = curve.branchesAt(beta);
      if (!br) throw new Error("separator " + curve.name + " has no branch at column "
        + x + " (beta = " + beta + "): a curve that does not span the width leaves its "
        + "visual-pairing contract undefined and must not be used as a separator");
      for (let y = 0; y < NY; y += 1) {
        const yy = 2 * (y + 0.5) / NY - 1;
        side[y * NX + x] += scale * (yy < br[0] ? 0 : (yy > br[1] ? 2 : 1));
      }
    }
  });
  return side;
}

// One column block of the quadtree, sampled independently of every other: the
// recursion never reads outside its own block, which is what lets the blocks be
// farmed out to separate processes.  Returns run-length rows over the block's
// own columns, against the WORKER's local key table -- the parent remaps.
function blockRuns(api, model, blockIndices) {
  const per = api.kernelOf(model).per;
  const keys = [], index = {};
  const idOf = k => { if(!(k in index)){index[k]=keys.length;keys.push(k);} return index[k]; };
  const cache = new Map();
  const NX = BASE_X * (1 << DEPTH), NY = BASE_Y * (1 << DEPTH);
  const at = (bx, by) => {                       // bx, by in [0,1]
    const k = bx.toFixed(7)+","+by.toFixed(7);
    if (cache.has(k)) return cache.get(k);
    const m = api.massesAt(Math.max(-1+1e-9, Math.min(1-1e-9, 2*by-1)));
    const v = idOf(api.censusSignature(model, {beta: per*bx, s0: m.s0, s1: m.s1}));
    cache.set(k, v);
    return v;
  };
  // A block the certified separator runs through is refined even when its five
  // probes agree.  The flood fill is about to cut along that curve, so the cells
  // on either side of the cut must carry their OWN census rather than a coarse
  // block's; otherwise a cell painted from a block that straddles the curve
  // becomes an island of its own the moment the cut is made.  This does not
  // change the finite-grid flood-fill count, which only the separator can fix; it keeps the
  // painted census honest where the cut lands.  Measured: without it the
  // depth-6 and depth-7 maps grow four 2-cell islands, with it they do not.
  const curves = SEPARATOR_CURVES[model] || [];
  const onSeparator = (x0, y0, w, h) => {
    const yb = 2*y0/NY - 1, yt = 2*(y0+h)/NY - 1;
    for (const cv of curves) {
      const s = branchSpan(cv, per*x0/NX, per*(x0+w)/NX, per);
      for (let k = 0; k < 2; k += 1) if (s[0][k] <= yt && s[1][k] >= yb) return true;
    }
    return false;
  };
  const cw = NX / BASE_X, ch = NY / BASE_Y;
  const blocks = {};
  for (const i of blockIndices) {
    const sub = new Int16Array(cw * NY).fill(-1);
    const x00 = i * cw;
    const paint = (x0, y0, w, h, v) => {
      for (let y = y0; y < y0 + h; y++) sub.fill(v, y*cw + x0 - x00, y*cw + x0 - x00 + w);
    };
    const rec = (x0, y0, w, h, depth) => {
      const c = [at((x0+0.5*w)/NX, (y0+0.5*h)/NY),
                 at((x0+0.02*w)/NX, (y0+0.02*h)/NY),
                 at((x0+0.98*w)/NX, (y0+0.02*h)/NY),
                 at((x0+0.02*w)/NX, (y0+0.98*h)/NY),
                 at((x0+0.98*w)/NX, (y0+0.98*h)/NY)];
      if ((c.every(v => v === c[0]) && !onSeparator(x0,y0,w,h)) || depth === 0 || w === 1 || h === 1) {
        paint(x0,y0,w,h,c[0]); return;
      }
      const hw = w >> 1, hh = h >> 1;
      rec(x0, y0, hw, hh, depth-1);           rec(x0+hw, y0, w-hw, hh, depth-1);
      rec(x0, y0+hh, hw, h-hh, depth-1);      rec(x0+hw, y0+hh, w-hw, h-hh, depth-1);
    };
    for (let j = 0; j < BASE_Y; j++) rec(x00, j*ch, cw, ch, DEPTH);
    const rows = [];
    for (let y = 0; y < NY; y++) {
      const r = []; let x0 = 0, v0 = sub[y*cw];
      for (let x = 1; x <= cw; x++) {
        const v = x < cw ? sub[y*cw+x] : -2;
        if (v !== v0) { r.push([x0, x, v0]); x0 = x; v0 = v; }
      }
      rows.push(r);
    }
    blocks[i] = rows;
  }
  return {keys: keys, blocks: blocks, samples: cache.size};
}

// Reassemble the sharded blocks into one grid, then find the pieces.
function buildModel(model, shards) {
  const NX = BASE_X * (1 << DEPTH), NY = BASE_Y * (1 << DEPTH), cw = NX / BASE_X;
  const per = shards.per;
  const keys = [], index = {};
  const idOf = k => { if(!(k in index)){index[k]=keys.length;keys.push(k);} return index[k]; };
  const grid = new Int16Array(NX * NY).fill(-1);
  for (const sh of shards.parts) {
    const remap = sh.keys.map(idOf);
    for (const i of Object.keys(sh.blocks)) {
      const x00 = (+i) * cw, rows = sh.blocks[i];
      for (let y = 0; y < NY; y++) {
        for (const r of rows[y]) grid.fill(remap[r[2]], y*NX + x00 + r[0], y*NX + x00 + r[1]);
      }
    }
  }
  if (grid.includes(-1)) throw new Error("shards did not cover the grid");

  // pieces = finite-grid connected components of equal signature (4-neighbour
  // union-find), with the presentation separator cutting cells that a grid
  // cannot tell apart.  These are not components of an exact level set.
  const side = separatorSides(model, per, NX, NY);
  const joined = (i,j) => grid[i]===grid[j] && (side===null || side[i]===side[j]);
  const parent = new Int32Array(NX*NY); for (let i=0;i<parent.length;i++) parent[i]=i;
  const find = a => { while (parent[a]!==a) { parent[a]=parent[parent[a]]; a=parent[a]; } return a; };
  const uni = (a,b) => { a=find(a); b=find(b); if(a!==b) parent[b]=a; };
  let cut = 0;
  for (let y=0;y<NY;y++) for (let x=0;x<NX;x++) {
    const i=y*NX+x;
    if (x+1<NX && grid[i]===grid[i+1]) { if (joined(i,i+1)) uni(i,i+1); else cut++; }
    if (y+1<NY && grid[i]===grid[i+NX]) { if (joined(i,i+NX)) uni(i,i+NX); else cut++; }
  }
  const comp = new Map();
  for (let i=0;i<NX*NY;i++) {
    const r=find(i);
    if(!comp.has(r)) comp.set(r,{id:grid[i],cells:0});
    comp.get(r).cells++;
  }

  // representative = the cell FARTHEST from any differing neighbour, so a click
  // lands deep inside the piece rather than near an edge the map cannot resolve
  const dist = new Int32Array(NX*NY).fill(-1);
  const q = [];
  for (let y=0;y<NY;y++) for (let x=0;x<NX;x++) {
    const i=y*NX+x;
    let edge = x===0||y===0||x===NX-1||y===NY-1;
    if(!edge && (!joined(i,i-1)||!joined(i,i+1)||!joined(i,i-NX)||!joined(i,i+NX))) edge=true;
    if(edge){ dist[i]=0; q.push(i); }
  }
  for (let h=0; h<q.length; h++) {
    const i=q[h], x=i%NX, y=(i-x)/NX, nb=[];
    if(x>0)nb.push(i-1); if(x<NX-1)nb.push(i+1); if(y>0)nb.push(i-NX); if(y<NY-1)nb.push(i+NX);
    for (const j of nb) if (dist[j]<0 && joined(i,j)) { dist[j]=dist[i]+1; q.push(j); }
  }
  const best = new Map();
  for (let i=0;i<NX*NY;i++) {
    const r=find(i);
    if(!best.has(r)||dist[i]>best.get(r).d) best.set(r,{d:dist[i],x:i%NX,y:(i-(i%NX))/NX});
  }

  // Every finite-grid component becomes a piece, indexed, and the run-length
  // rows reference that INDEX rather than the type: two pieces of one type must stay separable
  // so a click can resolve to the piece under the pointer.  `big` marks the ones
  // worth a dot and a click target.
  const MIN_CELLS = Math.max(4, Math.round(NX*NY/40000));
  const pieces = [], pieceOf = new Map();
  for (const [r,c] of comp) {
    const b = best.get(r);
    const beta = per*(b.x+0.5)/NX, yy = 2*(b.y+0.5)/NY - 1;
    const m = api.massesAt(yy);
    pieceOf.set(r, pieces.length);
    pieces.push({id:c.id, cells:c.cells, big:c.cells >= MIN_CELLS,
                 rep:[+beta.toFixed(6), +m.s0.toFixed(6), +m.s1.toFixed(6)],
                 at:[+((b.x+0.5)/NX).toFixed(6), +((b.y+0.5)/NY).toFixed(6)]});
  }

  // Outline each piece by walking its boundary EDGES and chaining them into
  // closed loops, HERE rather than in the browser: at this refinement the grid
  // is a million cells and the walk is not something a page load should do.
  // The pieces tile, so a shared boundary is one polyline emitted twice from
  // opposite sides -- the mosaic has no seams and no lanes.
  const pgrid = new Int32Array(NX*NY);
  for (let i=0;i<NX*NY;i++) pgrid[i] = pieceOf.get(find(i));
  const W = NX+1;                                   // lattice of corner points
  const edges = new Map();                          // piece -> Map(from -> [to])
  const addEdge = (p,ax,ay,bx,by) => {
    let m = edges.get(p); if (!m) { m = new Map(); edges.set(p,m); }
    const a = ay*W+ax, b = by*W+bx;
    const l = m.get(a); if (l) l.push(b); else m.set(a,[b]);
  };
  for (let y=0;y<NY;y++) for (let x=0;x<NX;x++) {
    const i=y*NX+x, p=pgrid[i];
    // counter-clockwise in grid coordinates, so loops close consistently
    if (x===0    || pgrid[i-1]  !== p) addEdge(p, x,   y+1, x,   y  );
    if (x===NX-1 || pgrid[i+1]  !== p) addEdge(p, x+1, y,   x+1, y+1);
    if (y===0    || pgrid[i-NX] !== p) addEdge(p, x,   y,   x+1, y  );
    if (y===NY-1 || pgrid[i+NX] !== p) addEdge(p, x+1, y+1, x,   y+1);
  }
  // Drop the interior points of a straight run: the outlines are rectilinear,
  // so this cuts the path size by orders of magnitude and changes nothing.
  const simplify = loop => {
    const out=[];
    for (let i=0;i<loop.length;i++) {
      const a=loop[(i-1+loop.length)%loop.length], b=loop[i], c=loop[(i+1)%loop.length];
      if ((a[0]===b[0]&&b[0]===c[0])||(a[1]===b[1]&&b[1]===c[1])) continue;
      out.push(b);
    }
    return out.length>2 ? out : loop;
  };
  let pts = 0;
  for (const [p,m] of edges) {
    const loops = [];
    for (const start of Array.from(m.keys())) {
      while (m.get(start) && m.get(start).length) {
        const loop=[]; let cur=start;
        while (true) {
          const nexts = m.get(cur);
          if (!nexts || !nexts.length) break;
          const nxt = nexts.pop();
          loop.push([cur%W, (cur-(cur%W))/W]);
          cur = nxt;
          if (cur === start) break;
        }
        if (loop.length>2) loops.push(simplify(loop));
      }
    }
    pieces[p].loops = loops;
    loops.forEach(l => { pts += l.length; });
  }
  pieces.forEach(q => { if (!q.loops) q.loops = []; });
  return {per:per, nx:NX, ny:NY, keys:keys, pieces:pieces,
          samples:shards.samples, finiteGridComponents:comp.size, outlinePoints:pts,
          separatorCuts:cut};
}

function readAsset() {
  if (!fs.existsSync(ASSET_PATH)) return null;
  const m = fs.readFileSync(ASSET_PATH,"utf8").match(/window\.CensusMap\s*=\s*(\{[\s\S]*\});\s*$/);
  if (!m) return null;
  try { return JSON.parse(m[1]); } catch (e) { return null; }
}

// Farm the column blocks out to worker processes: a classification costs
// 8-16 ms, so a publishable depth is tens of minutes on one core and a couple
// of minutes across the box.  Blocks are interleaved (i % N === j) rather than
// contiguous, because the cost per block is very uneven -- boundaries cluster.
function runShards(model) {
  const jobs = [];
  for (let j = 0; j < WORKERS; j += 1) {
    const idx = [];
    for (let i = j; i < BASE_X; i += WORKERS) idx.push(i);
    if (idx.length) jobs.push({j: j, idx: idx});
  }
  return Promise.all(jobs.map(job => new Promise((resolve, reject) => {
    const ch = spawn(process.execPath,
      [__filename, "--worker", model, String(DEPTH), job.idx.join(",")],
      {stdio: ["ignore", "pipe", "pipe"]});
    const buf = [], err = [];
    ch.stdout.on("data", d => buf.push(d));
    ch.stderr.on("data", d => err.push(d));
    ch.on("close", code => {
      if (code !== 0) return reject(new Error(`shard ${job.j} exited ${code}: ${Buffer.concat(err)}`));
      const part = JSON.parse(Buffer.concat(buf).toString("utf8"));
      process.stdout.write(`  ${model}: shard ${job.j + 1}/${jobs.length} done `
        + `(blocks ${job.idx.join(",")}), ${part.samples} calls\n`);
      resolve(part);
    });
  }))).then(parts => ({per: api.kernelOf(model).per, parts: parts,
                       samples: parts.reduce((a, b) => a + b.samples, 0)}));
}

function worker() {
  const k = process.argv.indexOf("--worker");
  const model = process.argv[k+1];
  const idx = process.argv[k+3].split(",").map(Number);
  loadClassification();
  const r = blockRuns(api, model, idx);
  process.stdout.write(JSON.stringify(r));
}

async function generate() {
  loadClassification();
  const asset = {depth: DEPTH, models: {}};
  for (const model of MODELS) {
    const d = buildModel(model, await runShards(model));
    asset.models[model] = d;
    const big = d.pieces.filter(p=>p.big).length;
    console.log(`  ${model}: ${d.keys.length} types, ${big} displayed groups `
      + `(of ${d.pieces.length} finite-grid components), `
      + `${d.samples} classification calls, grid ${d.nx}x${d.ny}, ${d.outlinePoints} outline points, `
      + `${d.separatorCuts} same-census joins cut by the certified separator`);
    d.keys.forEach((k,i)=>console.log(`     [${i}] ${k||"(nothing)"}`));
  }
  const text =
    "// GENERATED by precompute_census_map.js -- do not edit.\n" +
    "// Finite-grid flood-fill rendering of the two-student census (minima only,\n" +
    "// student-weight sign kept), with one representative per rendered group.\n" +
    "// Regenerate with:\n" +
    "//   node precompute_census_map.js\n" +
    "window.CensusMap = " + JSON.stringify(asset) + ";\n";
  writeAsset(text);
  console.log(`wrote census-map.js (${text.length.toLocaleString("en-US")} bytes, depth ${DEPTH})`);
}

function writeAsset(text, assetPath = ASSET_PATH, fileSystem = fs) {
  // Keep the temporary file beside its destination: rename is then an atomic
  // replacement on the destination filesystem, and any failed write leaves
  // the previously published asset untouched.
  const directory = path.dirname(assetPath);
  const base = path.basename(assetPath);
  let temporary = null;
  let descriptor = null;
  try {
    for (let attempt = 0; attempt < 100; attempt += 1) {
      temporary = path.join(directory,
        `.${base}.stage-${process.pid}-${Date.now().toString(36)}-${attempt}`);
      try {
        descriptor = fileSystem.openSync(temporary, "wx", 0o600);
        break;
      } catch (error) {
        if (error && error.code === "EEXIST") continue;
        throw error;
      }
    }
    if (descriptor === null) throw new Error(`could not reserve temporary asset beside ${assetPath}`);
    fileSystem.writeFileSync(descriptor, text, "utf8");
    fileSystem.fsyncSync(descriptor);
    fileSystem.closeSync(descriptor);
    descriptor = null;
    fileSystem.renameSync(temporary, assetPath);
    temporary = null;
  } finally {
    if (descriptor !== null) fileSystem.closeSync(descriptor);
    if (temporary !== null) {
      try { fileSystem.unlinkSync(temporary); }
      catch (error) { if (!error || error.code !== "ENOENT") throw error; }
    }
  }
}

function check() {
  const a = readAsset();
  if (!a) { console.error("census-map.js: missing or unparsable; run node precompute_census_map.js"); return 1; }
  for (const m of ["centered","noncentered"]) {
    const d = a.models && a.models[m];
    if (!d || !Array.isArray(d.pieces) || !d.pieces.every(q => Array.isArray(q.loops))) {
      console.error("census-map.js: model "+m+" missing or malformed"); return 1;
    }
    if (!Number.isInteger(d.finiteGridComponents)
        || d.finiteGridComponents !== d.pieces.length
        || Object.prototype.hasOwnProperty.call(d, "components")) {
      console.error("census-map.js: model "+m
        + " must label its renderer count as finiteGridComponents"); return 1;
    }
  }
  console.log("census-map.js: present and well formed");
  return 0;
}

if (require.main === module) {
  if (process.argv.includes("--worker")) { worker(); }
  else if (process.argv.includes("--check")) { process.exit(check()); }
  else { generate().catch(e => { console.error(e); process.exit(1); }); }
}

module.exports = {check, readAsset, writeAsset};
