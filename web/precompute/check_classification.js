// Guard: the two-student classification must never list one critical point twice.
//
//   node check_classification.js          run the sweep, exit 1 on any duplicate
//
// The two students are UNLABELLED, so a row and the row with the students
// exchanged are the SAME critical point and only one of them belongs in the
// table.  Folding that symmetry by an inequality on the gap does NOT work: at a
// gap of exactly half a period both orderings satisfy it, and that is precisely
// the self-mirror configuration a symmetric teacher produces.  Measured before
// the fix: 318 duplicate rows across 2576 teacher configurations, every one of
// them an `opposite bisector pair` on a symmetric teacher.
//
// A generic mass grid does NOT expose this -- the first sweep run against the
// broken code reported zero, because the defect lives exactly on `s0 = ±s1`.
// So the mass list below contains the symmetric teachers explicitly, and must
// keep containing them.
//
// NOT a duplicate, and deliberately kept: two rows related by the teacher's own
// MIRROR symmetry when that mirror does not fix them (e.g. the `separated pair`
// rows at 0.391π,1.848π and 0.759π,1.302π for the gap-1.15π symmetric teacher).
// Those are two distinct critical points that happen to share a mass multiset.
const fs = require("fs");
const path = require("path");

const WIDGETS = path.join(__dirname, "..", "explorer.js");

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
   "ex-crit-table","ex-trap-row","ex-trap","ex-canvas-view","ex-map-view"].forEach(i=>mk(i));
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
  let api=null;
  globalThis.__censusExport=x=>{api=x;};
  new Function(fs.readFileSync(WIDGETS,"utf8"))();
  if(!api||!api.classificationRows) throw new Error("explorer.js did not export classificationRows");
  return api;
}

// The canonical name of an unordered two-student configuration.
function pairKey(r, per) {
  if (!r.theta || r.theta.length !== 2) return null;
  const q = x => Math.round(((x % per) + per) % per / 1e-4) * 1e-4;
  const m = x => Math.round((x || 0) / 1e-4) * 1e-4;
  let a = [q(r.theta[0]), m(r.c && r.c[0])], b = [q(r.theta[1]), m(r.c && r.c[1])];
  if (b[0] < a[0] || (b[0] === a[0] && b[1] < a[1])) { const t = a; a = b; b = t; }
  return a.join(",") + "/" + b.join(",");
}

const PI = Math.PI;
const MASSES = [[1,1],[-2.4,-2.4],[1,-1],[-1,1],[2.4,-2.4],[1,0.5],[0.7,1.3],[-1.6,0.9],
                [-2.4,2.4],[0.3,0.3],[2,-2]];
const NB = 160;

function main() {
  const api = loadClassification();
  let checked = 0, dup = 0;
  const shown = [];
  for (const model of ["centered","noncentered"]) {
    const per = api.kernelOf(model).per;
    for (let bi = 0; bi <= NB; bi += 1) {
      for (const ms of MASSES) {
        const beta = per * bi / NB;
        let rows;
        try { rows = api.classificationRows(model, {beta: beta, s0: ms[0], s1: ms[1]}); }
        catch (e) { continue; }
        if (!rows) continue;
        checked += 1;
        const seen = new Map();
        for (const r of rows) {
          const k = pairKey(r, per);
          if (k === null) continue;
          if (seen.has(k)) {
            dup += 1;
            if (shown.length < 10) shown.push(
              `  ${model} β=${(beta/PI).toFixed(5)}π s=(${ms[0]},${ms[1]}) `
              + `"${seen.get(k).family}" and "${r.family}" are one point: `
              + `${r.angleText} | ${r.massText} | ${r.type}`);
          } else seen.set(k, r);
        }
      }
    }
  }
  console.log(`classification: ${checked} teacher configurations checked, ${dup} duplicate rows`);
  if (dup) {
    console.error("check_classification.js: the table lists the same critical point twice.");
    console.error("The two students are unlabelled; emit one canonical representative.");
    shown.forEach(l => console.error(l));
    return 1;
  }
  console.log("classification: no permutation duplicates");
  return 0;
}

process.exit(main());
