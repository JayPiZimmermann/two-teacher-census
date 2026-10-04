// DOM shim loading the SHIPPED website/widgets.js (read only) and returning its
// classification export.  Same shim as website/precompute_census_map.js, so the
// certificate compares against exactly the code the page runs.
const fs = require("fs");
const path = require("path");

const WIDGETS = path.join(__dirname, "..", "..", "website", "widgets.js");

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
  let api = null;
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
  if(!api) throw new Error("widgets.js did not export the classification hook");
  return api;
}

module.exports = {loadClassification: loadClassification};
