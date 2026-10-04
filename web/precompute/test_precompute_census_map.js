// Focused publication tests: no census computation or site build is run here.
const assert = require("assert");
const childProcess = require("child_process");
const fs = require("fs");
const os = require("os");
const path = require("path");

const generator = require("./precompute_census_map.js");
const temporary = fs.mkdtempSync(path.join(os.tmpdir(), "census-map-publication-"));
const asset = path.join(temporary, "census-map.js");

try {
  fs.writeFileSync(asset, "old asset\n", "utf8");
  const failingFs = Object.assign({}, fs, {
    renameSync() { throw new Error("injected rename failure"); }
  });
  assert.throws(() => generator.writeAsset("new asset\n", asset, failingFs),
                /injected rename failure/);
  assert.strictEqual(fs.readFileSync(asset, "utf8"), "old asset\n");
  assert.deepStrictEqual(fs.readdirSync(temporary), ["census-map.js"]);

  generator.writeAsset("new asset\n", asset);
  assert.strictEqual(fs.readFileSync(asset, "utf8"), "new asset\n");
  assert.deepStrictEqual(fs.readdirSync(temporary), ["census-map.js"]);

  const validAsset = "window.CensusMap = " + JSON.stringify({
    models: {centered: {pieces: [{loops: []}], finiteGridComponents: 1},
             noncentered: {pieces: [{loops: []}], finiteGridComponents: 1}}
  }) + ";\n";
  fs.writeFileSync(asset, validAsset, "utf8");
  const before = fs.readFileSync(asset, "utf8");
  const result = childProcess.spawnSync(process.execPath,
    [path.join(__dirname, "precompute_census_map.js"), "--check", "--out", asset],
    {encoding: "utf8"});
  assert.strictEqual(result.status, 0, result.stderr);
  assert.strictEqual(fs.readFileSync(asset, "utf8"), before);
  assert.deepStrictEqual(fs.readdirSync(temporary), ["census-map.js"]);
} finally {
  fs.rmSync(temporary, {recursive: true, force: true});
}
