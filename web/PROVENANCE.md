# Provenance

Extracted on 2026-10-03 from the research repository
`JayPiZimmermann/Skip-Connections-Avoid-Spurious-Local-Minima` at commit
`2a4c2eac24e8397625fc0c9599b35707fe958f96`:

| here | there |
|---|---|
| `web/explorer.js` | `website/widgets.js`, the shared SVG helpers (lines 7–66) and `initLandscapeExplorer` (lines 149–2653), unchanged |
| `web/census-map.js` | `website/census-map.js` (generated asset, depth 8) |
| `web/precompute/precompute_census_map.js` | `website/precompute_census_map.js`; the only edits point it at `../explorer.js` and `../census-map.js` |
| `web/precompute/check_classification.js` | `website/check_classification.js`, same repointing |
| `web/precompute/test_precompute_census_map.js` | `website/tests/test_precompute_census_map.js`, same repointing; the synthetic asset in the test gained the `finiteGridComponents` field the checker requires |
| `web/index.html` | new; the explorer markup follows `website/components/landscape-explorer.fragment`, the theorem statements follow `website/sections/05-plain-relu.md` and `website/sections/appendix-c.md` |
| `web/vendor/katex` | `website/vendor/katex` (MIT, see its LICENSE) |
| `certificates/` | `certificates/census_map_centered`, `certificates/census_map_noncentered`, `certificates/COMPLETENESS.md`; files above 3 MB omitted, see `certificates/LARGE_ARTIFACTS.md` |

## Regenerating and checking

```bash
cd web/precompute
node check_classification.js                 # no critical point listed twice
node precompute_census_map.js --check        # census-map.js matches explorer.js
node precompute_census_map.js --depth 8      # regenerate (minutes, all cores)
node test_precompute_census_map.js           # publication tests of the writer
```

## Previewing

```bash
python3 -m http.server --directory web 8765   # then open http://127.0.0.1:8765/
```
