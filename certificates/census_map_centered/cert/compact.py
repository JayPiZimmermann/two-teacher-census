import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
"""Compact the raw certificate into a distributable certificate.json.

The raw run emits every certified curve box (~92k) and every wide box (~102k),
which is 263 MB.  Nothing is discarded silently: the full run is kept
gzip-compressed as certificate_full.json.gz, and the compact file keeps
  - all summary statistics computed over the FULL set (measures, hulls, counts),
  - an evenly spaced sample of the curve boxes,
  - a summary plus a sample of the wide/unresolved boxes.
"""
import gzip
import json
import shutil

D = _HERE + "/"
SAMPLE_CURVE = 400
SAMPLE_WIDE = 60

# The raw certifier records paths into its full ``boxes`` arrays.  Compaction
# replaces those arrays by samples, so keep the human-facing ledger and the
# distributable-file inventory synchronized here rather than relying on a
# one-off edit to a previously compacted artifact.
FILES = {
    "centered.py": "kernel module: float mirror of widgets.js + mpmath.iv interval forms",
    "ref_widgets.js": "reference generator: literal widgets.js functions (identifier-extracted, guarded) -> ref_samples.txt",
    "ref_libm.js": "libm attribution harness (node's own sin/cos per sample)",
    "validate.py": "validation: widgets.js mirror, libm attribution, reduction proof (sympy + intervals), linear identity",
    "validation.json": "the latest validation run (regenerate: node ref_widgets.js && python3 validate.py)",
    "certify.py": "the certification",
    "compact.py": "downsamples the raw output into certificate.json",
    "make_md.py": "regenerates CERTIFICATE.md from certificate.json + validation.json",
    "certificate.json": "the compact machine-readable certificate (this file)",
    "certificate_full.json.gz": "the complete run, every certified and wide box (not in version control; certify.py regenerates it)",
    "CERTIFICATE.md": "human-readable appendix",
    "run.log": "console log of the certification run",
}

with open(D + "certificate.json") as fh:
    c = json.load(fh)


def sample(lst, k):
    if len(lst) <= k:
        return list(lst), len(lst)
    step = len(lst) / float(k)
    return [lst[int(i * step)] for i in range(k)], len(lst)


def measure(lst, key="beta_box"):
    tot = 0.0
    for r in lst:
        b = r[key]
        tot += float(b[1]) - float(b[0])
    return tot


for name in ("torque_lens", "potential_curve"):
    cur = c["curves"][name]
    boxes = cur.pop("boxes")
    wide = cur.pop("wide_boxes")
    s, n = sample(boxes, SAMPLE_CURVE)
    cur["n_certified_boxes"] = n
    cur["certified_boxes_sample"] = s
    cur["certified_boxes_sample_note"] = (
        "evenly spaced sample of the %d certified boxes; the summary numbers "
        "above are computed over all of them.  The complete list is in "
        "certificate_full.json.gz" % n)
    ws, wn = sample(wide, SAMPLE_WIDE)
    cur["wide_boxes_summary"] = {
        "n": wn,
        "beta_measure": measure(wide),
        "max_width": max([r["width"] for r in wide]) if wide else 0.0,
        "reasons": {r: sum(1 for x in wide if x.get("reason") == r)
                    for r in sorted(set(x.get("reason") for x in wide))},
        "sample": ws}

for piece in c["pieces_1D"]:
    if piece["id"] in ("1D-1", "1D-2"):
        piece["certified_boxes"] = (
            "certificate.json: curves.torque_lens.certified_boxes_sample "
            "(full list: certificate_full.json.gz)")
    elif piece["id"] in ("1D-3", "1D-4"):
        piece["certified_boxes"] = (
            "certificate.json: curves.potential_curve.certified_boxes_sample "
            "(full list: certificate_full.json.gz)")
    elif piece["id"] == "1D-5":
        piece["note"] = (
            "coincides with the degenerate stratum 1D-9 (second teacher "
            "massless).  Here s1 = h(0) = 0 exactly and s0 != 0, so atan2 "
            "gives +-pi/2 and the fold gives y = 0; the per-beta table reports "
            "an enclosure of width < 1e-27 around 0.")
    elif piece["id"] == "1D-6":
        piece["note"] = (
            "coincides with the degenerate stratum 1D-10.  In the per-beta "
            "table the y field of this root is reported as null: s0 = -h(0) "
            "= 0 EXACTLY, so atan2(s0,s1) sits on the seam of the (-,-) fold "
            "and the interval routine correctly declines to pick a branch.  "
            "The value y = +1 follows from the exact vanishing of s0 together "
            "with ratioCoord's convention a <= 0 -> a + pi, not from a "
            "numerical evaluation.")

c["curves"]["weight_wronskian_image"]["per_beta_table_evidence"] = (
    "every row of per_beta_table has Wwgt with n_interior_roots = 0 and "
    "exactly two kink roots t = 0 and t = beta (mod pi)")
c["files"] = FILES

# unresolved: keep the beta collars in full, summarise the rest
un = c["unresolved_boxes"]
collars = [u for u in un if u["kind"] == "beta_collar_uncertified"]
rest = [u for u in un if u["kind"] != "beta_collar_uncertified"]
kinds = sorted(set(u["kind"] for u in rest))
c["unresolved_boxes"] = {
    "beta_collars_uncertified": collars,
    "other_kinds": {k: {"n": sum(1 for u in rest if u["kind"] == k),
                        "max_width": max([u["width"] for u in rest
                                          if u["kind"] == k] or [0.0]),
                        "total_beta_measure": sum(u["width"] for u in rest
                                                  if u["kind"] == k),
                        "sample": sample([u for u in rest if u["kind"] == k],
                                         SAMPLE_WIDE)[0]}
                    for k in kinds},
    "total_entries": len(un)}

with open(D + "certificate_compact.json", "w") as fh:
    json.dump(c, fh, indent=1)
print("compact written")

# archive the full one
with open(D + "certificate.json", "rb") as fi, gzip.open(
        D + "certificate_full.json.gz", "wb", compresslevel=6) as fo:
    shutil.copyfileobj(fi, fo, 1 << 22)
shutil.move(D + "certificate_compact.json", D + "certificate.json")
print("done")
