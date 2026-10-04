"""Merge the per-worker fold-freeness subtrees into ONE certificate.

The parallelism must not cost soundness, so the merge checks the two things
that could: that the parts' index ranges tile `[0, n_parts)` exactly with no
gap and no overlap, and that every part carries the SAME spec (so they are
statements about the same region under the same partition and the same test
ordering).  Only then are their bit streams concatenated in cell-index order,
which is exactly the order the checker walks.

A face is complete iff the partition covers, every part's ranges tile, and no
subtree contains an undecided leaf -- the last of which the checker decides,
not this script.

Run: `python3 foldfree_merge.py <tag>`
"""
import base64
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    tag = sys.argv[1]
    parts = sorted(glob.glob(os.path.join(HERE, "foldfree_%s_p*.json" % tag)))
    if not parts:
        print("no parts found for tag %s" % tag)
        return 1
    docs = [json.load(open(p)) for p in parts]
    spec = docs[0]["spec"]
    for d in docs[1:]:
        if d["spec"] != spec:
            print("PARTS DISAGREE ON THE SPEC -- refusing to merge")
            return 1
    docs.sort(key=lambda d: d["lo_idx"])
    nxt = 0
    for d in docs:
        if d["lo_idx"] != nxt:
            print("PARTITION GAP OR OVERLAP: expected cell %d, part starts at "
                  "%d" % (nxt, d["lo_idx"]))
            return 1
        nxt = d["hi_idx"]
    if nxt != spec["n_parts"]:
        print("PARTITION INCOMPLETE: parts cover [0,%d) of %d cells"
              % (nxt, spec["n_parts"]))
        return 1
    bits = []
    for d in docs:
        raw = base64.b64decode(d["bits_b64"])
        b = []
        for byte in raw:
            for k in range(7, -1, -1):
                b.append((byte >> k) & 1)
        bits.extend(b[:d["n_nodes"]])
    pad = (-len(bits)) % 8
    packed = bytearray()
    for i in range(0, len(bits) + pad, 8):
        byte = 0
        for k in range(8):
            byte = (byte << 1) | (bits[i + k] if i + k < len(bits) else 0)
        packed.append(byte)
    out = os.path.join(HERE, "foldfree_%s.json" % tag)
    with open(out, "w") as fh:
        json.dump({"spec": spec, "n_nodes": len(bits),
                   "bits_b64": base64.b64encode(bytes(packed)).decode()},
                  fh, separators=(",", ":"))
    print("merged %d parts covering all %d cells -> %s (%d nodes, %.1f KB)"
          % (len(docs), spec["n_parts"], os.path.basename(out), len(bits),
             len(packed) / 1024.0))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
