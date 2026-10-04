"""Replay a ``minimum_dfs.py`` selected-minimum certificate.

The checker trusts only the spec and bitstream.  It regenerates the root
partition and every child box, re-derives each stored interval proof method,
verifies structural coverage, rejects leftover/truncated bits, and treats
every ``UNDECIDED`` leaf as an incomplete wall candidate.  Proof methods are
checked directly rather than through the search priority order, so adding a
new earlier test cannot invalidate a still-correct older leaf.

Usage: python3 minimum_dfs_check.py minimum_dfs_<tag>_pNNN.json
"""
import base64
import json
import os
import sys
import time

import minimum_dfs as D

HERE = os.path.dirname(os.path.abspath(__file__))


class Reader:
    def __init__(self, encoded, n_bits):
        self.data = base64.b64decode(encoded)
        self.n_bits = n_bits
        self.pos = 0

    def read(self, width=1):
        value = 0
        for _ in range(width):
            if self.pos >= self.n_bits:
                raise ValueError("bitstream exhausted")
            byte = self.data[self.pos >> 3]
            value = (value << 1) | ((byte >> (7 - (self.pos & 7))) & 1)
            self.pos += 1
        return value


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    path = sys.argv[1]
    if not os.path.isabs(path):
        path = os.path.join(HERE, path)
    with open(path) as handle:
        doc = json.load(handle)
    spec = doc["spec"]
    lo_idx, hi_idx = doc["lo_idx"], doc["hi_idx"]
    if not (D.partition_covers(spec) and 0 <= lo_idx < hi_idx
            and hi_idx <= spec["n_parts"] and 0 < spec["delta"]
            and spec["delta"] < spec["dmax"] < D.mpmath.pi):
        print("REPLAY FAILED: inadmissible spec or partition")
        return 1
    reader = Reader(doc["bits_b64"], doc["n_nodes"])
    roots = D.root_cells(spec)
    counts = {name: 0 for name in D.VERDICTS}
    bad = []
    wall_boxes = []
    leaf_volume = 0.0
    leaves_seen = 0
    started = time.time()

    def walk(box, depth=0):
        nonlocal leaf_volume, leaves_seen
        if depth > 500:
            raise ValueError("tree depth exceeds 500")
        if reader.read() == 1:
            walk(D.child(spec, box, 0), depth + 1)
            walk(D.child(spec, box, 1), depth + 1)
            return
        code = reader.read(3)
        claimed = D.NAME.get(code)
        if claimed is None:
            bad.append(("unknown verdict code", code, box))
            return
        # Re-derive the stored proof method directly.  Search priority is not
        # part of the mathematical claim and can change when a sharper test is
        # added ahead of an older one.
        got = D.verdict(spec, *box, target=claimed)
        if claimed == "UNDECIDED":
            width = max(box[2 * i + 1] - box[2 * i] for i in range(4))
            valid = got is None and width < spec["minw"]
            if valid and len(wall_boxes) < 50:
                wall_boxes.append(box)
        else:
            valid = got == claimed
        if not valid and len(bad) < 30:
            bad.append((claimed, got, box))
        counts[claimed] += 1
        leaves_seen += 1
        leaf_volume += D.volume(box)
        if leaves_seen % 5000 == 0:
            print("  replayed %7d leaves, %7d/%d bits, %.0fs" %
                  (leaves_seen, reader.pos, reader.n_bits,
                   time.time() - started), flush=True)

    try:
        for idx in range(lo_idx, hi_idx):
            walk(roots[idx])
    except (ValueError, IndexError) as exc:
        print("REPLAY FAILED:", exc)
        return 1
    expected_volume = sum(D.volume(roots[i]) for i in range(lo_idx, hi_idx))
    tolerance = 2e-12 * max(1.0, abs(expected_volume))
    print("certificate :", os.path.basename(path))
    print("cells       : [%d,%d) of %d" %
          (lo_idx, hi_idx, spec["n_parts"]))
    print("bits        : %d claimed, %d consumed" %
          (doc["n_nodes"], reader.pos))
    print("leaves      :", json.dumps(counts, sort_keys=True))
    print("coverage    : %.16e leaf volume / %.16e root volume" %
          (leaf_volume, expected_volume))
    if wall_boxes:
        print("wall boxes  : first %d replayed candidates" % len(wall_boxes))
        for box in wall_boxes[:10]:
            print("  beta=[%.12g,%.12g] y=[%.12g,%.12g] "
                  "s=[%.12g,%.12g] D=[%.12g,%.12g]" % box)
    if reader.pos != doc["n_nodes"]:
        bad.append(("node count mismatch", reader.pos, doc["n_nodes"]))
    if abs(leaf_volume - expected_volume) > tolerance:
        bad.append(("coverage mismatch", leaf_volume, expected_volume))
    if bad:
        print("REPLAY FAILED: %d mismatches" % len(bad))
        for item in bad[:10]:
            print(" ", item)
        return 1
    if counts["UNDECIDED"]:
        print("CERTIFICATE VALID BUT INCOMPLETE: %d PSD-wall candidates"
              % counts["UNDECIDED"])
        return 1
    if spec.get("scope") == "census":
        # NEG_T00 is not admissible here, so no zero was set aside for being
        # a saddle: this is the fold-freeness `hfold` of
        # `censusAngleMapJ_zeroCount_eq_witness`, not the selected-minimum
        # regularity, and the two must not read alike in a log.
        print("CERTIFICATE VALID AND COMPLETE: every leaf recomputes and "
              "EVERY angle-map zero has nonzero Schur determinant throughout "
              "the named region; no saddle fold is admitted.")
    else:
        print("CERTIFICATE VALID AND COMPLETE: every leaf recomputes and the "
              "selected-root determinant is nonzero throughout the named "
              "region; unselected saddle folds are allowed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
