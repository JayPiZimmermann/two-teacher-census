"""Replay a selected-zero boundary-exclusion certificate."""
import base64
import json
import os
import sys

import minimum_boundary_dfs as B

HERE = os.path.dirname(os.path.abspath(__file__))


class Reader:
    def __init__(self, encoded, n_bits):
        self.data = base64.b64decode(encoded)
        self.n_bits = n_bits
        self.pos = 0

    def read(self):
        if self.pos >= self.n_bits:
            raise ValueError("bitstream exhausted")
        byte = self.data[self.pos >> 3]
        value = (byte >> (7 - (self.pos & 7))) & 1
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
    if not (B.partition_covers(spec) and 0 <= lo_idx < hi_idx
            and hi_idx <= spec["n_parts"]):
        print("REPLAY FAILED: inadmissible spec")
        return 1
    reader = Reader(doc["bits_b64"], doc["n_nodes"])
    roots = B.root_cells(spec)
    leaves = 0
    bad = []

    def walk(box, depth=0):
        nonlocal leaves
        if depth > 500:
            raise ValueError("tree depth exceeds 500")
        if reader.read():
            walk(B.child(spec, box, 0), depth + 1)
            walk(B.child(spec, box, 1), depth + 1)
            return
        leaves += 1
        if not B.excludes(spec, box) and len(bad) < 20:
            bad.append(box)
        if leaves % 5000 == 0:
            print("  replayed %d boundary leaves" % leaves, flush=True)

    try:
        for idx in range(lo_idx, hi_idx):
            walk(roots[idx])
    except (ValueError, IndexError) as exc:
        print("REPLAY FAILED:", exc)
        return 1
    print("certificate :", os.path.basename(path))
    print("face        :", spec["face"])
    print("bits        : %d claimed, %d consumed" %
          (doc["n_nodes"], reader.pos))
    print("leaves      :", leaves)
    if reader.pos != doc["n_nodes"] or bad:
        print("REPLAY FAILED: %d invalid leaves" % len(bad))
        for box in bad[:5]:
            print(" ", box)
        return 1
    # Report what this artifact's scope actually proves.  A census forest
    # admits map exclusions only, so it establishes the strictly stronger
    # statement, and a message naming the selected set would understate it.
    if spec.get("scope") == "census":
        print("CERTIFICATE VALID AND COMPLETE: NO zero of the angle map lies "
              "on this physical census-strip face throughout the named "
              "teacher region.")
    else:
        print("CERTIFICATE VALID AND COMPLETE: the selected set misses this "
              "physical census-strip face throughout the named teacher "
              "region.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
