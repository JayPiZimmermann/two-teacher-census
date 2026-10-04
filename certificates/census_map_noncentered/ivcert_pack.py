"""Pack a separated-enumeration certificate from path-keyed JSON to a bitstring.

The certificate is a bisection tree (see `ivcert_check.py`).  Storing it as a
dictionary of path strings is ~40 bytes per leaf; the tree needs FIVE BITS.
This converts one to the other and verifies the round trip, so the packed file
carries exactly the same object.

THE ENCODING, in full — a referee must be able to write the reader from this
paragraph alone, which is why it is also copied into README.md.

* The tree is walked in PREORDER, root cells in index order, and within a node
  the child `0` (lower half) before the child `1` (upper half).
* At every node exactly one bit is written first: **`1` = internal**, its two
  children follow immediately in the stream; **`0` = leaf**, and then three
  more bits give the verdict.
* The split is never encoded because it is determined by the box: bisect the
  WIDER side (`s`-side when `s1 - s0 >= d1 - d0`), child `0` is the lower
  half.  A reader that regenerates the root cells from the spec therefore
  reconstructs every box exactly, with no coordinate in the file.
* Verdict codes, three bits:

      0  PREFILTER          the sign law excludes the box
      1  EXCL_PLAIN         plain interval enclosure of G0 or G1 misses zero
      2  EXCL_CENTERED      mean-value enclosure of G0 or G1 misses zero
      3  KRAWCZYK_EMPTY     K(X) cap X empty
      4  KRAWCZYK_UNIQUE    K(X) subset int(X): exactly one solution
      5  KRAWCZYK_INFLATED  the same on an inflated box
      6  UNDECIDED_MINWIDTH no verdict, minimum width reached
      7  UNDECIDED_OTHER    no verdict, any other reason

  Codes 6 and 7 are not distinguished further because the checker treats ANY
  undecided leaf the same way: it invalidates the completeness claim and the
  replay exits non-zero.  The unpacked reason is in `signlaw_bb_<tag>.json`,
  which is prospecting, not evidence.
* Bits are packed MSB-first into bytes, the final byte zero-padded, and the
  byte string base64-encoded as `bits_b64`.  `n_nodes` records how many bits
  of the last byte are meaningful, so the padding can never be misread as
  extra tree.

Usage:  python3 ivcert_pack.py ivcert_<tag>.json    -> ivcert_<tag>.bits.json
"""
import base64
import json
import os
import sys

VERDICTS = ["PREFILTER", "EXCL_PLAIN", "EXCL_CENTERED", "KRAWCZYK_EMPTY",
            "KRAWCZYK_UNIQUE", "KRAWCZYK_INFLATED", "UNDECIDED_MINWIDTH",
            "UNDECIDED_OTHER"]
CODE = {v: i for i, v in enumerate(VERDICTS)}
HERE = os.path.dirname(os.path.abspath(__file__))


def _code_of(tag):
    if tag in CODE:
        return CODE[tag]
    return CODE["UNDECIDED_OTHER"] if tag.startswith("UNDECIDED") else None


class BitWriter(object):
    def __init__(self):
        self.bits = []

    def w(self, value, width):
        for i in range(width - 1, -1, -1):
            self.bits.append((value >> i) & 1)

    def to_bytes(self):
        out = bytearray()
        for i in range(0, len(self.bits), 8):
            byte = 0
            for j in range(8):
                byte = (byte << 1) | (self.bits[i + j]
                                      if i + j < len(self.bits) else 0)
            out.append(byte)
        return bytes(out)


class BitReader(object):
    def __init__(self, data, n_bits):
        self.data = data
        self.n = n_bits
        self.i = 0

    def r(self, width=1):
        v = 0
        for _ in range(width):
            if self.i >= self.n:
                raise ValueError("bitstream exhausted")
            byte = self.data[self.i >> 3]
            v = (v << 1) | ((byte >> (7 - (self.i & 7))) & 1)
            self.i += 1
        return v


def pack(doc):
    """Path-keyed leaves -> preorder bitstream."""
    by_root = {}
    for key, v in doc["leaves"].items():
        ri, p = key.split(":", 1)
        by_root.setdefault(int(ri), {})[p] = v
    bw = BitWriter()

    def walk(ri, p, depth):
        v = by_root.get(ri, {}).get(p)
        if v is None:
            if depth > 400:
                raise ValueError("no leaf at root %d path %s" % (ri, p))
            bw.w(1, 1)
            walk(ri, p + "0", depth + 1)
            walk(ri, p + "1", depth + 1)
            return
        code = _code_of(v)
        if code is None:
            raise ValueError("unknown verdict %r" % v)
        bw.w(0, 1)
        bw.w(code, 3)

    sys.setrecursionlimit(20000)
    for ri in range(doc["spec"]["n_root_cells"]):
        walk(ri, "", 0)
    return bw


def unpack(bits_b64, n_bits, n_root_cells):
    """Preorder bitstream -> path-keyed leaves (the exact inverse)."""
    br = BitReader(base64.b64decode(bits_b64), n_bits)
    leaves = {}

    def walk(ri, p):
        if br.r(1) == 1:
            walk(ri, p + "0")
            walk(ri, p + "1")
            return
        leaves["%d:%s" % (ri, p)] = VERDICTS[br.r(3)]

    sys.setrecursionlimit(20000)
    for ri in range(n_root_cells):
        walk(ri, "")
    return leaves


def main():
    src = sys.argv[1]
    path = src if os.path.isabs(src) else os.path.join(HERE, src)
    doc = json.load(open(path))
    bw = pack(doc)
    raw = bw.to_bytes()
    b64 = base64.b64encode(raw).decode()

    # ROUND TRIP, always: the packed file must decode to the SAME object, or
    # it is not the same certificate.
    back = unpack(b64, len(bw.bits), doc["spec"]["n_root_cells"])
    orig = {k: (v if v in CODE else "UNDECIDED_OTHER")
            for k, v in doc["leaves"].items()}
    if back != orig:
        print("ROUND TRIP FAILED: %d vs %d leaves" % (len(back), len(orig)))
        sys.exit(1)

    out = {"spec": doc["spec"],
           "format": "preorder bisection tree; 1 bit per node "
                     "(1=internal, 0=leaf), leaves followed by a 3-bit "
                     "verdict code; split = bisect the wider side, child 0 "
                     "= lower half; bits packed MSB-first, base64.  Verdict "
                     "codes: " + ", ".join("%d=%s" % (i, v)
                                           for i, v in enumerate(VERDICTS)),
           "n_nodes": len(bw.bits),
           "bits_b64": b64}
    dst = path.replace(".json", ".bits.json")
    with open(dst, "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    print("%s -> %s" % (os.path.basename(path), os.path.basename(dst)))
    print("  leaves %d, nodes %d, %.1f KB -> %.1f KB (%.0fx smaller), "
          "round trip OK"
          % (len(orig), len(bw.bits), os.path.getsize(path) / 1024.0,
             os.path.getsize(dst) / 1024.0,
             os.path.getsize(path) / max(1.0, float(os.path.getsize(dst)))))


if __name__ == "__main__":
    main()
