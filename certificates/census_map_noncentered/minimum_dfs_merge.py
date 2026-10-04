"""Merge selected-minimum DFS shards into one replayable certificate.

Every input must have the same exact spec, contain no undecided leaf, and the
input ranges must tile ``[0, n_parts)`` without a gap or overlap.  Bitstreams
are concatenated in increasing root-cell order, which is exactly the order
``minimum_dfs_check.py`` walks.  The merged artifact is still replayed by that
checker; this script adds no numerical trust.

Usage: python3 minimum_dfs_merge.py <tag>
       python3 minimum_dfs_merge.py <tag> --output <path>
"""
import argparse
import base64
import glob
import json
import os

import minimum_dfs as D


HERE = os.path.dirname(os.path.abspath(__file__))


def decode_bits(encoded, n_bits):
    raw = base64.b64decode(encoded, validate=True)
    if n_bits < 0 or n_bits > 8 * len(raw):
        raise ValueError("invalid bit count %d for %d encoded bytes" %
                         (n_bits, len(raw)))
    return [((raw[index >> 3] >> (7 - (index & 7))) & 1)
            for index in range(n_bits)]


def encode_bits(bits):
    padding = (-len(bits)) % 8
    packed = bytearray()
    for start in range(0, len(bits) + padding, 8):
        byte = 0
        for offset in range(8):
            index = start + offset
            byte = (byte << 1) | (bits[index] if index < len(bits) else 0)
        packed.append(byte)
    return base64.b64encode(bytes(packed)).decode()


def merge_docs(docs):
    if not docs:
        raise ValueError("no shard documents")
    spec = docs[0].get("spec")
    if not isinstance(spec, dict) or not D.partition_covers(spec):
        raise ValueError("first shard has an inadmissible spec")
    for doc in docs[1:]:
        if doc.get("spec") != spec:
            raise ValueError("shards disagree on the exact spec")
    ordered = sorted(docs, key=lambda doc: doc.get("lo_idx", -1))
    next_cell = 0
    bits = []
    counts = {name: 0 for name in D.VERDICTS}
    visited_nodes = 0
    unresolved_volume = 0.0
    seconds = 0.0
    ranges = []
    for doc in ordered:
        lo_idx = doc.get("lo_idx")
        hi_idx = doc.get("hi_idx")
        if not isinstance(lo_idx, int) or not isinstance(hi_idx, int) \
                or lo_idx != next_cell or not lo_idx < hi_idx \
                or hi_idx > spec["n_parts"]:
            raise ValueError("partition gap/overlap at expected cell %d: %r" %
                             (next_cell, (lo_idx, hi_idx)))
        shard_stats = doc.get("stats")
        shard_counts = shard_stats.get("counts", {}) \
            if isinstance(shard_stats, dict) else {}
        if set(shard_counts) != set(D.VERDICTS):
            raise ValueError("shard [%d,%d) has incomplete verdict stats" %
                             (lo_idx, hi_idx))
        if shard_counts.get("UNDECIDED"):
            raise ValueError("shard [%d,%d) contains undecided leaves" %
                             (lo_idx, hi_idx))
        bits.extend(decode_bits(doc.get("bits_b64", ""),
                                doc.get("n_nodes", -1)))
        for name in D.VERDICTS:
            counts[name] += shard_counts[name]
        visited_nodes += shard_stats.get("visited_nodes", 0)
        unresolved_volume += shard_stats.get("unresolved_volume", 0.0)
        seconds += shard_stats.get("seconds", 0.0)
        ranges.append([lo_idx, hi_idx])
        next_cell = hi_idx
    if next_cell != spec["n_parts"]:
        raise ValueError("partition incomplete: covered [0,%d) of %d cells" %
                         (next_cell, spec["n_parts"]))
    return {"spec": spec, "lo_idx": 0, "hi_idx": spec["n_parts"],
            "n_nodes": len(bits), "bits_b64": encode_bits(bits),
            "stats": {"visited_nodes": visited_nodes, "counts": counts,
                      "unresolved_volume": unresolved_volume,
                      "seconds": seconds, "merged_ranges": ranges}}


def load_docs(paths):
    docs = []
    for path in paths:
        with open(path) as handle:
            docs.append(json.load(handle))
    return docs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tag")
    parser.add_argument("--output")
    args = parser.parse_args()
    pattern = os.path.join(HERE, "minimum_dfs_%s_p*.json" % args.tag)
    paths = sorted(glob.glob(pattern))
    if not paths:
        print("MERGE FAILED: no shards found for tag %s" % args.tag)
        return 1
    try:
        merged = merge_docs(load_docs(paths))
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print("MERGE FAILED:", exc)
        return 1
    output = args.output or os.path.join(
        HERE, "minimum_dfs_%s_merged.json" % args.tag)
    if not os.path.isabs(output):
        output = os.path.join(HERE, output)
    temporary = output + ".tmp"
    with open(temporary, "w") as handle:
        json.dump(merged, handle, separators=(",", ":"))
    os.replace(temporary, output)
    print("merged %d shards covering all %d cells -> %s (%d bits)" %
          (len(paths), merged["spec"]["n_parts"],
           os.path.basename(output), merged["n_nodes"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
