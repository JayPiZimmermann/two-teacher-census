"""Exact parallel replay of an ordinary selected-minimum DFS forest.

The artifact and its semantics are unchanged.  This checker cuts each stored
root at a deterministic tree depth (or earlier when the stored tree already
has a leaf), assigns the resulting prefix-free cover to fresh worker
processes, and recomposes their results.

The coordinator strictly decodes and parses the ENTIRE bitstream, rejects
truncation, extra bytes/bits, nonzero padding, malformed prefix structure,
unknown verdicts, UNDECIDED leaves, incompatible specs, and stored statistics
that disagree with the tree.  Every worker reopens the same SHA-256-identified
artifact, reconstructs the same prefix map, regenerates its boxes via
``minimum_dfs.child``, and uses the same
``minimum_dfs.verdict(spec, *box, target=claimed)`` proof-method check as the
serial checker.  The coordinator requires every prefix exactly once and
rechecks total nodes, verdict counts, and geometric coverage.

The bounded exact atom cache is installed independently in each worker.  Use
``--no-cache`` to run the identical parallel decomposition without it.

Usage::

  python3 minimum_dfs_check_parallel.py ARTIFACT --depth 8 --workers 16

When called by ``minimum_atlas_check.py``, set ``MINIMUM_REPLAY_WORKERS`` to
control the nested worker count because the atlas invokes checker scripts with
only the artifact argument.
"""
import argparse
import base64
import concurrent.futures
import hashlib
import json
import multiprocessing
import os

import minimum_dfs as D
import minimum_dfs_cached as C


HERE = os.path.dirname(os.path.abspath(__file__))


class Reader:
    def __init__(self, bits, start=0, end=None):
        self.bits = bits
        self.pos = start
        self.end = len(bits) if end is None else end

    def read(self, width=1):
        if self.pos + width > self.end:
            raise ValueError("bitstream is truncated")
        value = 0
        for _ in range(width):
            value = (value << 1) | self.bits[self.pos]
            self.pos += 1
        return value


class BitString:
    """A checked, compact MSB-first view of the encoded tree bits."""

    def __init__(self, raw, n_bits):
        self.raw = raw
        self.n_bits = n_bits

    def __len__(self):
        return self.n_bits

    def __getitem__(self, index):
        if not 0 <= index < self.n_bits:
            raise IndexError(index)
        return (self.raw[index >> 3] >> (7 - (index & 7))) & 1


def decode_bits(encoded, n_bits):
    raw = base64.b64decode(encoded, validate=True)
    if type(n_bits) is not int or n_bits < 0:
        raise ValueError("invalid bit count")
    expected = (n_bits + 7) // 8
    if len(raw) != expected:
        raise ValueError("bitstream has truncated or extra encoded bytes")
    padding = (-n_bits) % 8
    if raw and padding and raw[-1] & ((1 << padding) - 1):
        raise ValueError("bitstream has nonzero padding bits")
    return BitString(raw, n_bits)


def validate_document(document):
    try:
        spec = document["spec"]
        lo_idx, hi_idx = document["lo_idx"], document["hi_idx"]
    except (KeyError, TypeError) as exc:
        raise ValueError("document is missing its forest specification") \
            from exc
    if not (isinstance(spec, dict) and D.partition_covers(spec)
            and type(lo_idx) is int and type(hi_idx) is int
            and 0 <= lo_idx < hi_idx <= spec["n_parts"]
            and 0 < spec["delta"] < spec["dmax"] < D.mpmath.pi):
        raise ValueError("inadmissible spec or partition")
    return spec, lo_idx, hi_idx


def empty_counts():
    return {name: 0 for name in D.VERDICTS}


def add_counts(target, source):
    for name in D.VERDICTS:
        target[name] += source[name]


def unit_key(cell, path):
    return "%d:%s" % (cell, path)


def scan_node(reader, absolute_depth):
    """Parse one complete subtree without doing interval arithmetic."""
    if absolute_depth > 500:
        raise ValueError("tree depth exceeds 500")
    nodes = 1
    counts = empty_counts()
    if reader.read() == 1:
        left_nodes, left_counts = scan_node(reader, absolute_depth + 1)
        right_nodes, right_counts = scan_node(reader, absolute_depth + 1)
        nodes += left_nodes + right_nodes
        add_counts(counts, left_counts)
        add_counts(counts, right_counts)
        return nodes, counts
    code = reader.read(3)
    name = D.NAME.get(code)
    if name is None:
        raise ValueError("unknown verdict code")
    counts[name] += 1
    return nodes, counts


def scan_document(document, frontier_depth):
    """Return the exact prefix-free replay units and structural totals."""
    if type(frontier_depth) is not int or not 0 <= frontier_depth <= 32:
        raise ValueError("frontier depth must lie in [0,32]")
    spec, lo_idx, hi_idx = validate_document(document)
    bits = decode_bits(document.get("bits_b64", ""),
                       document.get("n_nodes"))
    reader = Reader(bits)
    units = {}
    prefix_nodes = 0

    def cut(cell, path, level):
        nonlocal prefix_nodes
        start = reader.pos
        marker = reader.read()
        # A leaf before the requested depth is itself one member of the
        # prefix-free cover.  Otherwise every internal prefix node contributes
        # exactly its two deterministic children.
        if marker == 0 or level == frontier_depth:
            reader.pos = start
            nodes, counts = scan_node(reader, level)
            end = reader.pos
            key = unit_key(cell, path)
            if key in units:
                raise ValueError("duplicate replay unit")
            units[key] = {"cell": cell, "path": path,
                          "start": start, "end": end,
                          "nodes": nodes, "counts": counts}
            return
        prefix_nodes += 1
        cut(cell, path + "0", level + 1)
        cut(cell, path + "1", level + 1)

    for cell in range(lo_idx, hi_idx):
        cut(cell, "", 0)
    if reader.pos != len(bits):
        raise ValueError("bitstream has extra tree bits")

    structural_counts = empty_counts()
    structural_nodes = prefix_nodes
    for unit in units.values():
        structural_nodes += unit["nodes"]
        add_counts(structural_counts, unit["counts"])
    stats = document.get("stats")
    stored_counts = stats.get("counts") if isinstance(stats, dict) else None
    counts_are_integers = isinstance(stored_counts, dict) \
        and set(stored_counts) == set(D.VERDICTS) \
        and all(type(value) is int and value >= 0
                for value in stored_counts.values())
    if not isinstance(stats, dict) \
            or type(stats.get("visited_nodes")) is not int \
            or not counts_are_integers \
            or stats.get("visited_nodes") != structural_nodes \
            or stored_counts != structural_counts \
            or type(stats.get("unresolved_volume")) not in (int, float) \
            or stats.get("unresolved_volume") != 0.0:
        raise ValueError("stored statistics do not recompose from the tree")
    if structural_counts["UNDECIDED"]:
        raise ValueError("certificate contains UNDECIDED leaves")

    roots = D.root_cells(spec)
    prefix_volume = 0.0
    for unit in units.values():
        box = roots[unit["cell"]]
        for bit in unit["path"]:
            box = D.child(spec, box, int(bit))
        unit["box"] = box
        prefix_volume += D.volume(box)
    root_volume = sum(D.volume(roots[cell])
                      for cell in range(lo_idx, hi_idx))
    tolerance = 2e-12 * max(1.0, abs(root_volume))
    if abs(prefix_volume - root_volume) > tolerance:
        raise ValueError("prefix units do not cover the stored roots")
    return {"spec": spec, "lo_idx": lo_idx, "hi_idx": hi_idx,
            "bits": bits, "units": units,
            "prefix_nodes": prefix_nodes,
            "structural_nodes": structural_nodes,
            "structural_counts": structural_counts,
            "root_volume": root_volume}


def load_hashed(path, expected_hash=None):
    with open(path, "rb") as handle:
        raw = handle.read()
    digest = hashlib.sha256(raw).hexdigest()
    if expected_hash is not None and digest != expected_hash:
        raise ValueError("artifact changed during parallel replay")
    return json.loads(raw), digest


def replay_unit(scanned, unit):
    spec, bits = scanned["spec"], scanned["bits"]
    reader = Reader(bits, unit["start"], unit["end"])
    counts = empty_counts()
    nodes = 0
    leaf_volume = 0.0

    def walk(box, absolute_depth):
        nonlocal nodes, leaf_volume
        if absolute_depth > 500:
            raise ValueError("tree depth exceeds 500")
        nodes += 1
        if reader.read() == 1:
            walk(D.child(spec, box, 0), absolute_depth + 1)
            walk(D.child(spec, box, 1), absolute_depth + 1)
            return
        code = reader.read(3)
        claimed = D.NAME.get(code)
        if claimed is None:
            raise ValueError("unknown verdict code")
        if claimed == "UNDECIDED":
            raise ValueError("certificate contains an UNDECIDED leaf")
        got = D.verdict(spec, *box, target=claimed)
        if got != claimed:
            raise ValueError("leaf method mismatch: claimed %s, got %s" %
                             (claimed, got))
        counts[claimed] += 1
        leaf_volume += D.volume(box)

    walk(unit["box"], len(unit["path"]))
    if reader.pos != unit["end"]:
        raise ValueError("unit replay did not consume its exact bit slice")
    if nodes != unit["nodes"] or counts != unit["counts"]:
        raise ValueError("unit numerical replay disagrees with structural scan")
    expected_volume = D.volume(unit["box"])
    tolerance = 2e-12 * max(1.0, abs(expected_volume))
    if abs(leaf_volume - expected_volume) > tolerance:
        raise ValueError("unit leaves do not cover their prefix box")
    return nodes, counts, leaf_volume


def replay_group(path, digest, frontier_depth, keys, use_cache):
    """Fresh-process worker: independently reload, rescan, and replay keys."""
    document, actual = load_hashed(path, digest)
    if actual != digest:
        raise ValueError("worker artifact hash mismatch")
    scanned = scan_document(document, frontier_depth)
    if use_cache:
        C.install_atom_cache()
    units = scanned["units"]
    if len(set(keys)) != len(keys) or any(key not in units for key in keys):
        raise ValueError("worker received duplicate or unknown replay units")
    counts = empty_counts()
    nodes = 0
    leaf_volume = 0.0
    for key in keys:
        part_nodes, part_counts, part_volume = replay_unit(
            scanned, units[key])
        nodes += part_nodes
        add_counts(counts, part_counts)
        leaf_volume += part_volume
    return {"keys": keys, "nodes": nodes, "counts": counts,
            "leaf_volume": leaf_volume, "pid": os.getpid()}


def balanced_groups(units, workers):
    """Assign every unit once, greedily balancing structural node counts."""
    if workers < 1:
        raise ValueError("workers must be positive")
    count = min(workers, len(units))
    if count == 0:
        raise ValueError("forest has no replay units")
    groups = [[] for _ in range(count)]
    loads = [0] * count
    for key, unit in sorted(units.items(),
                            key=lambda item: (-item[1]["nodes"], item[0])):
        index = min(range(count), key=lambda i: (loads[i], i))
        groups[index].append(key)
        loads[index] += unit["nodes"]
    return groups


def validate_assignment(groups, expected):
    flat = [key for group in groups for key in group]
    if len(flat) != len(set(flat)):
        raise ValueError("replay assignment duplicates a prefix")
    if set(flat) != set(expected):
        raise ValueError("replay assignment misses or invents a prefix")


def parallel_replay(path, frontier_depth, workers, use_cache=True):
    document, digest = load_hashed(path)
    scanned = scan_document(document, frontier_depth)
    groups = balanced_groups(scanned["units"], workers)
    validate_assignment(groups, scanned["units"])
    results = []
    context = multiprocessing.get_context("spawn")
    with concurrent.futures.ProcessPoolExecutor(
            max_workers=len(groups), mp_context=context) as pool:
        pending = [pool.submit(replay_group, path, digest, frontier_depth,
                               group, use_cache)
                   for group in groups]
        for future in concurrent.futures.as_completed(pending):
            results.append(future.result())
    validate_assignment([result["keys"] for result in results],
                        scanned["units"])
    counts = empty_counts()
    subtree_nodes = 0
    leaf_volume = 0.0
    for result in results:
        subtree_nodes += result["nodes"]
        leaf_volume += result["leaf_volume"]
        add_counts(counts, result["counts"])
    total_nodes = scanned["prefix_nodes"] + subtree_nodes
    tolerance = 2e-12 * max(1.0, abs(scanned["root_volume"]))
    if total_nodes != scanned["structural_nodes"] \
            or counts != scanned["structural_counts"] \
            or abs(leaf_volume - scanned["root_volume"]) > tolerance:
        raise ValueError("worker results do not recompose the ordinary forest")
    _, final_digest = load_hashed(path, digest)
    if final_digest != digest:
        raise ValueError("artifact changed after parallel replay")
    return {"document": document, "digest": digest, "scanned": scanned,
            "groups": groups, "results": results,
            "counts": counts, "nodes": total_nodes,
            "leaf_volume": leaf_volume}


def default_workers():
    raw = os.environ.get("MINIMUM_REPLAY_WORKERS", "8")
    try:
        return int(raw)
    except ValueError as exc:
        raise ValueError("MINIMUM_REPLAY_WORKERS is not an integer") from exc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact")
    parser.add_argument("--depth", type=int, default=8)
    try:
        worker_default = default_workers()
    except ValueError as exc:
        parser.error(str(exc))
    parser.add_argument("--workers", type=int, default=worker_default)
    parser.add_argument("--no-cache", action="store_true")
    args = parser.parse_args()
    if args.workers < 1 or not 0 <= args.depth <= 32:
        parser.error("workers must be positive and depth must lie in [0,32]")
    path = args.artifact if os.path.isabs(args.artifact) \
        else os.path.join(HERE, args.artifact)
    try:
        result = parallel_replay(
            path, args.depth, args.workers, not args.no_cache)
    except Exception as exc:
        print("PARALLEL REPLAY FAILED:", exc)
        return 1
    scanned = result["scanned"]
    print("certificate :", os.path.basename(path))
    print("cells       : [%d,%d) of %d" %
          (scanned["lo_idx"], scanned["hi_idx"],
           scanned["spec"]["n_parts"]))
    print("sha256      :", result["digest"])
    print("bits        :", result["document"]["n_nodes"])
    print("units       : %d at depth <= %d over %d fresh workers" %
          (len(scanned["units"]), args.depth, len(result["results"])))
    print("leaves      :", json.dumps(result["counts"], sort_keys=True))
    print("coverage    : %.16e leaf volume / %.16e root volume" %
          (result["leaf_volume"], scanned["root_volume"]))
    scope = scanned["spec"].get("scope", "selected")
    if scope == "census":
        claim = "EVERY angle-map zero is fold-free"
    else:
        claim = "every selected zero is regular"
    print("PARALLEL CERTIFICATE VALID AND COMPLETE: %s; every prefix and "
          "leaf recomputed%s." %
          (claim, " with exact atom caching" if not args.no_cache else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
