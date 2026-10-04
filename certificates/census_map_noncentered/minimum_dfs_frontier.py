"""Checker-compatible DFS below one deterministic within-cell frontier.

The ordinary selected-minimum forest has two independent notions of
partitioning:

* ``n_parts`` partitions the full-period ``u`` coordinate into root cells;
* the deterministic binary tree below each root partitions all four
  coordinates.

``minimum_shard_escalate.py`` gives each root cell its own node budget.  A
pathological root can still exhaust that budget even though all of its
children are individually tractable.  This module gives a deterministic
subtree (named by a binary path below one root cell) its own budget.  A full
fixed-depth frontier can then be merged back into one completely ordinary
``minimum_dfs`` shard.  The unchanged ``minimum_dfs_check.py`` replays that
shard from its root and is the final authority.

The optional map-first search policy changes neither the leaf predicates nor
the certificate format.  Above ``full_width`` it tries only angle-map
exclusions, postponing Krawczyk/selector classification until the box is
finer.  This is useful in the high-angle p61 region, where every observed
leaf was a map exclusion but the failed parents repeatedly paid for selector
work.  A postponed test can only cause another deterministic split; it can
never turn an unproved box into a leaf.

Worker usage::

  python3 minimum_dfs_frontier.py TEMPLATE CELL PATH BUDGET OUTPUT \
      [--full-width 0.001] [--no-cache]

``PATH`` is a string of ``0`` and ``1`` (use ``-`` for the empty path).
``TEMPLATE`` is any ordinary shard carrying the exact target spec; its bits
are not consumed.
"""
import argparse
import base64
import json
import os
import time

import minimum_dfs as D
import minimum_dfs_cached as C


FORMAT = "minimum-dfs-frontier-v1"


def decode_bits(encoded, n_bits):
    raw = base64.b64decode(encoded, validate=True)
    if not isinstance(n_bits, int) or n_bits < 0 or n_bits > 8 * len(raw):
        raise ValueError("invalid bit count")
    expected_bytes = (n_bits + 7) // 8
    if len(raw) != expected_bytes:
        raise ValueError("bitstream has extra encoded bytes")
    padding = (-n_bits) % 8
    if raw and padding and raw[-1] & ((1 << padding) - 1):
        raise ValueError("bitstream has nonzero padding bits")
    return [((raw[index >> 3] >> (7 - (index & 7))) & 1)
            for index in range(n_bits)]


def box_at_path(spec, cell, path):
    """Regenerate a subtree root using only the stored spec and its path."""
    if not D.partition_covers(spec) or not 0 <= cell < spec["n_parts"]:
        raise ValueError("inadmissible spec or root cell")
    if any(bit not in "01" for bit in path):
        raise ValueError("frontier path is not binary")
    box = D.root_cells(spec)[cell]
    for bit in path:
        box = D.child(spec, box, int(bit))
    return box


def map_first_verdict(spec, *box, full_width=0.001):
    """Use only existing proof methods, delaying costly selector tests.

    Once a box is at most ``full_width`` in every coordinate, delegate to the
    ordinary verdict function.  Hence this policy retains all of the original
    termination mechanisms near a genuine selected root.  On larger boxes it
    evaluates the four map exclusions in an order suited to the pathological
    high-angle cells.  Replay does not trust this ordering: it re-derives the
    named leaf method directly with :func:`minimum_dfs.verdict`.
    """
    widths = tuple(box[2 * i + 1] - box[2 * i] for i in range(4))
    if max(widths) <= full_width:
        return D.verdict(spec, *box)
    bl, bh, yl, yh, sl, sh, dl, dh = box
    B, Y = D.mk(bl, bh), D.mk(yl, yh)
    S, gap = D.mk(sl, sh), D.mk(dl, dh)
    try:
        # In the observed p61 trees this is overwhelmingly the successful
        # leaf method.  Trying it first avoids work whose result would not be
        # stored in the certificate.
        if max(widths) <= spec.get("row_mv_maxw", 0.1):
            centered = D.rowform_mean_value(spec, *box)
            if centered is not None and (
                    D.X.sgn(centered[0]) or D.X.sgn(centered[1])):
                return "MAP_CENTERED"
        row0, row1 = D.coordinate_rowform(spec, B, Y, S, gap)
        if D.X.sgn(row0) or D.X.sgn(row1):
            return "MAP_PLAIN"
        num0, num1 = D.coordinate_mass_numerators(spec, B, Y, S, gap)
        if D.X.sgn(num0) and D.X.sgn(num1) \
                and D.signlaw_kills(spec, B, Y, S, gap):
            return "MAP_AUX"
        collar = D.collar_rowform(spec, B, Y, S, gap)
        if collar is not None and (
                D.X.sgn(collar[0]) or D.X.sgn(collar[1])):
            return "MAP_COLLAR"
    except Exception:
        # The ordinary generator also treats an interval-method exception as
        # failure of that method and subdivides the box.
        pass
    return None


def run_box(spec, box, budget, full_width=0.001):
    """Run the ordinary preorder DFS below one already-regenerated box."""
    bits = []
    counts = {name: 0 for name in D.VERDICTS}
    nodes = 0
    unresolved_volume = 0.0
    stack = [box]
    started = time.time()
    while stack:
        current = stack.pop()
        nodes += 1
        if nodes > budget:
            raise RuntimeError("budget exhausted in frontier subtree")
        result = map_first_verdict(
            spec, *current, full_width=full_width)
        if result is None:
            width = max(current[2 * i + 1] - current[2 * i]
                        for i in range(4))
            if width < spec["minw"]:
                result = "UNDECIDED"
                unresolved_volume += D.volume(current)
            else:
                bits.append(1)
                stack.append(D.child(spec, current, 1))
                stack.append(D.child(spec, current, 0))
                continue
        bits.append(0)
        code = D.CODE[result]
        bits.extend(((code >> 2) & 1, (code >> 1) & 1, code & 1))
        counts[result] += 1
    return {
        "n_nodes": len(bits),
        "bits_b64": D.encode(bits),
        "stats": {
            "visited_nodes": nodes,
            "counts": counts,
            "unresolved_volume": unresolved_volume,
            "seconds": time.time() - started,
        },
    }


def make_fragment(spec, cell, path, budget, full_width=0.001):
    subtree = run_box(spec, box_at_path(spec, cell, path), budget,
                      full_width=full_width)
    return {
        "format": FORMAT,
        "spec": spec,
        "cell_idx": cell,
        "path": path,
        "n_nodes": subtree["n_nodes"],
        "bits_b64": subtree["bits_b64"],
        "stats": subtree["stats"],
        "search": {"policy": "map-first-then-full",
                   "full_width": full_width,
                   "exact_atom_cache": C.D.J.atoms is C.cached_atoms},
    }


def audit_fragment(spec, cell, path, fragment):
    """Check one fragment's exact tree, statistics, and geometric coverage."""
    bits = decode_bits(fragment.get("bits_b64", ""),
                       fragment.get("n_nodes"))
    position = 0
    nodes = 0
    counts = {name: 0 for name in D.VERDICTS}
    leaf_volume = 0.0

    def read(width=1):
        nonlocal position
        if position + width > len(bits):
            raise ValueError("frontier bitstream is truncated")
        value = 0
        for _ in range(width):
            value = (value << 1) | bits[position]
            position += 1
        return value

    def walk(box, depth=0):
        nonlocal nodes, leaf_volume
        if depth > 500:
            raise ValueError("frontier tree depth exceeds 500")
        nodes += 1
        if read() == 1:
            walk(D.child(spec, box, 0), depth + 1)
            walk(D.child(spec, box, 1), depth + 1)
            return
        code = read(3)
        name = D.NAME.get(code)
        if name is None:
            raise ValueError("frontier leaf has an unknown verdict")
        counts[name] += 1
        leaf_volume += D.volume(box)

    root = box_at_path(spec, cell, path)
    walk(root)
    if position != len(bits):
        raise ValueError("frontier bitstream has extra tree bits")
    stats = fragment.get("stats", {})
    if stats.get("counts") != counts:
        raise ValueError("frontier verdict statistics do not match its tree")
    if stats.get("visited_nodes") != nodes:
        raise ValueError("frontier node statistics do not match its tree")
    if counts["UNDECIDED"] or stats.get("unresolved_volume") != 0.0:
        raise ValueError("frontier fragment contains undecided leaves")
    expected_volume = D.volume(root)
    tolerance = 2e-12 * max(1.0, abs(expected_volume))
    if abs(leaf_volume - expected_volume) > tolerance:
        raise ValueError("frontier fragment does not cover its subtree")
    return bits, counts, nodes, leaf_volume


def merge_fragments(spec, cell, depth, fragments):
    """Join a complete, possibly refined frontier into an ordinary shard.

    ``depth`` is the initial frontier depth.  Every fragment path must be at
    least that long, but any one initial path may be replaced recursively by
    its two children (and either child may be refined again).  Thus the path
    set must be a prefix-free complete binary cover, not necessarily a
    uniform-depth set.
    """
    if depth < 0 or not D.partition_covers(spec) \
            or not 0 <= cell < spec["n_parts"]:
        raise ValueError("inadmissible merge target")
    by_path = {}
    counts = {name: 0 for name in D.VERDICTS}
    subtree_nodes = 0
    unresolved_volume = 0.0
    leaf_volume = 0.0
    seconds = 0.0
    for fragment in fragments:
        if fragment.get("format") != FORMAT \
                or fragment.get("spec") != spec \
                or fragment.get("cell_idx") != cell:
            raise ValueError("frontier fragment metadata mismatch")
        path = fragment.get("path")
        if not isinstance(path, str) or any(bit not in "01" for bit in path):
            raise ValueError("frontier fragment has an invalid path")
        if len(path) < depth:
            raise ValueError("frontier path is shallower than its base depth")
        if path in by_path:
            raise ValueError("duplicate frontier path %r" % path)
        stats = fragment.get("stats", {})
        fragment_bits, part_counts, part_nodes, part_volume = audit_fragment(
            spec, cell, path, fragment)
        by_path[path] = fragment_bits
        for name in D.VERDICTS:
            counts[name] += part_counts[name]
        subtree_nodes += part_nodes
        unresolved_volume += stats["unresolved_volume"]
        leaf_volume += part_volume
        seconds += stats.get("seconds", 0.0)
    paths = set(by_path)
    for path in paths:
        if any(path[:stop] in paths for stop in range(len(path))):
            raise ValueError("frontier paths overlap by prefix at %r" % path)
    internal_prefixes = {path[:stop]
                         for path in paths for stop in range(len(path))}
    prefix_nodes = 0

    def tree(prefix):
        nonlocal prefix_nodes
        if prefix in by_path:
            return by_path[prefix]
        if prefix not in internal_prefixes:
            raise ValueError("frontier is incomplete below path %r" % prefix)
        prefix_nodes += 1
        return [1] + tree(prefix + "0") + tree(prefix + "1")

    bits = tree("")
    root_volume = D.volume(D.root_cells(spec)[cell])
    tolerance = 2e-12 * max(1.0, abs(root_volume))
    if abs(leaf_volume - root_volume) > tolerance:
        raise ValueError("frontier paths do not cover the coarse root cell")
    visited_nodes = prefix_nodes + subtree_nodes
    return {
        "spec": spec,
        "lo_idx": cell,
        "hi_idx": cell + 1,
        "n_nodes": len(bits),
        "bits_b64": D.encode(bits),
        "stats": {
            "visited_nodes": visited_nodes,
            "counts": counts,
            "unresolved_volume": unresolved_volume,
            "seconds": seconds,
            "frontier_base_depth": depth,
            "frontier_max_depth": max(map(len, paths)),
            "frontier_fragments": len(paths),
        },
    }


def load_spec(path):
    with open(path) as handle:
        document = json.load(handle)
    spec = document.get("spec")
    if not isinstance(spec, dict) or not D.partition_covers(spec):
        raise ValueError("template has no admissible minimum-DFS spec")
    return spec


def write_json_atomic(path, document):
    temporary = path + ".tmp.%d" % os.getpid()
    with open(temporary, "w") as handle:
        json.dump(document, handle, separators=(",", ":"))
    os.replace(temporary, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("template")
    parser.add_argument("cell", type=int)
    parser.add_argument("path")
    parser.add_argument("budget", type=int)
    parser.add_argument("output")
    parser.add_argument("--full-width", type=float, default=0.001)
    parser.add_argument("--no-cache", action="store_true")
    args = parser.parse_args()
    path = "" if args.path == "-" else args.path
    if args.budget < 1 or args.full_width < 0:
        parser.error("budget must be positive and full-width nonnegative")
    try:
        spec = load_spec(args.template)
        if not args.no_cache:
            C.install_atom_cache()
        fragment = make_fragment(spec, args.cell, path, args.budget,
                                 full_width=args.full_width)
    except (OSError, ValueError, RuntimeError) as exc:
        print(exc)
        return 2
    write_json_atomic(args.output, fragment)
    print("wrote frontier cell %d path %s: %d visited nodes, %d bits" %
          (args.cell, path or "-", fragment["stats"]["visited_nodes"],
           fragment["n_nodes"]))
    if fragment["stats"]["counts"]["UNDECIDED"]:
        print("INCOMPLETE: surviving boxes are PSD-wall candidates")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
