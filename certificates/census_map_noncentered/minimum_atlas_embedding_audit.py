"""Read-only preflight for embedding recovered DFS roots in an atlas.

The recovery drivers publish ordinary ``n_parts=64`` shard files.  This tool
checks the exact path the atlas sweep will use without writing a merged
artifact:

* every matching canonical shard is strictly decoded and structurally
  recomposed with the parallel checker's scanner;
* filenames, scopes, boxes, partitions, statistics, and the full root-cell
  cover are checked exactly;
* :func:`minimum_atlas_sweep.usable_shards` and ``missing_ranges`` must see
  the same cover;
* :func:`minimum_dfs_merge.merge_docs` is run in memory and its result is
  strictly scanned; and
* any already-written merged forest must equal that in-memory result, while
  atlas discovery must select a structurally valid full forest.

Thus a depth-10 frontier result or a post-sweep adopted stage is accepted for
the same reason as any other ordinary shard: only its canonical ordinary JSON
schema matters.  No certificate, journal, or manifest is created or changed.

By default the CLI audits every non-journaled box in the requested target.
Use ``--tag`` to audit one or more exact tile tags while recovery is live.
"""
import argparse
import glob
import json
import os

import minimum_atlas_sweep as S
import minimum_dfs_check_parallel as P
import minimum_dfs_merge as M


class AuditError(ValueError):
    pass


def load_document(path):
    try:
        with open(path) as handle:
            document = json.load(handle)
    except (OSError, ValueError) as exc:
        raise AuditError("cannot load %s: %s" %
                         (os.path.basename(path), exc)) from exc
    if not isinstance(document, dict):
        raise AuditError("%s is not a JSON object" % os.path.basename(path))
    return document


def exact_missing_ranges(ranges, n_parts):
    """Return gaps, rejecting duplicate or overlapping claimed ranges."""
    gaps = []
    reach = 0
    for lo_idx, hi_idx in sorted(ranges):
        if lo_idx < reach:
            raise AuditError("duplicate/overlapping shard at [%d,%d)" %
                             (lo_idx, hi_idx))
        if lo_idx > reach:
            gaps.append((reach, lo_idx))
        reach = hi_idx
    if reach < n_parts:
        gaps.append((reach, n_parts))
    return gaps


def validate_expected_shard(path, document, tile, scope, n_parts):
    """Apply the sweep's exact tile test plus strict canonical naming."""
    try:
        scanned = P.scan_document(document, 0)
        spec = scanned["spec"]
        bounds = S.key_of_bounds(*(spec.get(name) for name in
                                   ("b0", "b1", "y0", "y1")))
    except Exception as exc:
        raise AuditError("strict structural scan failed for %s: %s" %
                         (os.path.basename(path), exc)) from exc
    if bounds != S.tile_key(tile) or spec.get("coord") != "relative" \
            or not S.same_strip(spec) or spec.get("n_parts") != n_parts \
            or S.spec_scope(spec) != scope:
        raise AuditError("tile/spec mismatch for %s" % os.path.basename(path))
    tag = S.tile_id(tile, scope)
    expected_name = ("minimum_dfs_%s_p%03d.json" %
                     (tag, scanned["lo_idx"]))
    if os.path.basename(path) != expected_name:
        raise AuditError("noncanonical shard name %s (expected %s)" %
                         (os.path.basename(path), expected_name))
    return scanned


def audit_recovered_forest(tile, scope="census", directory=None,
                           n_parts=64):
    """Strictly dry-merge one canonical forest without filesystem writes."""
    directory = S.HERE if directory is None else os.path.abspath(directory)
    tag = S.tile_id(tile, scope)
    pattern = os.path.join(directory, "minimum_dfs_%s_p*.json" % tag)
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise AuditError("%s has no canonical shards" % tag)
    documents = []
    ranges = []
    for path in paths:
        if os.path.islink(path) or not os.path.isfile(path):
            raise AuditError("shard is not a regular nonsymlink: %s" %
                             os.path.basename(path))
        document = load_document(path)
        scanned = validate_expected_shard(
            path, document, tile, scope, n_parts)
        documents.append(document)
        ranges.append((scanned["lo_idx"], scanned["hi_idx"]))
    gaps = exact_missing_ranges(ranges, n_parts)
    if gaps:
        raise AuditError("%s incomplete; missing %s" % (tag, gaps))

    # Compare directly to the live sweep's discovery surface when auditing its
    # own directory.  A malformed shard that usable_shards would skip is still
    # rejected above because minimum_dfs_merge would glob it later.
    if os.path.realpath(directory) == os.path.realpath(S.HERE):
        usable = S.usable_shards(tile, scope, n_parts, tag)
        if usable != sorted(ranges):
            raise AuditError("sweep usable_shards disagrees with strict cover")
        missing = S.missing_ranges(tile, scope, n_parts, tag)
        if missing:
            raise AuditError("sweep still sees missing ranges %s" % missing)

    try:
        merged = M.merge_docs(documents)
        merged_scan = P.scan_document(merged, 0)
    except Exception as exc:
        raise AuditError("exact dry merge failed for %s: %s" %
                         (tag, exc)) from exc
    if (merged_scan["lo_idx"], merged_scan["hi_idx"]) != (0, n_parts):
        raise AuditError("dry merge did not cover the full root partition")

    merged_path = os.path.join(
        directory, "minimum_dfs_%s_merged.json" % tag)
    if os.path.lexists(merged_path):
        if os.path.islink(merged_path) or not os.path.isfile(merged_path):
            raise AuditError("merged path is not a regular nonsymlink")
        existing = load_document(merged_path)
        try:
            P.scan_document(existing, 0)
        except Exception as exc:
            raise AuditError("existing merged forest is malformed: %s" % exc) \
                from exc
        if existing != merged:
            raise AuditError("existing merged forest is stale or differs "
                             "from canonical shards")

    discovered = None
    if os.path.realpath(directory) == os.path.realpath(S.HERE):
        regularity, _ = S.discover(scope)
        discovered = regularity.get(S.tile_key(tile))
        if discovered is not None:
            try:
                selected = load_document(discovered)
                selected_scan = P.scan_document(selected, 0)
            except Exception as exc:
                raise AuditError("atlas discovery selected an invalid full "
                                 "forest: %s" % exc) from exc
            if (selected_scan["lo_idx"], selected_scan["hi_idx"]) != \
                    (0, n_parts):
                raise AuditError("atlas discovery selected an incomplete forest")

    return {
        "tag": tag,
        "paths": paths,
        "ranges": sorted(ranges),
        "merged_nodes": merged_scan["structural_nodes"],
        "merged_counts": merged_scan["structural_counts"],
        "merged_present": os.path.exists(merged_path),
        "discovered": discovered,
    }


def target_tiles(y_min, beta_min):
    return S.selected_tiles(y_min=y_min, beta_min=beta_min)


def load_state(path):
    document = load_document(path)
    if document.get("format") != "minimum-atlas-sweep-v1" \
            or document.get("strip") != S.STRIP \
            or not isinstance(document.get("passed"), dict) \
            or not isinstance(document.get("failed"), dict):
        raise AuditError("state is not an exact minimum-atlas sweep journal")
    return document


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--state", default="minimum_atlas_f4_sweep_census_state.json")
    parser.add_argument("--scope", choices=("selected", "census"),
                        default="census")
    parser.add_argument("--y-min", default="-0.45")
    parser.add_argument("--beta-min", default="0.5")
    parser.add_argument("--tag", action="append", default=[],
                        help="audit this exact target tile tag; repeatable")
    args = parser.parse_args()
    state_path = args.state if os.path.isabs(args.state) else \
        os.path.join(S.HERE, args.state)
    try:
        state = load_state(state_path)
        tiles = target_tiles(args.y_min, args.beta_min)
        by_tag = {S.tile_id(tile, args.scope): tile for tile in tiles}
        if args.tag:
            unknown = sorted(set(args.tag) - set(by_tag))
            if unknown:
                raise AuditError("tags outside the exact target: %s" % unknown)
            chosen = [by_tag[tag] for tag in args.tag]
        else:
            chosen = [tile for tile in tiles
                      if S.state_key(tile) not in state["passed"]]
        if not chosen:
            print("EMBEDDING PREFLIGHT: no non-journaled target boxes")
            return 0
        failures = []
        for tile in chosen:
            tag = S.tile_id(tile, args.scope)
            try:
                result = audit_recovered_forest(tile, args.scope)
            except AuditError as exc:
                failures.append((tag, str(exc)))
                print("NOT READY %-52s %s" % (tag, exc))
            else:
                print("READY     %-52s %d shards, %d nodes%s" %
                      (tag, len(result["paths"]), result["merged_nodes"],
                       ", merged/discovered" if result["discovered"]
                       else ", dry-merge exact"))
        if failures:
            print("EMBEDDING PREFLIGHT INCOMPLETE: %d/%d boxes not ready" %
                  (len(failures), len(chosen)))
            return 1
        print("EMBEDDING PREFLIGHT GREEN: %d recovered boxes are exact "
              "merge inputs" % len(chosen))
        return 0
    except AuditError as exc:
        print("EMBEDDING PREFLIGHT FAILED:", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
