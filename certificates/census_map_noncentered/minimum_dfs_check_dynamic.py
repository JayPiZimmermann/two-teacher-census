"""Exact dynamically scheduled replay of an ordinary minimum-DFS forest.

This is a separate replay entry point; it does not change the ordinary
certificate format or :mod:`minimum_dfs_check_parallel`.  It uses the same
strict whole-document scan, SHA-256 identity check, deterministic prefix
cover, leaf-method recomputation, assignment audit, and final recomposition
as that checker.  The difference is scheduling:

* the default cut is depth 10, producing smaller replay units;
* units are submitted longest-first as many small queued jobs instead of one
  static node-balanced group per worker; and
* every spawned process loads and rescans the artifact once in its initializer
  and then handles as many queued jobs as the executor assigns to it.  Its
  exact atom cache therefore also survives across those jobs.

The executor's shared queue gives an idle process the next unit, which avoids
pinning a misleadingly cheap-looking but arithmetically expensive prefix to a
single static group.  ``--batch-size`` may amortize queue overhead, but the
default of one is the most responsive schedule for heterogeneous forests.

Usage::

  MINIMUM_REPLAY_WORKERS=4 python3 minimum_dfs_check_dynamic.py ARTIFACT

``MINIMUM_REPLAY_DEPTH`` (default 10) and ``MINIMUM_REPLAY_BATCH`` (default 1)
control calls made by atlas/frontier drivers that pass only the artifact.
"""
import argparse
import concurrent.futures
import json
import multiprocessing
import os

import minimum_dfs_cached as C
import minimum_dfs_check_parallel as P


HERE = os.path.dirname(os.path.abspath(__file__))

# Initialized independently in every spawned replay worker.  Keeping the
# decoded scan process-local lets one worker replay multiple queued groups
# without reopening and rescanning a large artifact for every group.
_WORKER_SCANNED = None
_WORKER_DIGEST = None


def environment_integer(name, default):
    raw = os.environ.get(name, str(default))
    try:
        return int(raw)
    except ValueError as exc:
        raise ValueError("%s is not an integer" % name) from exc


def queued_groups(units, batch_size=1):
    """Return a deterministic longest-first queue covering every unit once."""
    if type(batch_size) is not int or batch_size < 1:
        raise ValueError("batch size must be a positive integer")
    ordered = [key for key, _ in sorted(
        units.items(), key=lambda item: (-item[1]["nodes"], item[0]))]
    if not ordered:
        raise ValueError("forest has no replay units")
    return [ordered[index:index + batch_size]
            for index in range(0, len(ordered), batch_size)]


def initialize_worker(path, digest, frontier_depth, use_cache):
    """Load, hash, and structurally scan once for all jobs in this process."""
    global _WORKER_SCANNED, _WORKER_DIGEST
    document, actual = P.load_hashed(path, digest)
    if actual != digest:
        raise ValueError("worker artifact hash mismatch")
    _WORKER_SCANNED = P.scan_document(document, frontier_depth)
    _WORKER_DIGEST = digest
    if use_cache:
        C.install_atom_cache()


def replay_queued_group(keys):
    """Replay one small job using the worker's persistent decoded scan."""
    if _WORKER_SCANNED is None or _WORKER_DIGEST is None:
        raise RuntimeError("dynamic replay worker was not initialized")
    units = _WORKER_SCANNED["units"]
    if not keys or len(set(keys)) != len(keys) \
            or any(key not in units for key in keys):
        raise ValueError("worker received empty, duplicate, or unknown units")
    counts = P.empty_counts()
    nodes = 0
    leaf_volume = 0.0
    for key in keys:
        part_nodes, part_counts, part_volume = P.replay_unit(
            _WORKER_SCANNED, units[key])
        nodes += part_nodes
        P.add_counts(counts, part_counts)
        leaf_volume += part_volume
    return {"keys": list(keys), "nodes": nodes, "counts": counts,
            "leaf_volume": leaf_volume, "pid": os.getpid(),
            "digest": _WORKER_DIGEST}


def recompose(scanned, results, digest):
    """Audit dynamic assignment and recover the ordinary forest totals."""
    if any(result.get("digest") != digest for result in results):
        raise ValueError("worker replay used the wrong artifact hash")
    P.validate_assignment([result["keys"] for result in results],
                          scanned["units"])
    counts = P.empty_counts()
    subtree_nodes = 0
    leaf_volume = 0.0
    for result in results:
        subtree_nodes += result["nodes"]
        leaf_volume += result["leaf_volume"]
        P.add_counts(counts, result["counts"])
    total_nodes = scanned["prefix_nodes"] + subtree_nodes
    tolerance = 2e-12 * max(1.0, abs(scanned["root_volume"]))
    if total_nodes != scanned["structural_nodes"] \
            or counts != scanned["structural_counts"] \
            or abs(leaf_volume - scanned["root_volume"]) > tolerance:
        raise ValueError("worker results do not recompose the ordinary forest")
    return {"counts": counts, "nodes": total_nodes,
            "leaf_volume": leaf_volume}


def dynamic_replay(path, frontier_depth, workers, batch_size=1,
                   use_cache=True):
    if workers < 1:
        raise ValueError("workers must be positive")
    document, digest = P.load_hashed(path)
    scanned = P.scan_document(document, frontier_depth)
    groups = queued_groups(scanned["units"], batch_size)
    P.validate_assignment(groups, scanned["units"])
    results = []
    worker_count = min(workers, len(groups))
    context = multiprocessing.get_context("spawn")
    with concurrent.futures.ProcessPoolExecutor(
            max_workers=worker_count, mp_context=context,
            initializer=initialize_worker,
            initargs=(path, digest, frontier_depth, use_cache)) as pool:
        pending = [pool.submit(replay_queued_group, group)
                   for group in groups]
        for future in concurrent.futures.as_completed(pending):
            results.append(future.result())
    totals = recompose(scanned, results, digest)
    _, final_digest = P.load_hashed(path, digest)
    if final_digest != digest:
        raise ValueError("artifact changed after dynamic replay")
    return {"document": document, "digest": digest, "scanned": scanned,
            "groups": groups, "results": results,
            "worker_count": worker_count, **totals}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact")
    try:
        worker_default = environment_integer("MINIMUM_REPLAY_WORKERS", 8)
        depth_default = environment_integer("MINIMUM_REPLAY_DEPTH", 10)
        batch_default = environment_integer("MINIMUM_REPLAY_BATCH", 1)
    except ValueError as exc:
        parser.error(str(exc))
    parser.add_argument("--workers", type=int, default=worker_default)
    parser.add_argument("--depth", type=int, default=depth_default)
    parser.add_argument("--batch-size", type=int, default=batch_default)
    parser.add_argument("--no-cache", action="store_true")
    args = parser.parse_args()
    if args.workers < 1 or not 0 <= args.depth <= 32 \
            or args.batch_size < 1:
        parser.error("workers/batch-size must be positive and depth must "
                     "lie in [0,32]")
    path = args.artifact if os.path.isabs(args.artifact) \
        else os.path.join(HERE, args.artifact)
    try:
        result = dynamic_replay(
            path, args.depth, args.workers, args.batch_size,
            not args.no_cache)
    except Exception as exc:
        print("DYNAMIC REPLAY FAILED:", exc)
        return 1
    scanned = result["scanned"]
    pids = {item["pid"] for item in result["results"]}
    print("certificate :", os.path.basename(path))
    print("cells       : [%d,%d) of %d" %
          (scanned["lo_idx"], scanned["hi_idx"],
           scanned["spec"]["n_parts"]))
    print("sha256      :", result["digest"])
    print("bits        :", result["document"]["n_nodes"])
    print("units       : %d at depth <= %d in %d queued groups over %d/%d "
          "persistent workers" %
          (len(scanned["units"]), args.depth, len(result["groups"]),
           len(pids), result["worker_count"]))
    print("leaves      :", json.dumps(result["counts"], sort_keys=True))
    print("coverage    : %.16e leaf volume / %.16e root volume" %
          (result["leaf_volume"], scanned["root_volume"]))
    scope = scanned["spec"].get("scope", "selected")
    claim = ("EVERY angle-map zero is fold-free" if scope == "census"
             else "every selected zero is regular")
    print("DYNAMIC CERTIFICATE VALID AND COMPLETE: %s; every prefix and "
          "leaf recomputed%s." %
          (claim, " with exact atom caching" if not args.no_cache else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
