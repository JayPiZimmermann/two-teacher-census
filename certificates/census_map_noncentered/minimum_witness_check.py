"""Replay a selected-minimum or full-census witness-count certificate.

The manifest joins two independently checkable coverings at one teacher:

* a small-gap ``minimum_dfs`` forest whose every leaf excludes the angle map;
* a complete ``ivcert`` enumeration above that collar.

The checker replays both trees, keeps precisely the enumerated roots in the
named physical strip, certifies their family deduplication, and evaluates the
two Cramer selectors on every surviving root enclosure.  Thus the reported
selected count is derived from interval signs rather than from a stored
``schur`` or census label.  In ``census`` scope the claim is instead the total
number of certified angle-map zero classes in the target strip.

Usage: python3 minimum_witness_check.py minimum_witness_<tag>.json
"""
import json
import os
import sys

import mpmath

import census_cert as X
import face_bb_signlaw as F
import ivcert_check as I
import minimum_dfs as M
from minimum_dfs_check import Reader

HERE = os.path.dirname(os.path.abspath(__file__))
MAP_VERDICTS = {"MAP_PLAIN", "MAP_CENTERED", "MAP_AUX", "MAP_COLLAR"}


def local_path(path):
    return path if os.path.isabs(path) else os.path.join(HERE, path)


def replay_collar(path, teacher, strip, enum_delta):
    """Replay the map-only small-gap covering at the witness teacher."""
    with open(local_path(path)) as handle:
        doc = json.load(handle)
    spec = doc["spec"]
    beta, y = teacher["beta"], teacher["y"]
    expected = {
        "b0": beta, "b1": beta, "y0": y, "y1": y,
        "seam": strip["seam"], "delta": strip["delta"],
        "dmax": enum_delta, "coord": "absolute",
    }
    if any(spec.get(k) != v for k, v in expected.items()):
        raise ValueError("collar spec does not match the witness manifest")
    if not M.partition_covers(spec):
        raise ValueError("collar root partition does not cover its spec")
    if doc["lo_idx"] != 0 or doc["hi_idx"] != spec["n_parts"]:
        raise ValueError("collar certificate is not a full forest")

    reader = Reader(doc["bits_b64"], doc["n_nodes"])
    roots = M.root_cells(spec)
    leaves = 0
    counts = {name: 0 for name in MAP_VERDICTS}

    def walk(box, depth=0):
        nonlocal leaves
        if depth > 500:
            raise ValueError("collar tree depth exceeds 500")
        if reader.read() == 1:
            walk(M.child(spec, box, 0), depth + 1)
            walk(M.child(spec, box, 1), depth + 1)
            return
        claimed = M.NAME.get(reader.read(3))
        got = M.verdict(spec, *box, target=claimed)
        if claimed not in MAP_VERDICTS or got != claimed:
            raise ValueError("collar leaf does not replay map-only: %s / %s"
                             % (claimed, got))
        counts[claimed] += 1
        leaves += 1

    for root in roots:
        walk(root)
    if reader.pos != doc["n_nodes"]:
        raise ValueError("collar bitstream has unconsumed or truncated nodes")
    return leaves, counts


def unpack_enumeration(doc):
    spec = doc["spec"]
    if "bits_b64" not in doc:
        return doc["leaves"]
    import ivcert_pack as P
    return P.unpack(doc["bits_b64"], doc["n_nodes"],
                    spec["n_root_cells"])


def replay_enumeration(path, teacher):
    """Replay every enumeration leaf and return certified root enclosures."""
    with open(local_path(path)) as handle:
        doc = json.load(handle)
    spec = doc["spec"]
    if spec["beta"] != teacher["beta"] or spec["y"] != teacher["y"] \
            or spec.get("beta_hi") is not None \
            or spec.get("y_hi") is not None:
        raise ValueError("enumeration is not at the manifest's point teacher")
    leaves = unpack_enumeration(doc)
    Tc, roots = I.root_cells(spec)
    if len(roots) != spec["n_root_cells"]:
        raise ValueError("enumeration root cells disagree with its spec")
    by_root = {}
    prefixes = {}
    for key, verdict in leaves.items():
        ri_text, path_text = key.split(":", 1)
        ri = int(ri_text)
        by_root.setdefault(ri, {})[path_text] = verdict
        pref = prefixes.setdefault(ri, set())
        for length in range(len(path_text)):
            pref.add(path_text[:length])

    beta_lo = float(I.J.lo(Tc.B))
    beta_hi = float(I.J.hi(Tc.B))
    Ybox = I.J.I(str(spec["y"]))
    s1box = X.masses_at(Ybox)[1]
    positive = float(I.J.lo(s1box)) > 0.0
    mixed = float(I.J.hi(s1box)) < 0.0
    if not (positive or mixed):
        raise ValueError("enumeration teacher is on the mass-zero wall")

    used = set()
    enclosures = []
    verified = 0

    def walk(ri, box, path_text, depth=0):
        nonlocal verified
        if depth > 500:
            raise ValueError("enumeration tree depth exceeds 500")
        verdict = by_root.get(ri, {}).get(path_text)
        if verdict is None:
            if path_text not in prefixes.get(ri, set()):
                raise ValueError("enumeration tree has a missing child")
            walk(ri, I.child(box, 0), path_text + "0", depth + 1)
            walk(ri, I.child(box, 1), path_text + "1", depth + 1)
            return
        if verdict.startswith("UNDECIDED"):
            raise ValueError("enumeration contains an undecided leaf")
        ok, encl = I.verify_leaf(Tc, box, verdict, beta_lo, beta_hi, positive)
        if not ok:
            raise ValueError("enumeration leaf does not replay: %d:%s %s"
                             % (ri, path_text, verdict))
        used.add("%d:%s" % (ri, path_text))
        verified += 1
        if verified % 10000 == 0:
            print("  replayed %d enumeration leaves" % verified, flush=True)
        if encl is not None:
            S, D = encl
            enclosures.append({
                "th1_mid": float(I.J.mid(S)),
                "th0_mid": float(I.J.mid(S + D)),
                "D_mid": float(I.J.mid(D)),
                "width": max(I.J.width(S), I.J.width(D)),
                "encl": [I.J.lo(S), I.J.hi(S), I.J.lo(D), I.J.hi(D)],
                "schur": "unused",
            })

    sys.setrecursionlimit(10000)
    for ri, root in enumerate(roots):
        walk(ri, root, "")
    if used != set(leaves):
        raise ValueError("enumeration contains unreachable leaf records")
    if verified != spec["n_leaves"]:
        raise ValueError("enumeration leaf count disagrees with its spec")
    return spec, Tc, enclosures, verified


def in_target_strip(enclosures, enum_spec, strip):
    """Keep roots in the target strip and reject all boundary ambiguity."""
    seam = mpmath.mpf(repr(strip["seam"]))
    seam_hi = seam + I.J.hi(2 * X.PIv)
    dmax = mpmath.mpf(repr(strip["dmax"]))
    enum_delta = mpmath.mpf(repr(enum_spec["delta"]))
    inside = []
    outside = 0
    for sol in enclosures:
        slo, shi, dlo, dhi = sol["encl"]
        if not (seam < slo and shi < seam_hi and enum_delta < dlo):
            raise ValueError("root enclosure touches an enumeration boundary")
        if dhi < dmax:
            inside.append(sol)
        elif dmax < dlo:
            outside += 1
        else:
            raise ValueError("root enclosure touches the target dmax face")
    return inside, outside


def selector_classification(classes, teacher):
    """Classify every root by strict interval signs of the Lean selectors."""
    B = M.mk(teacher["beta"], teacher["beta"])
    Y = M.mk(teacher["y"], teacher["y"])
    selected = 0
    rows = []
    for root in classes:
        slo, shi, dlo, dhi = root["encl"]
        S, D = M.mk(slo, shi), M.mk(dlo, dhi)
        t00, _, _, det = M.cramer_selectors(B, Y, S, D, False)
        s00, sdet = X.sgn(t00), X.sgn(det)
        if s00 > 0 and sdet > 0:
            label = "SELECTED"
            selected += 1
        elif s00 < 0 or sdet < 0:
            label = "UNSELECTED"
        else:
            raise ValueError("selector sign is unresolved on a root enclosure")
        rows.append((root, label, s00, sdet))
    return selected, rows


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    with open(local_path(sys.argv[1])) as handle:
        manifest = json.load(handle)
    if manifest.get("format") != "minimum-witness-v1":
        print("REPLAY FAILED: unknown manifest format")
        return 1
    teacher = manifest["teacher"]
    strip = manifest["strip"]
    try:
        with open(local_path(manifest["enumeration"])) as handle:
            enum_delta = json.load(handle)["spec"]["delta"]
        if not (0 < strip["delta"] < enum_delta < strip["dmax"]
                < mpmath.pi):
            raise ValueError("manifest strip/collar ordering is inadmissible")
        collar_leaves, collar_counts = replay_collar(
            manifest["collar"], teacher, strip, enum_delta)
        enum_spec, Tc, enclosures, enum_leaves = replay_enumeration(
            manifest["enumeration"], teacher)
        if enum_spec["seam"] != strip["seam"]:
            raise ValueError("enumeration seam differs from target strip")
        target, outside = in_target_strip(enclosures, enum_spec, strip)
        dedup = F.certified_dedup(Tc, target, teacher["beta"])
        if not dedup["disjointness_certified"] \
                or dedup["overlapping_class_pairs"]:
            raise ValueError("target root classes are not disjoint")
        if dedup["n_merges"]:
            raise ValueError("target strip has duplicate labelled root boxes")
        selected, rows = selector_classification(dedup["classes"], teacher)
        scope = manifest.get("scope", "selected")
        if scope == "selected":
            count = selected
            expected = manifest["expected_selected"]
            count_name = "selected witness count"
        elif scope == "census":
            count = dedup["n_classes"]
            expected = manifest["expected_zeros"]
            count_name = "angle-map witness zero count"
        else:
            raise ValueError("unknown witness scope %s" % scope)
        if count != expected:
            raise ValueError("%s %d differs from expected %d"
                             % (count_name, count, expected))
    except (KeyError, ValueError, IndexError) as exc:
        print("REPLAY FAILED:", exc)
        return 1

    print("manifest     : %s" % os.path.basename(sys.argv[1]))
    print("teacher      : beta=%s y=%s" % (teacher["beta"], teacher["y"]))
    print("strip        : seam=%s D=[%s,%s]" %
          (strip["seam"], strip["delta"], strip["dmax"]))
    print("small collar : %d map-only leaves %s" %
          (collar_leaves, json.dumps(collar_counts, sort_keys=True)))
    print("enumeration  : %d leaves, %d raw root enclosures" %
          (enum_leaves, len(enclosures)))
    print("target roots : %d raw, %d certified classes; %d fold copies outside"
          % (len(target), dedup["n_classes"], outside))
    for root, label, s00, sdet in sorted(rows,
                                         key=lambda item: item[0]["D_mid"]):
        print("  th1=%.9f D=%.9f  T00 sign=%+d det sign=%+d  %s" %
              (root["th1_mid"], root["D_mid"], s00, sdet, label))
    print("scope        : %s" % scope)
    print("CERTIFICATE VALID AND COMPLETE: %s = %d" % (count_name, count))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
