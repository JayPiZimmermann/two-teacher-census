"""Compare the CERTIFIED census against the SHIPPED censusSignature on a set of
points, and classify every mismatch by its cause."""
import glob
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load(tag=""):
    pts = []
    for f in sorted(glob.glob(os.path.join(HERE, "band%s_*.json" % tag))):
        with open(f) as fh:
            pts += json.load(fh)
    pts = [p for p in pts if "census" in p]
    pts.sort(key=lambda p: (p["beta"], p["y"]))
    return pts


def shipped(pts):
    q = json.dumps([[p["beta"], p["y"]] for p in pts])
    r = subprocess.run(["node", os.path.join(HERE, "bandJ.js"), q],
                       capture_output=True, text=True, check=True)
    return json.loads(r.stdout)


def count_traps(c):
    it = c.split(" | ")
    return {
        "coincident": sum(1 for x in it if x.startswith("coincident:trap")),
        "separate": sum(1 for x in it if x.startswith("separate:trap")),
    }


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else ""
    pts = load(tag)
    sh = shipped(pts)
    n = len(pts)
    match = 0
    mismatch = []
    for p, s in zip(pts, sh):
        if p["census"] == s["census"]:
            match += 1
        else:
            a, b = count_traps(p["census"]), count_traps(s["census"])
            cause = []
            if a["coincident"] != b["coincident"]:
                cause.append("coincident_%d_vs_%d" % (a["coincident"],
                                                      b["coincident"]))
            if a["separate"] != b["separate"]:
                cause.append("separate_%d_vs_%d" % (a["separate"], b["separate"]))
            mismatch.append({"beta": p["beta"], "y": p["y"],
                             "cert": p["census"], "ship": s["census"],
                             "cause": "+".join(cause) or "other",
                             "n_torque_roots": p["n_torque_roots"],
                             "certified": p["certified"],
                             "sep_D": p["sep_D"], "sep_schur": p["sep_schur"]})
    import collections
    print("points %d   match %d   mismatch %d" % (n, match, len(mismatch)))
    print("certified flag true on %d of %d"
          % (sum(1 for p in pts if p["certified"]), n))
    print("distinct CERTIFIED census values: %d"
          % len({p["census"] for p in pts}))
    print("distinct SHIPPED   census values: %d"
          % len({s["census"] for s in sh}))
    print("mismatch causes:", dict(collections.Counter(
        m["cause"] for m in mismatch)))
    print("mismatch beta histogram:", dict(collections.Counter(
        round(m["beta"], 3) for m in mismatch)))
    for m in mismatch[:12]:
        print("  beta=%.6f y=%+.4f  %s" % (m["beta"], m["y"], m["cause"]))
        print("     CERT: %s" % m["cert"])
        print("     SHIP: %s" % m["ship"])
    out = {"n_points": n, "n_match": match, "n_mismatch": len(mismatch),
           "distinct_certified": sorted({p["census"] for p in pts}),
           "distinct_shipped": sorted({s["census"] for s in sh}),
           "mismatch_causes": dict(collections.Counter(m["cause"]
                                                       for m in mismatch)),
           "mismatches": mismatch,
           "points": [dict(p, shipped=s["census"]) for p, s in zip(pts, sh)]}
    with open(os.path.join(HERE, "band_compare%s.json" % tag), "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote band_compare%s.json" % tag)


if __name__ == "__main__":
    main()
