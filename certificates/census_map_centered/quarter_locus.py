"""Enumerate the quarter locus {G = 0} and its image in the map, vectorised."""
import json, math, numpy as np, mpmath
PI = math.pi
def fH(t):
    x = np.mod(t, PI); return (PI/2 - x)*np.sin(x)
def fG(b, t):
    return fH(t-b)*fH(t+PI/2) - fH(t)*fH(t+PI/2-b)
def fy(b, t):
    s0, s1 = -fH(t-b), fH(t)
    a = np.arctan2(s0, s1); a = np.mod(a, PI)
    a = np.where(a <= 0, a + PI, a)
    return 2*a/PI - 1

NB, NT = 400, 40000
ys = {}
per_beta = []
for i in range(1, NB):
    b = PI*i/NB
    t = np.linspace(0, PI, NT+1)
    g = fG(b, t)
    k = np.flatnonzero(g[:-1]*g[1:] < 0)
    lo, hi = t[k].copy(), t[k+1].copy()
    glo = fG(b, lo)
    for _ in range(60):                      # vectorised bisection
        m = 0.5*(lo+hi); gm = fG(b, m)
        take = (glo*gm <= 0)
        hi = np.where(take, m, hi)
        lo = np.where(take, lo, m)
        glo = np.where(take, glo, gm)
    roots = list(0.5*(lo+hi)) + [0.0, PI/2]
    yv = [float(fy(b, r)) for r in roots]
    per_beta.append({"beta": b, "n_sign_change_zeros": int(len(k)),
                     "y": sorted(set(round(v, 7) for v in yv))})
    for v in yv:
        key = round(v, 7); ys[key] = ys.get(key, 0) + 1

print("quarter locus over %d betas (%d t-samples each):" % (NB-1, NT))
for k2, v in sorted(ys.items()):
    print("   y = %+.7f   hit %d times" % (k2, v))
print("distinct y values:", sorted(ys))
nz = sorted(set(p["n_sign_change_zeros"] for p in per_beta))
print("sign-change zero counts of G per beta:", nz)

mpmath.mp.dps = 50
def mH(t):
    x = mpmath.fmod(t, mpmath.pi)
    if x < 0: x += mpmath.pi
    return (mpmath.pi/2 - x)*mpmath.sin(x)
def mG(b, t): return mH(t-b)*mH(t+mpmath.pi/2) - mH(t)*mH(t+mpmath.pi/2-b)
b = mpmath.mpf("1.234567891")
print("\nmpmath dps=50 confirmation at beta =", mpmath.nstr(b, 12))
for nm, t in [("beta/2", b/2), ("beta/2+pi/2", b/2+mpmath.pi/2), ("beta", b),
              ("beta+pi/2", b+mpmath.pi/2), ("0", mpmath.mpf(0)), ("pi/2", mpmath.pi/2)]:
    print("   t = %-14s G = %s" % (nm, mpmath.nstr(mG(b, t), 8)))

json.dump({"scan": {"nbeta": NB-1, "nt": NT},
           "y_values": {str(k2): v for k2, v in sorted(ys.items())},
           "distinct_y": sorted(ys),
           "sign_change_zero_counts": nz,
           "per_beta_sample": per_beta[::40]},
          open("quarter_locus.json", "w"), indent=1)
print("wrote quarter_locus.json")
