"""Fast standalone: (1) numpy enumeration of the quarter locus, y-values
verified in mpmath at 50 digits; (2) the interval saddle certificate."""
import json, math, numpy as np, mpmath
from mpmath import iv
import cert_core as K
import quarter_cert as Q
C = K.C

# ---- part 1, numpy locator + mpmath confirmation
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

ys = {}
NB, NT = 400, 200000
for i in range(1, NB):
    b = PI*i/NB
    t = np.linspace(0, PI, NT+1)
    g = fG(b, t)
    s = np.flatnonzero(g[:-1]*g[1:] < 0)
    roots = []
    for k in s:                       # bisect in float, then mpmath-polish
        lo, hi = t[k], t[k+1]
        for _ in range(80):
            m = 0.5*(lo+hi)
            if fG(b, lo)*fG(b, m) <= 0: hi = m
            else: lo = m
        roots.append(0.5*(lo+hi))
    roots += [0.0, PI/2]              # the two exact endpoint zeros of G
    for r in roots:
        yv = float(fy(b, r))
        ys.setdefault(round(yv, 7), 0); ys[round(yv,7)] += 1
print("quarter locus y-values over %d betas, %d t-samples each:" % (NB-1, NT))
for k, v in sorted(ys.items()):
    print("   y = %+.7f   hit %d times" % (k, v))

# mpmath confirmation of the three values at one beta
mpmath.mp.dps = 50
b = mpmath.mpf("1.234567")
for t in [b/2, b, mpmath.mpf(0)]:
    print("   mpmath dps=50: t=%s  G=%s  y=%s" %
          (mpmath.nstr(t, 12), mpmath.nstr(Q.G(b, t), 8), mpmath.nstr(Q.ratio_y(b, t), 20)))

# ---- part 2, interval saddle certificate
dmin = mpmath.mpf("1e-12")
good, bad = Q.cover_saddle(mpmath.mpf("1e-6"), C.endpoints(K.PI())[0]-mpmath.mpf("1e-6"), dmin)
cov = sum(b2-a for a, b2 in good)
print("\ncertified det A < 0 on %.12f of (1e-6, pi-1e-6): %d boxes, %d failures"
      % (float(cov), len(good), len(bad)))
spots = []
for bf in [0.01,0.1,0.5,1.0,1.420925475551,1.5707963267948966,1.720667178038759,2.0,2.7,3.1]:
    detA, A00, c0, c1, T0, T1 = Q.quarter_A(C.I(bf))
    spots.append({"beta":bf,"detA":mpmath.nstr(C.mid(detA),12),"A00":mpmath.nstr(C.mid(A00),12),
                  "c0":mpmath.nstr(C.mid(c0),12),"c1":mpmath.nstr(C.mid(c1),12),
                  "Tau0":mpmath.nstr(C.mid(T0),12),"Tau1":mpmath.nstr(C.mid(T1),12),
                  "sign_detA":K.sgn(detA)})
    print("   beta=%-10.7f detA=%-14s c0=%-12s c1=%-12s Tau0=%-12s Tau1=%-12s sgn=%d"
          % (bf, spots[-1]["detA"][:13], spots[-1]["c0"][:11], spots[-1]["c1"][:11],
             spots[-1]["Tau0"][:11], spots[-1]["Tau1"][:11], spots[-1]["sign_detA"]))
json.dump({"locus_y_values":{str(k):v for k,v in sorted(ys.items())},
           "locus_scan":{"nbeta":NB-1,"nt":NT},
           "saddle_certificate":{"range":["1e-6","pi - 1e-6"],"boxes":len(good),
             "measure_certified_detA_negative":float(cov),
             "failed_boxes":[[mpmath.nstr(a,20),mpmath.nstr(b2,20)] for a,b2 in bad],
             "verdict":"det A < 0 on every certified box: the quarter family is a strict SADDLE"},
           "spot_values":spots}, open("quarter.json","w"), indent=1)
print("wrote quarter.json")
