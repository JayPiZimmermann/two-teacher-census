"""Sensitivity + fine-grid check of the mixed-sign separated search."""
import math, json, numpy as np, mpmath, separated as S
S.mp_setup(60)
PI = math.pi
res = {"sensitivity_control": [], "fine_mixed": [], "edge_mixed": []}

# 1. sensitivity: the locator must find the known positive-teacher beam family
for (b, y) in [(0.6,-0.2),(1.2,-0.3),(1.9,-0.5),(2.7,-0.8),(1.5708,-0.5)]:
    psi=(y+1)*PI/2; s0,s1=math.sin(psi),math.cos(psi)
    for N in (240, 800):
        c = S.scan_teacher(b, s0, s1, N=N)
        res["sensitivity_control"].append({"beta":b,"y":y,"N":N,"found":len(c),
            "pts":[[round(t0,9),round(t1,9),round(c0,9),round(c1,9)] for t0,t1,c0,c1 in c]})
    print("control beta=%.4f y=%+.2f  N=240:%d  N=800:%d" % (b,y,
        res["sensitivity_control"][-2]["found"], res["sensitivity_control"][-1]["found"]), flush=True)

# 2. fine grid on mixed-sign teachers
tot=0
for b in [0.2,0.5,0.9,1.3,1.5708,1.9,2.3,2.7,3.0]:
    for y in [0.02,0.1,0.25,0.4,0.5,0.6,0.75,0.9,0.98]:
        psi=(y+1)*PI/2; s0,s1=math.sin(psi),math.cos(psi)
        c = S.scan_teacher(b, s0, s1, N=800, rtol=0.6)
        tot += len(c)
        if c:
            res["fine_mixed"].append({"beta":b,"y":y,"found":len(c),
                "pts":[[t0,t1,c0,c1] for t0,t1,c0,c1 in c]})
    print("fine sweep beta=%.4f done, running total %d" % (b,tot), flush=True)
res["fine_mixed_total"] = tot

# 3. edges of the mixed sector (approaching y=0+ and y=1-)
tot2=0
for b in [0.4,1.2,2.0,2.8]:
    for y in [1e-3,1e-2,0.999,0.99]:
        psi=(y+1)*PI/2; s0,s1=math.sin(psi),math.cos(psi)
        c = S.scan_teacher(b, s0, s1, N=600, rtol=0.6)
        tot2 += len(c)
        if c: res["edge_mixed"].append({"beta":b,"y":y,"found":len(c)})
res["edge_mixed_total"] = tot2
print("edge sweep total", tot2)
json.dump(res, open("sep_fine.json","w"), indent=1)
print("control totals:", [r["found"] for r in res["sensitivity_control"]])
print("fine mixed total candidates:", tot, " edge:", tot2)
