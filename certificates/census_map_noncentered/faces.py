"""Pointwise certified witnesses grouped by seven proposed labels.

Agreement among these finite witnesses does not certify that a label is one
global face or that its census is constant between the sampled teachers.
The historical ``face`` record key is a proposed-label identifier, not a
global-topology assertion. ``CERTIFICATE.md`` (status 2026-08-24) is the
authoritative scope ledger.
"""
import json, os, subprocess, sys
import census_cert as X
HERE = os.path.dirname(os.path.abspath(__file__))
FACES = [
 {"id":"F1","name":"same-sign sector y>0","region":"0<beta<pi, 0<y<1",
  "rep":(0.75,0.40),"extra":[(1.50,0.60),(2.60,0.20),(0.50,0.80)]},
 {"id":"F2","name":"mixed sector, band below y=0, OUTSIDE the torque lens",
  "region":"0<beta<pi, w_a(beta)<y<0, outside lens","rep":(0.75,-0.02),
  "extra":[(0.30,-0.02),(1.50,-0.04),(2.20,-0.06)]},
 {"id":"F3","name":"mixed sector, band below y=0, INSIDE the torque lens",
  "region":"beta_c<beta<pi, w_a(beta)<y<0, inside lens","rep":(3.10,-0.06),
  "extra":[(3.10,-0.10),(3.13,-0.08),(3.05,-0.09)]},
 {"id":"F4","name":"mixed sector, middle band, OUTSIDE the torque lens",
  "region":"0<beta<pi, -1-w_a<y<w_a, outside lens","rep":(0.75,-0.40),
  "extra":[(0.30,-0.50),(1.50,-0.60),(2.00,-0.30)]},
 {"id":"F5","name":"mixed sector, middle band, INSIDE the torque lens",
  "region":"beta*_1<beta<pi, inside lens","rep":(2.60,-0.50),
  "extra":[(2.80,-0.50),(3.00,-0.40),(2.40,-0.50)]},
 {"id":"F6","name":"mixed sector, band above y=-1, OUTSIDE the torque lens",
  "region":"0<beta<pi, -1<y<-1-w_a, outside lens","rep":(0.75,-0.98),
  "extra":[(0.30,-0.98),(1.50,-0.96),(2.20,-0.94)]},
 {"id":"F7","name":"mixed sector, band above y=-1, INSIDE the torque lens",
  "region":"beta_c<beta<pi, -1<y<-1-w_a, inside lens","rep":(3.10,-0.94),
  "extra":[(3.10,-0.90),(3.13,-0.92),(3.05,-0.91)]},
]
slot, nslots = int(sys.argv[1]), int(sys.argv[2])
jobs = []
for f in FACES:
    jobs.append((f["id"], f["rep"], True))
    for e in f["extra"]:
        jobs.append((f["id"], e, False))
mine = [j for n, j in enumerate(jobs) if n % nslots == slot]
out = []
for fid, (b, y), isrep in mine:
    try:
        r = X.certified_census(b, y, N=160)
        out.append({"face": fid, "beta": b, "y": y, "is_rep": isrep,
                    "census": r["census"], "n_torque_roots": r["n_torque_roots"],
                    "n_separated": r["n_separated"],
                    "sep_D": [round(s["D_mid"],6) for s in r["separated"]],
                    "sep_schur": [s["schur"] for s in r["separated"]],
                    "certified": r["certified_modulo_locator_completeness"]})
    except Exception as e:
        out.append({"face": fid, "beta": b, "y": y, "error": str(e)})
    with open(os.path.join(HERE, "faces_%d.json" % slot), "w") as fh:
        json.dump(out, fh)
print("slot %d done %d" % (slot, len(out)))
