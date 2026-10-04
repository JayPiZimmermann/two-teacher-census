"""Regenerate the README's replayed-witness table FROM THE REPLAY LOGS.

A row is emitted only when its replay reported VALID AND COMPLETE with zero
failed leaves, zero undecided leaves and certified class disjointness; the
family count, leaf count and boundary margin are read out of the log rather
than typed.  Witnesses whose replay has not finished are listed and omitted.

Run from this directory: `python3 fill_witness_table.py`
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ORDER = [
 ("F2", "C_F2_0p3_m0p02",  "(0.3, −0.02)"),
 ("F2", "C_F2_0p75_m0p02", "(0.75, −0.02)"),
 ("F2", "C_F2_1p5_m0p04",  "(1.5, −0.04)"),
 ("F2", "C_F2_2p2_m0p06",  "(2.2, −0.06)"),
 ("F4", "C_F4_0p3_m0p5",   "(0.3, −0.5)"),
 ("F4", "C_F4_0p75_m0p4",  "(0.75, −0.4)"),
 ("F4", "C_F4_1p5_m0p6",   "(1.5, −0.6)"),
 ("F4", "C_F4_2_m0p3",     "(2.0, −0.3)"),
 ("F5", "C_F5_2p4_m0p5",   "(2.4, −0.5)"),
 ("F5", "C_F5_2p6_m0p5",   "(2.6, −0.5)"),
 ("F5", "C_F5_2p8_m0p5",   "(2.8, −0.5)"),
 ("F5", "C_F5_3_m0p42",    "(3.0, −0.42)"),
 ("F5", "C_F5_3_m0p4",     "(3.0, −0.4)"),
 ("F5", "C_F5_3_m0p5",     "(3.0, −0.5)"),
 ("F5", "C_F5_3_m0p58",    "(3.0, −0.58)"),
 ("F5", "C_F5_3_m0p60",    "(3.0, −0.60)"),
 ("F3", "C_F3_3p08_m0p08",  "(3.08, −0.08)"),
 ("F6", "C_F6_0p75_m0p98", "(0.75, −0.98)"),
 ("F6", "C_F6_1p5_m0p96",  "(1.5, −0.96)"),
 ("F6", "C_F6_2p2_m0p94",  "(2.2, −0.94)"),
 ("F7", "C_F7_3p08_m0p92",  "(3.08, −0.92)"),
]
rows, bad = [], []
for face, tag, label in ORDER:
    p = os.path.join(HERE, "replay_%s.log" % tag)
    txt = open(p).read() if os.path.exists(p) else ""
    if "CERTIFICATE VALID AND COMPLETE" not in txt:
        bad.append((tag, "not valid-and-complete"))
        continue
    mf = re.search(r"SEPARATED FAMILIES: (\d+)", txt)
    mm = re.search(r"faces = ([0-9.e+-]+)", txt)
    ml = re.search(r"leaves\s+: (\d+) claimed, (\d+) verified, (\d+) FAILED, (\d+) undecided", txt)
    md = re.search(r"disjointness certified: (\w+)", txt)
    if not (mf and mm and ml and md):
        bad.append((tag, "missing lines")); continue
    if ml.group(3) != "0" or ml.group(4) != "0" or md.group(1) != "True":
        bad.append((tag, "failed/undecided/overlap")); continue
    rows.append("| %s | `%s` | %s | %s | %.2e |"
                % (face, label, mf.group(1), format(int(ml.group(1)), ",").replace(",", " "),
                   float(mm.group(1))))
if bad:
    print("NOT YET REPLAYED (omitted from the table):", bad)
table = ("| face | `(beta, y)` | families | leaves | boundary margin |\n"
         "|---|---|---|---|---|\n" + "\n".join(rows))
p = os.path.join(HERE, "README.md")
s = open(p).read()
start = s.index("| face | `(beta, y)` | families |")
end = s.index("\n\n### F5 is not one count-homogeneous face")
s = s[:start] + table + s[end:]
open(p, "w").write(s)
print(table)
print("\nrows emitted:", len(rows))
