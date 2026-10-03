#!/usr/bin/env python3
"""Wo koennte ein zusaetzliches Histidin einen Schalter erzeugen?

Ein pKa-Sprung beim Binden entsteht nur, wenn das Histidin (a) ein Carboxylat
des Ziels erreichen kann und (b) dabei dem Wasser entzogen wird. Dieses Skript
listet fuer jede Binderposition beides auf, damit die Mutationskandidaten
nicht geraten werden.

Aufruf:  python scan_his_sites.py <design.pdb> [offset]
         offset: Designrest i -> reifer Rest i+offset   (Lauf B: 309)
"""
import sys, copy, warnings
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
warnings.filterwarnings("ignore")

PDB = sys.argv[1]
OFF = int(sys.argv[2]) if len(sys.argv) > 2 else 309
ACID_O = {"OD1", "OD2", "OE1", "OE2"}
THREE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
         "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
         "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
# in der Maus abweichende Positionen (reife Nummerierung)
MOUSE_DIFF = {306,324,337,340,353,359,366,369,388,390,418,443}

s = PDBParser(QUIET=True).get_structure("x", PDB)[0]
A, B = s["A"], s["B"]
seq = "".join(THREE.get(r.get_resname(), "X") for r in B if r.id[0] == " ")
print(f"{PDB}\nBinder {len(seq)} aa\n{seq}\n")

# --- Zielcarboxylate in Reichweite ---
batoms = np.array([a.coord for a in B.get_atoms()])
acid = []
for r in A:
    if r.get_resname() not in ("ASP", "GLU") or r.id[0] != " ": continue
    ox = [a for a in r if a.get_name() in ACID_O]
    if not ox: continue
    d = min(np.linalg.norm(batoms - a.coord, axis=1).min() for a in ox)
    if d < 10:
        acid.append((r.id[1] + OFF, r.get_resname(), ox, d))
acid.sort(key=lambda x: x[3])
print("Saure Zielreste innerhalb 10 A des Binders:")
for n, nm, ox, d in acid:
    tag = "  (in der Maus abweichend!)" if n in MOUSE_DIFF else ""
    print(f"  {nm}{n}  naechstes Binderatom {d:.2f} A{tag}")
print()

# --- Vergrabenheit jeder Binderposition ---
sr = ShrakeRupley()
sr.compute(s, level="R")
cplx = {r.id[1]: r.sasa for r in B if r.id[0] == " "}
free = copy.deepcopy(s)
for ch in list(free):
    if ch.id != "B": free.detach_child(ch.id)
sr.compute(free, level="R")
alone = {r.id[1]: r.sasa for r in free["B"] if r.id[0] == " "}

# --- Kandidatenliste ---
rows = []
for r in B:
    if r.id[0] != " ": continue
    i = r.id[1]
    cb = r["CB"] if "CB" in r else (r["CA"] if "CA" in r else None)
    if cb is None: continue
    best = None
    for n, nm, ox, _ in acid:
        d = min(np.linalg.norm(cb.coord - a.coord) for a in ox)
        if best is None or d < best[0]: best = (d, f"{nm}{n}")
    if best is None: continue
    bur = 100 * (1 - cplx.get(i, 0) / max(alone.get(i, 1e-9), 1e-9))
    rows.append((best[0], i, THREE.get(r.get_resname(), "X"), best[1], bur,
                 alone.get(i, 0), cplx.get(i, 0)))

rows.sort()
print("Binderpositionen, sortiert nach Reichweite CB -> naechstes Carboxylat")
print("Ein Histidin-Ring reicht vom CB aus etwa 2.5-4.5 A weit.\n")
print(f"{'Pos':>4} {'AS':>3} {'CB->COO':>8} {'Partner':>8} {'vergraben':>10} "
      f"{'SASA frei':>10} {'im Komplex':>11}")
for d, i, aa, part, bur, sf, sc in rows:
    if d > 11: continue
    flag = ""
    if aa != "H" and 3.0 <= d <= 7.5 and bur >= 40: flag = "  <== KANDIDAT"
    if aa == "H": flag = "  (schon His)"
    if aa in "GP": flag += "  [G/P, Ruckgrat]"
    print(f"{i:4d} {aa:>3} {d:8.2f} {part:>8} {bur:9.0f}% {sf:10.1f} {sc:11.1f}{flag}")

print("\nKandidat = nicht-His, CB 3.0-7.5 A vom Carboxylat, mindestens 40 % vergraben.")
print("Positionen an mausdivergenten Zielresten sind oben markiert und zu meiden.")
