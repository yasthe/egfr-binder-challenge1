#!/usr/bin/env python3
"""Zuschnitt des 6ARU-Ziels fuer den Patch-C-Lauf.
Sucht den kleinsten Ausschnitt, dessen Schnittkanten weit genug von den
Ankerresten entfernt liegen, und prueft ihn, bevor etwas gestartet wird."""
import sys, warnings, numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
warnings.filterwarnings("ignore")

PDB     = "6ARU.pdb"
ANCHOR  = [323, 344, 355]          # Patch C
EXPECT  = {320:"GLU",323:"ASP",344:"ASP",355:"ASP",468:"SER",467:"ILE"}
MAXLEN  = 175                      # Reste Ziel; + Binder <= ~230 passen in 12 GB
MINEDGE = 18.0                     # A, Mindestabstand Anker -> Schnittkante

st = PDBParser(QUIET=True).get_structure("x", PDB)[0]
A  = st["A"]
res = {r.id[1]: r for r in A if r.id[0] == " "}
nums = sorted(res)
print(f"Kette A: {len(nums)} Reste, {nums[0]}-{nums[-1]}")

# --- 1. Nummerierung verifizieren -----------------------------------------
bad = [f"{n}: erwartet {w}, gefunden {res[n].get_resname() if n in res else 'FEHLT'}"
       for n, w in EXPECT.items() if n not in res or res[n].get_resname() != w]
if bad:
    print("NUMMERIERUNG STIMMT NICHT:"); [print("  ", b) for b in bad]; sys.exit(1)
print("Nummerierung verifiziert (reife Zaehlung, wie 1YY9).")

co = {n: np.array([a.coord for a in r]) for n, r in res.items()}
def dmin(i, j):
    return np.linalg.norm(co[i][:, None, :] - co[j][None, :, :], axis=2).min()

# --- 2. Zuschnitt suchen ---------------------------------------------------
best = None
for s in range(280, min(ANCHOR) + 1):
    if s not in res: continue
    for e in range(max(ANCHOR), 500):
        if e not in res: continue
        span = [n for n in nums if s <= n <= e]
        if len(span) > MAXLEN: break
        edge = [n for n in span if n <= s + 4 or n >= e - 4]
        sc = min(dmin(a, n) for a in ANCHOR for n in edge)
        if best is None or (sc, -len(span)) > (best[0], -best[3]):
            best = (sc, s, e, len(span))
sc, S, E, L = best
print(f"\nGewaehlter Ausschnitt: {S}-{E}  ({L} Reste)")
print(f"Kleinster Abstand Anker -> Schnittkante: {sc:.1f} A  (Ziel >= {MINEDGE})")
if sc < MINEDGE:
    print("WARNUNG: Kante zu nah am Patch. Nicht starten, erst melden.")

# --- 3. Luecken im Ausschnitt ---------------------------------------------
span = [n for n in nums if S <= n <= E]
gaps = [(span[i], span[i+1]) for i in range(len(span)-1) if span[i+1] != span[i]+1]
print(f"Luecken im Ausschnitt: {gaps if gaps else 'keine'}")
if gaps: print("  -> jede Luecke ist eine zusaetzliche kuenstliche Kante.")

# --- 4. Verdachtszone: was lag an den Kanten? (Checkliste Punkt 3) --------
cut = [n for n in nums if n < S or n > E]
near = sorted({n for n in cut for a in span[:5] + span[-5:] if dmin(n, a) < 10})
print(f"Weggeschnittene Reste innerhalb 10 A der Kanten: {len(near)}")
print("  ", near[:30], "..." if len(near) > 30 else "")

# --- 5. Zugaenglichkeit der Anker: voll vs. getrimmt ----------------------
sr = ShrakeRupley()
sr.compute(st, level="R")
full = {n: res[n].sasa for n in ANCHOR}
import copy
tr = copy.deepcopy(st)
for ch in list(tr):
    if ch.id != "A": tr.detach_child(ch.id)
for r in [r for r in tr["A"] if r.id[0] != " " or not (S <= r.id[1] <= E)]:
    tr["A"].detach_child(r.id)
sr.compute(tr, level="R")
trm = {r.id[1]: r.sasa for r in tr["A"]}
print("\nSASA der Anker (A^2):")
for n in ANCHOR:
    print(f"  {res[n].get_resname()}{n}:  voll+Fab {full[n]:6.1f}   getrimmt {trm.get(n,0):6.1f}")

# --- 6. freie Flaeche im Umkreis, am VOLLSTAENDIGEN Rezeptor --------------
sr.compute(st, level="R")
env = [r for r in A if r.id[0] == " " and
       min(dmin(r.id[1], a) for a in ANCHOR) < 15]
print(f"\nFreie Flaeche im 15-A-Umkreis (vollstaendiger Rezeptor, ohne Fab-Blockade "
      f"gerechnet mit Fab): {sum(r.sasa for r in env):.0f} A^2 ueber {len(env)} Reste")
print("  Faustwert: ein 50-Reste-Binder braucht 800-1200 A^2 zusammenhaengend.")

# --- 7. schreiben ----------------------------------------------------------
out = "EGFR_6aru_patchC.pdb"
with open(PDB) as fh, open(out, "w") as o:
    for ln in fh:
        if ln.startswith("ATOM") and ln[21] == "A" and ln[16] in " A":
            if S <= int(ln[22:26]) <= E:
                o.write(ln[:16] + " " + ln[17:])
    o.write("END\n")
print(f"\ngeschrieben: {out}")
print(f"OFFSET fuer spaetere Analysen: BindCraft-Rest i  ->  reifer Rest i+{S-1}")
print(f"HOTSPOTS fuer die JSON: " + ",".join(f"A{a}" for a in ANCHOR))
