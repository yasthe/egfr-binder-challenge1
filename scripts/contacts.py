#!/usr/bin/env python3
"""Hat der Lauf den Patch getroffen, und beruehrt er Rest 359?

Checkliste Punkt 5. Hotspots ziehen, aber fesseln nicht - in Lauf B1 und B2
landete der Binder woanders, ohne dass die Metriken das verrieten.

Aufruf:  python contacts.py "designs_C/Accepted/*.pdb" [offset]
         offset weglassen = 299 (Lauf C, Ziel 300-450)
         Lauf B2:  python contacts.py "designs_B2/.../*.pdb" 309
"""
import glob, sys, collections, numpy as np

DES   = sys.argv[1] if len(sys.argv) > 1 else "designs_C/Accepted/*.pdb"
OFF   = int(sys.argv[2]) if len(sys.argv) > 2 else 299
CUT   = 4.5
ANCH  = [320, 323, 344, 355]
VETO  = 359                       # in der Maus Arginin -> Ausschluss
PEN   = [353, 324, 418, 443]                # Abwertung

def chains(path):
    A, B = collections.defaultdict(list), []
    for ln in open(path):
        if not ln.startswith("ATOM") or ln[16] not in " A": continue
        xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        if ln[21] == "A": A[int(ln[22:26]) + OFF].append(xyz)
        elif ln[21] == "B": B.append(xyz)
    return A, np.array(B)

files = sorted(glob.glob(DES))
if not files: sys.exit(f"keine Dateien unter {DES}")
print(f"{len(files)} Designs, Offset {OFF}, Kontaktgrenze {CUT} A\n")

tally = collections.Counter(); rows = []
for f in files:
    A, B = chains(f)
    if not len(B): continue
    touched = set()
    for n, pts in A.items():
        P = np.array(pts)
        if np.linalg.norm(P[:, None, :] - B[None, :, :], axis=2).min() < CUT:
            touched.add(n)
    tally.update(touched)
    hits = [a for a in ANCH if a in touched]
    rows.append((len(hits), VETO in touched, [p for p in PEN if p in touched],
                 sorted(touched), f.split("/")[-1]))

print("--- Wie oft wurde welcher Zielrest beruehrt ---")
for n, c in sorted(tally.items()):
    mark = ""
    if n in ANCH: mark = "  ANKER"
    if n == VETO: mark = "  <== VETO (Maus: Arg)"
    if n in PEN:  mark = "  <== Abwertung"
    if mark or c > len(files) * 0.3:
        print(f"  {n}: {c:4d} / {len(files)}  ({100*c/len(files):4.0f} %){mark}")

print("\n--- Pro Design ---")
print(f"{'Anker':>5} {'Veto359':>8} {'Abwert.':>8}  Datei")
for h, v, p, t, n in sorted(rows, key=lambda r: (-r[0], r[1])):
    print(f"{h:5d} {'JA' if v else '-':>8} {','.join(map(str,p)) or '-':>8}  {n}")

good = [r for r in rows if r[0] >= 2 and not r[1]]
print(f"\nMindestens 2 Anker und kein Kontakt zu 359: {len(good)} / {len(rows)}")
print("Das ist der Filter, der vor allen Metriken greift.")
