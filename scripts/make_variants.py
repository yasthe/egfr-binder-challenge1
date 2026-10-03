#!/usr/bin/env python3
"""Histidin-Varianten des Kandidaten bauen und neu packen.

Das Rueckgrat bleibt unveraendert; mutiert werden nur Seitenketten, und die
Umgebung wird in 8 A neu gepackt. Grundlage ist scan_his_sites.py:

  S21  3.79 A zu D344, 95 % vergraben
  Q43  5.75 A zu D323, 72 % vergraben
  Q24  6.22 A zu D344, 61 % vergraben
  M11  7.19 A zu D355, 87 % vergraben

Aufruf:  python make_variants.py <design.pdb>
Ergebnis: variants/<name>.pdb  plus eine Tabelle der Rosetta-Energien.
"""
import os, sys, itertools, warnings
warnings.filterwarnings("ignore")

BASE = sys.argv[1]
OUT  = "variants"
SITES = {11: "M", 21: "S", 24: "Q", 43: "Q"}     # Binder-Nummerierung, Kette B

from pyrosetta import init, pose_from_pdb, get_fa_scorefxn
from pyrosetta.toolbox import mutate_residue

init("-mute all -ex1 -ex2aro")
sfx = get_fa_scorefxn()
os.makedirs(OUT, exist_ok=True)

base = pose_from_pdb(BASE)
info = base.pdb_info()
e0 = sfx(base)
print(f"Ausgangsstruktur: {base.total_residue()} Reste, Score {e0:.1f}\n")

combos = []
for k in range(0, 4):
    for c in itertools.combinations(sorted(SITES), k + 1):
        combos.append(c)
combos = [()] + combos           # Original zuerst

print(f"{'Name':<22} {'Mutationen':<22} {'Score':>9} {'dScore':>8}")
rows = []
for c in combos:
    name = "parent" if not c else "_".join(f"{SITES[i]}{i}H" for i in c)
    p = base.clone()
    ok = True
    for i in c:
        idx = info.pdb2pose("B", i)
        if idx == 0:
            print(f"  {name}: Rest B{i} nicht gefunden"); ok = False; break
        if p.residue(idx).name1() != SITES[i]:
            print(f"  {name}: B{i} ist {p.residue(idx).name1()}, "
                  f"erwartet {SITES[i]}"); ok = False; break
        mutate_residue(p, idx, "H", pack_radius=8.0, pack_scorefxn=sfx)
    if not ok: continue
    e = sfx(p)
    path = os.path.join(OUT, f"{name}.pdb")
    p.dump_pdb(path)
    rows.append((e - e0, name, e))
    print(f"{name:<22} {','.join(str(i) for i in c) or '-':<22} {e:9.1f} {e-e0:+8.1f}")

print(f"\n{len(rows)} Varianten in {OUT}/")
print("\nNach Energieaenderung sortiert (positiv = schlechter gepackt):")
for d, n, e in sorted(rows):
    warn = "   <== Packungsproblem" if d > 15 else ""
    print(f"  {d:+8.1f}  {n}{warn}")
print("\nEine deutlich positive Energieaenderung heisst, dass das Histidin")
print("sterisch nicht hineinpasst - solche Varianten nicht einreichen,")
print("auch wenn PROPKA einen guten Schalter vorhersagt.")
