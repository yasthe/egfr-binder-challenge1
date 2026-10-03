#!/usr/bin/env python3
"""Sequenzen fuer das vorhandene Rueckgrat neu wuerfeln, mit Histidin-Bias.

Das Rueckgrat des validierten Kandidaten bleibt unveraendert; ProteinMPNN
entwirft die Binderoberflaeche neu und darf dabei die Umgebung eines
Histidins mitgestalten - der Unterschied zu den Punktmutanten, bei denen
das Histidin in eine fuer etwas anderes optimierte Nachbarschaft gesetzt
wurde.

Mehrere Biasstufen, weil ein zu starker Bias nur die Histidinzahl hochtreibt
statt Geometrie zu erzeugen. Stufe 0.0 laeuft als Kontrolle mit.

Aufruf:  python mpnn_his.py <design.pdb> [anzahl_pro_stufe] [top_pro_stufe]
Ergebnis: mpnn_variants/*.pdb  zum Weiterreichen an ph_scan.py und check_full.py
"""
import os, sys, warnings
import numpy as np
warnings.filterwarnings("ignore")

PDB   = sys.argv[1]
NUM   = int(sys.argv[2]) if len(sys.argv) > 2 else 40     # Sequenzen je Stufe
TOP   = int(sys.argv[3]) if len(sys.argv) > 3 else 8      # davon weiterverfolgt
BIAS  = [0.0, 0.5, 1.0, 1.5]
TEMP  = 0.2
OUT   = "mpnn_variants"
MINH, MAXH = 1, 3          # Histidine im Binder

import jax
if not hasattr(jax, "tree_map"):
    jax.tree_map = jax.tree_util.tree_map
from colabdesign.mpnn import mk_mpnn_model
from colabdesign.shared.model import aa_order

os.makedirs(OUT, exist_ok=True)
m = mk_mpnn_model(model_name="v_48_020")
m.prep_inputs(pdb_filename=PDB, chain="A,B", fix_pos="A", rm_aa="C")
tlen, blen = m._lengths
base_bias = m._inputs["bias"].copy()
print(f"Ziel {tlen} Reste, Binder {blen} Reste\n")

# --- Ausgangssequenz ---
from pyrosetta import init, pose_from_pdb, get_fa_scorefxn
from pyrosetta.toolbox import mutate_residue
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    RestrictToRepacking, InitializeFromCommandline)
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

init("-mute all -ex1 -ex2aro")
sfx = get_fa_scorefxn()
base = pose_from_pdb(PDB)
info = base.pdb_info()
parent = "".join(base.residue(info.pdb2pose("B", i)).name1()
                 for i in range(1, blen + 1))
e0 = sfx(base)
print(f"Ausgangssequenz ({len(parent)} aa, Score {e0:.1f}):\n{parent}\n")

tf = TaskFactory()
tf.push_back(InitializeFromCommandline())
tf.push_back(RestrictToRepacking())
packer = PackRotamersMover(sfx)
packer.task_factory(tf)

# --- sammeln ---
cand = {}
for b in BIAS:
    m._inputs["bias"] = base_bias.copy()
    if b:
        m._inputs["bias"][tlen:, aa_order["H"]] += b
    out = m.sample(num=NUM, batch=1, temperature=TEMP)
    rows = []
    for s, sc in zip(np.atleast_1d(out["seq"]), np.atleast_1d(out["score"])):
        seq = str(s).split("/")[-1]
        if len(seq) != blen or "C" in seq:
            continue
        nh = seq.count("H")
        if not (MINH <= nh <= MAXH):
            continue
        ident = sum(a == c for a, c in zip(parent, seq)) / blen * 100
        rows.append((float(sc), seq, nh, ident))
    rows.sort()
    print(f"Bias {b:.1f}: {len(rows)} von {NUM} brauchbar "
          f"(1-3 His, kein Cys), beste {min(TOP, len(rows))} weiterverfolgt")
    for i, (sc, seq, nh, ident) in enumerate(rows[:TOP]):
        name = f"b{str(b).replace('.','')}_{i:02d}"
        cand[name] = (seq, sc, nh, ident, b)
        print(f"   {name}  score {sc:.3f}  His {nh}  Identitaet {ident:.0f}%")
print()

# --- auffaedeln und packen ---
print(f"{'Name':<12} {'His':>4} {'Ident':>6} {'Mutationen':>11} {'Score':>9} {'dScore':>8}")
res = []
for name, (seq, sc, nh, ident, b) in cand.items():
    p = base.clone()
    nmut = 0
    for i, (a, c) in enumerate(zip(parent, seq), start=1):
        if a == c:
            continue
        idx = info.pdb2pose("B", i)
        mutate_residue(p, idx, c, pack_radius=0.0, pack_scorefxn=sfx)
        nmut += 1
    packer.apply(p)
    e = sfx(p)
    p.dump_pdb(os.path.join(OUT, f"{name}.pdb"))
    res.append((e - e0, name, nh, ident, nmut))
    print(f"{name:<12} {nh:4d} {ident:5.0f}% {nmut:11d} {e:9.1f} {e-e0:+8.1f}")

print(f"\n{len(res)} Strukturen in {OUT}/")
print("\nNach Energieaenderung sortiert:")
for d, n, nh, ident, nmut in sorted(res):
    print(f"  {d:+8.1f}  {n}  ({nh} His, {ident:.0f}% identisch, {nmut} Mutationen)")
print("\nNaechste Schritte:")
print("  python ph_scan.py   \"mpnn_variants/*.pdb\"")
print("  python check_full.py \"mpnn_variants/*.pdb\" structures/EGFR_d3_trim.pdb 309")
print("\nAchtung: diese Sequenzen wurden nie von AlphaFold2 gefaltet.")
print("Was ueberlebt, muss am Ende neu vorhergesagt werden.")
