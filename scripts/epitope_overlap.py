#!/usr/bin/env python3
"""Schneidet das Epitop jedes Designs mit dem von Cetuximab beruehrten
Rezeptorbereich.

Eine sterische Ueberlappung mit dem Fab belegt Konkurrenz um die Bindung, nicht
ein gemeinsames Epitop. Dieses Skript bestimmt die vom Antikoerper beruehrten
Rezeptorreste und schneidet sie mit denen des Binders.

Aufruf:
  python epitope_overlap.py einreichung_top20.csv [--ref 6ARU.pdb]
         [--receptor A] [--antibody B,C] [--patchdir <Verzeichnis>]

Aus jedem Verzeichnis aufrufbar; Zielausschnitte und Referenzstruktur werden
in den ueblichen Ablageorten gesucht.

Ausgabe je Design: eigenes Epitop, Cetuximab-Epitop, Schnittmenge, Anteil.
"""
import os, sys
import numpy as np
import pandas as pd

CSV = sys.argv[1] if len(sys.argv) > 1 else "einreichung_top20.csv"
REF = sys.argv[sys.argv.index("--ref") + 1] if "--ref" in sys.argv else "6ARU.pdb"
if not os.path.exists(REF):
    for _d in (".", os.path.expanduser("~/egfr-ph-binder"), "structures", "targets"):
        if os.path.exists(os.path.join(_d, os.path.basename(REF))):
            REF = os.path.join(_d, os.path.basename(REF)); break
RECEPTOR = sys.argv[sys.argv.index("--receptor") + 1] if "--receptor" in sys.argv else "A"
ANTIBODY = (sys.argv[sys.argv.index("--antibody") + 1].split(",")
            if "--antibody" in sys.argv else ["B", "C"])
CUT = 4.5

# Zielausschnitt je Lauf. Gesucht wird in mehreren Verzeichnissen, damit das
# Skript sowohl aus dem Arbeitsverzeichnis als auch aus dem Repository laeuft.
PATCHNAMEN = {"B": "EGFR_d3_trim.pdb", "B2": "EGFR_d3_trim.pdb",
              "C": "EGFR_6aru_patchC.pdb", "F": "EGFR_6aru_patchF.pdb"}
SUCHPFADE = ["targets", ".", "structures",
             os.path.expanduser("~/egfr-ph-binder"),
             os.path.expanduser("~/egfr-ph-binder/structures")]
if "--patchdir" in sys.argv:
    SUCHPFADE.insert(0, sys.argv[sys.argv.index("--patchdir") + 1])

def patchpfad(lauf):
    n = PATCHNAMEN.get(lauf)
    if not n:
        return None
    for d in SUCHPFADE:
        p = os.path.join(d, n)
        if os.path.exists(p):
            return p
    return None

def rd(p, atom_only=True):
    r = {}
    for ln in open(p):
        if ln.startswith("ATOM") or (not atom_only and ln.startswith("HETATM")):
            if ln[16] not in (" ", "A"):
                continue
            r.setdefault((ln[21], int(ln[22:26])), []).append(
                (ln[12:16].strip(),
                 (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))))
    return r

def ch(r, c):
    return {k[1]: v for k, v in r.items() if k[0] == c}

def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc

ref = rd(REF)
rec = ch(ref, RECEPTOR)
if not rec:
    sys.exit(f"Kette {RECEPTOR} in {REF} nicht gefunden")

# --- Cetuximab-Epitop: Rezeptorreste nahe an den Antikoerperketten
ab = np.array([x for c in ANTIBODY for v in ch(ref, c).values() for _, x in v])
if ab.size == 0:
    sys.exit(f"Antikoerperketten {ANTIBODY} in {REF} nicht gefunden")
cetux = []
for i, at in sorted(rec.items()):
    a = np.array([x for _, x in at])
    if np.linalg.norm(a[:, None, :] - ab[None, :, :], axis=-1).min() <= CUT:
        cetux.append(i)
cetux = set(cetux)
print(f"Cetuximab-Epitop in {REF} (Kette {RECEPTOR}, Ketten {','.join(ANTIBODY)}, "
      f"{CUT} A): {len(cetux)} Reste")
print(f"  {sorted(cetux)}\n")

print(f"{'Design':<32} {'Epitop':>7} {'gemeinsam':>10} {'Anteil':>8}")
zeilen = []
for _, r in pd.read_csv(CSV).iterrows():
    p = r.get("PDB")
    if not isinstance(p, str) or not os.path.exists(p):
        continue
    des = rd(p)
    tgt, bnd = ch(des, "A"), ch(des, "B")
    lauf = str(r["Design"]).split("_")[1]
    pf = patchpfad(lauf)
    if not pf:
        print(f"{r['Design']:<32} Zielpatch {PATCHNAMEN.get(lauf)} nicht gefunden "
              f"(gesucht in: {', '.join(SUCHPFADE)})"); continue
    praw = rd(pf)
    pres = ch(praw, sorted({k[0] for k in praw})[0])
    dids, pids = sorted(tgt), sorted(pres)
    if len(dids) != len(pids):
        print(f"{r['Design']:<32} Patch passt nicht"); continue
    P = [dict(tgt[d]).get("CA") for d in dids]
    Q = [dict(pres[q]).get("CA") for q in pids]
    pair = [(a, b) for a, b in zip(P, Q) if a and b]
    R, t = kabsch(np.array([a for a, _ in pair]), np.array([b for _, b in pair]))
    b = np.array([x for v in bnd.values() for _, x in v]) @ R.T + t
    epi = set()
    for i, at in sorted(rec.items()):
        a = np.array([x for _, x in at])
        if np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1).min() <= CUT:
            epi.add(i)
    gem = epi & cetux
    anteil = 100.0 * len(gem) / len(epi) if epi else 0.0
    print(f"{r['Design']:<32} {len(epi):7d} {len(gem):10d} {anteil:7.0f}%")
    zeilen.append(dict(Design=r["Design"], Epitop=len(epi),
                       GemeinsamMitCetuximab=len(gem), Anteil=round(anteil, 1),
                       Reste=" ".join(map(str, sorted(epi)))))
if zeilen:
    ziel = "data/epitope_overlap.csv" if os.path.isdir("data") else "epitope_overlap.csv"
    pd.DataFrame(zeilen).to_csv(ziel, index=False)
    print(f"\ngeschrieben: {ziel}")
print("\nAnteil = Anteil des Design-Epitops, der auch von Cetuximab beruehrt "
      f"wird (Schwelle {CUT} A).")
