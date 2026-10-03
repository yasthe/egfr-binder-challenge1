#!/usr/bin/env python3
"""Jedes Design gegen die vollstaendige 6ARU-Ektodomaene zurueckpruefen.
Checkliste Punkt 6: nicht nur Clashes, auch Beruehrungsflaeche.
Liest den behaltenen Bereich direkt aus der Zieldatei, funktioniert also
auch bei Luecken im Ausschnitt."""
import glob, sys, numpy as np

REF    = "6ARU.pdb"
# Aufruf:
#   python check_full.py "<glob>" [zieldatei.pdb] [offset]
# Die Zieldatei ist das Ziel, gegen das DIESE Designs entworfen wurden.
# Ohne Angabe: EGFR_6aru_patchC.pdb (Lauf C).
#   Lauf C:  python check_full.py "designs_C/Accepted/*.pdb"
#   Lauf F:  python check_full.py "designs_F/Accepted/*.pdb" EGFR_6aru_patchF.pdb
#   Lauf B2: python check_full.py "designs_B2/Accepted/*.pdb" EGFR_d3_trim.pdb 309
TARGET = "EGFR_6aru_patchC.pdb"
DES    = "designs_C/Accepted/*.pdb" if len(sys.argv) < 2 else sys.argv[1]
FIXOFF = None
for a in sys.argv[2:]:
    if a.lstrip("-").isdigit():
        FIXOFF = int(a)
    else:
        TARGET = a
VDW    = {'C':1.70,'N':1.55,'O':1.52,'S':1.80}

def rd(path, chain=None, kept=None):
    out = []
    for ln in open(path):
        if not ln.startswith("ATOM"): continue
        if chain and ln[21] != chain: continue
        if ln[16] not in " A": continue
        n = int(ln[22:26]); e = (ln[76:78].strip() or ln[12:16].strip()[0])
        if kept is not None and n not in kept: continue
        out.append((n, ln[12:16].strip(), e,
                    float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return out

keptnums = sorted({a[0] for a in rd(TARGET, "A")})
OFF_MAP  = {i+1: n for i, n in enumerate(keptnums)}
if FIXOFF is not None:
    OFF_MAP = {i: i + FIXOFF for i in range(1, 1000)}
    print(f"fester Offset {FIXOFF}: Designrest i -> reifer Rest i+{FIXOFF}")
ref      = rd(REF, "A")
refCA    = {a[0]: np.array(a[3:]) for a in ref if a[1] == "CA"}
omit     = [a for a in ref if a[0] not in set(keptnums)]
Om       = np.array([a[3:] for a in omit]); Omr = np.array([VDW.get(a[2],1.7) for a in omit])
print(f"Ziel behaelt {len(keptnums)} Reste ({keptnums[0]}-{keptnums[-1]}); "
      f"geprueft wird gegen {len(omit)} Atome ausserhalb davon.")

def kabsch(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, W = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ W))
    R = V @ np.diag([1, 1, d]) @ W
    return R, Q.mean(0) - P.mean(0) @ R

files = sorted(glob.glob(DES))
print(f"{len(files)} Designs\n")
rows = []
for f in files:
    dA = rd(f, "A"); dB = rd(f, "B")
    P, Q = [], []
    for a in dA:
        if a[1] != "CA": continue
        m = OFF_MAP.get(a[0])
        if m in refCA:
            P.append(a[3:]); Q.append(refCA[m])
    if len(P) < 30: print(f"{f}: zu wenige Anker"); continue
    P, Q = np.array(P), np.array(Q)
    R, t = kabsch(P, Q)
    rms = np.sqrt(((P @ R + t - Q) ** 2).sum(1).mean())
    B  = np.array([a[3:] for a in dB]) @ R + t
    Br = np.array([VDW.get(a[2], 1.7) for a in dB])
    D  = np.linalg.norm(B[:, None, :] - Om[None, :, :], axis=2)
    ov = (Br[:, None] + Omr[None, :]) - D
    clash = int((ov > -0.4).sum())
    near5 = int((D.min(1) < 5.0).sum())            # Binderatome nahe am weggeschnittenen Teil
    rows.append((clash, near5, D.min(), rms, f))

rows.sort()
print(f"{'Clashes':>8} {'Atome<5A':>9} {'minDist':>8} {'RMSD':>6}  Datei")
for c, n5, md, rm, f in rows:
    flag = "  OK" if c == 0 and n5 <= 5 else ("  randnah" if c == 0 else "")
    print(f"{c:8d} {n5:9d} {md:8.2f} {rm:6.2f}  {f.split('/')[-1]}{flag}")
ok = [r for r in rows if r[0] == 0 and r[1] <= 5]
print(f"\nsauber: {len(ok)} / {len(rows)}")
print("Atome<5A deutlich ueber 5 heisst: keine Clashes, aber der Binder liegt")
print("trotzdem an einem Teil an, der beim Design nicht existierte.")
