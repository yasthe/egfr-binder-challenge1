#!/usr/bin/env python3
"""Abschlusspruefung der Einreichungskandidaten.

Prueft, was in den bisherigen Durchlaeufen NICHT geprueft wurde. Reine
CPU-Arbeit, keine Strukturvorhersage.

  1  Sequenz-Schwachstellen   Ladung, pI, Hydrophobizitaet, Cystein, Methionin,
                              Deamidierungsmotive, hydrophobe Strecken
  2  Aehnlichkeit untereinander   paarweise Identitaet - zwei fast gleiche
                              Designs belegen zwei Plaetze fuer ein Molekuel
  3  Grenzflaeche unabhaengig nachgerechnet   Rosetta InterfaceAnalyzer statt
                              der BindCraft-Zahlen aus dem Entwurf
  4  Binder allein            Energie und Rest-Hydrophobizitaet im ungebundenen
                              Zustand - Aggregationsrisiko
  5  Epitop                   welche EGFR-Reste beruehrt werden, Treffer der
                              vorgesehenen Hotspots
  6  Glykosylierung           liegt eine N-Glykosylierungsstelle im Epitop?
                              Der echte Rezeptor traegt dort Zucker.
  7  Sterik in mehreren Rezeptorzustaenden   dieselbe Pruefung gegen jede
                              vorhandene Referenzstruktur, nicht nur eine

Aufruf:
  python final_checks.py einreichung_top20.csv [--refs 6ARU.pdb,structures/1IVO.pdb,structures/1YY9.pdb]
  Zusatz --relax  rechnet die Grenzflaeche nach FastRelax nach (deutlich langsamer)
"""
import os, re, sys, glob, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")

CSV   = sys.argv[1] if len(sys.argv) > 1 else "einreichung_top20.csv"
REFS  = (sys.argv[sys.argv.index("--refs") + 1].split(",")
         if "--refs" in sys.argv else
         ["6ARU.pdb", "structures/1IVO.pdb", "structures/1YY9.pdb"])
RELAX = "--relax" in sys.argv
CLASH_CUT, NEAR_CUT, GLYC_CUT = 2.5, 5.0, 10.0

AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
       "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
       "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
KD = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,
      "H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,
      "T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
PKA = {"D":3.65,"E":4.25,"C":8.18,"Y":10.07,"H":6.00,"K":10.53,"R":12.48}
NEG, POS = "DECY", "HKR"

def pI(seq):
    def q(ph):
        n = 1/(1+10**(ph-9.69)) - 1/(1+10**(2.34-ph))
        for a in seq:
            if a in NEG: n -= 1/(1+10**(PKA[a]-ph))
            elif a in POS: n += 1/(1+10**(ph-PKA[a]))
        return n
    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo+hi)/2
        if q(mid) > 0: lo = mid
        else: hi = mid
    return (lo+hi)/2

def ladung(seq, ph=7.4):
    n = 1/(1+10**(ph-9.69)) - 1/(1+10**(2.34-ph))
    for a in seq:
        if a in NEG: n -= 1/(1+10**(PKA[a]-ph))
        elif a in POS: n += 1/(1+10**(ph-PKA[a]))
    return n

def ident(a, b):
    if len(a) == len(b):
        return 100.0*sum(x == y for x, y in zip(a, b))/len(a)
    try:
        from Bio import Align
        al = Align.PairwiseAligner(mode="global", open_gap_score=-11,
                                   extend_gap_score=-1)
        al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
        aln = al.align(a, b)[0]
        m = sum(sum(a[s1+k] == b[s2+k] for k in range(e1-s1))
                for (s1, e1), (s2, e2) in zip(aln.aligned[0], aln.aligned[1]))
        return 100.0*m/min(len(a), len(b))
    except Exception:
        return float("nan")

def lies_pdb(path, nur_atom=True):
    """-> {(kette, resnr): {'name':resname, 'atome':[(atomname, xyz)]}}"""
    res = {}
    for ln in open(path):
        if nur_atom and not ln.startswith("ATOM"):
            continue
        if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
            continue
        if ln[16] not in (" ", "A"):          # nur erste Konformation
            continue
        key = (ln[21], int(ln[22:26]))
        xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        r = res.setdefault(key, {"name": ln[17:20].strip(), "atome": []})
        r["atome"].append((ln[12:16].strip(), xyz))
    return res

def kette(res, ch):
    return {k[1]: v for k, v in res.items() if k[0] == ch}

def seq_von(resdict):
    return "".join(AA3.get(resdict[i]["name"], "X") for i in sorted(resdict))

def koords(resdict, nur=None):
    out = []
    for i in sorted(resdict):
        for an, xyz in resdict[i]["atome"]:
            if nur is None or an in nur:
                out.append((i, an, xyz))
    return out

def kabsch(P, Q):
    """Rotation+Translation, die P auf Q legt."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P-pc).T @ (Q-qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc

def sequons(seq, start):
    """N-X-S/T, X != P. -> Liste der Asn-Resnummern."""
    out = []
    for i in range(len(seq)-2):
        if seq[i] == "N" and seq[i+1] != "P" and seq[i+2] in "ST":
            out.append(start + i)
    return out

# ------------------------------------------------------------------ Eingabe
df = pd.read_csv(CSV)
if "PDB" not in df.columns:
    sys.exit("Spalte PDB fehlt - bitte die von rank_final.py erzeugte CSV nutzen.")
df = df[df["PDB"].notna()].reset_index(drop=True)
print(f"{len(df)} Kandidaten aus {CSV}\n")

hotspots_map = {}
for c in sorted(set(glob.glob("designs*/final_design_stats.csv") +
                    glob.glob("designs*/**/final_design_stats.csv", recursive=True))):
    try:
        t = pd.read_csv(c)
        for _, q in t.iterrows():
            h = q.get("Target_Hotspot")
            if isinstance(h, str):
                hotspots_map[str(q["Design"])] = [
                    int(x) for x in re.findall(r"[A-Za-z]?(\d+)", h)]
    except Exception:
        pass
print(f"Hotspot-Angaben fuer {len(hotspots_map)} Designs gefunden")

refs = {}
for p in REFS:
    if os.path.exists(p):
        refs[os.path.basename(p).replace(".pdb", "")] = lies_pdb(p, nur_atom=False)
        print(f"Referenz geladen: {p}")
    else:
        print(f"Referenz fehlt:   {p}  (uebersprungen)")
print()

# --------------------------------------------------------- 1 Sequenzprofil
print("=" * 100)
print("1  SEQUENZ-SCHWACHSTELLEN\n")
zeilen = []
for _, r in df.iterrows():
    s = r["Sequence"]
    hyd = [a for a in s if a in "AVILMFWC"]
    laeufe = max((len(m.group()) for m in re.finditer(r"[AVILMFWY]+", s)), default=0)
    zeilen.append(dict(
        Design=r["Design"], Len=len(s), pI=round(pI(s), 2),
        Ladung74=round(ladung(s), 1), GRAVY=round(sum(KD[a] for a in s)/len(s), 2),
        Hydrophob=f"{100*len(hyd)/len(s):.0f}%", Cys=s.count("C"), Met=s.count("M"),
        NGlyc=len(sequons(s, 1)), Deamid=len(re.findall(r"N[GS]", s)),
        MaxHydLauf=laeufe))
seq_tab = pd.DataFrame(zeilen)
print(seq_tab.to_string(index=False))
print("\n  Cys>0 = Risiko falscher Disulfide.  Deamid = NG/NS, Alterungsmotive.")
print("  MaxHydLauf > 5 = zusammenhaengender hydrophober Abschnitt, Aggregationsrisiko.")
print("  NGlyc im Binder ist bei zellfreier Herstellung folgenlos, bei Zellen nicht.\n")

# ------------------------------------------------ 2 Aehnlichkeit unter sich
print("=" * 100)
print("2  AEHNLICHKEIT DER KANDIDATEN UNTEREINANDER\n")
n = len(df)
M = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        M[i, j] = 100.0 if i == j else ident(df.loc[i, "Sequence"], df.loc[j, "Sequence"])
kurz = [d.replace("EGFRpH_", "") for d in df["Design"]]
print(pd.DataFrame(M, index=kurz, columns=[k[:14] for k in kurz]).to_string(
    float_format=lambda v: f"{v:.0f}"))
paare = [(kurz[i], kurz[j], M[i, j]) for i in range(n) for j in range(i+1, n)
         if M[i, j] >= 70]
if paare:
    print("\n  Sehr aehnliche Paare (>=70% Identitaet):")
    for a, b, v in paare:
        print(f"    {a}  /  {b}   {v:.0f}%")
    print("  Diese belegen zwei Plaetze, verhalten sich experimentell aber"
          " vermutlich gleich.")
else:
    print("\n  Keine zwei Kandidaten ueber 70% Identitaet - gute Streuung.")
print()

# ------------------------------- 3/4/5 Rosetta, Binder allein, Epitop
print("=" * 100)
print("3  GRENZFLAECHE UNABHAENGIG NACHGERECHNET   4  BINDER ALLEIN   5  EPITOP\n")
rosetta = {}
try:
    from pyrosetta import init, pose_from_pdb, get_fa_scorefxn, Pose
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
    init("-mute all -ex1 -ex2aro")
    sfx = get_fa_scorefxn()
    if RELAX:
        from pyrosetta.rosetta.protocols.relax import FastRelax
        fr = FastRelax(); fr.set_scorefxn(sfx); fr.max_iter(100)
    print(f"{'Design':<30} {'dG':>8} {'dSASA':>8} {'SC':>6} {'unsat':>6} "
          f"{'nRes':>5} {'Binder':>9} {'expHyd':>7}")
    for _, r in df.iterrows():
        p = pose_from_pdb(r["PDB"])
        if RELAX:
            fr.apply(p)
        ia = InterfaceAnalyzerMover(1)          # Sprung 1 = Grenzflaeche A|B
        ia.set_scorefunction(sfx)
        ia.set_pack_separated(True)
        ia.set_compute_interface_sc(True)
        ia.set_compute_packstat(True)
        ia.apply(p)
        d = dict(dG=ia.get_interface_dG(), dSASA=ia.get_interface_delta_sasa(),
                 nres=ia.get_num_interface_residues())
        try:    d["sc"] = ia.get_all_data().sc_value
        except Exception: d["sc"] = float("nan")
        try:    d["unsat"] = ia.get_interface_delta_hbond_unsat()
        except Exception: d["unsat"] = float("nan")
        # Binder allein
        info = p.pdb_info()
        bidx = [i for i in range(1, p.total_residue()+1) if info.chain(i) == "B"]
        binder = Pose(); binder.assign(p)
        binder = p.split_by_chain()[2] if p.num_chains() > 1 else None
        d["binder_score"] = sfx(binder) if binder else float("nan")
        d["binder_len"] = binder.total_residue() if binder else 0
        # exponierte Hydrophobizitaet des freien Binders
        from pyrosetta.rosetta.core.scoring.sasa import SasaCalc
        sc_ = SasaCalc(); sc_.calculate(binder)
        sasa = sc_.get_residue_sasa()
        tot = sum(sasa[i] for i in range(1, binder.total_residue()+1))
        hyd = sum(sasa[i] for i in range(1, binder.total_residue()+1)
                  if binder.residue(i).name1() in "AVILMFWC")
        d["exp_hyd"] = 100.0*hyd/tot if tot else float("nan")
        rosetta[r["Design"]] = d
        print(f"{r['Design'][:30]:<30} {d['dG']:8.1f} {d['dSASA']:8.0f} "
              f"{d['sc']:6.2f} {d['unsat']:6.1f} {d['nres']:5d} "
              f"{d['binder_score']:9.1f} {d['exp_hyd']:6.1f}%")
    print("\n  dG/dSASA unabhaengig von den BindCraft-Zahlen gerechnet. Grosse"
          "\n  Abweichungen nach oben waeren ein Warnsignal.")
    print("  expHyd = Anteil hydrophober Oberflaeche am freien Binder; ueber"
          " etwa 35% steigt das Aggregationsrisiko.\n")
except Exception as e:
    print(f"  uebersprungen ({type(e).__name__}: {e})\n")

# ------------------------------------- 5/6/7 Epitop, Glykane, Mehrzustand
print("=" * 100)
print("6  GLYKOSYLIERUNG IM EPITOP   7  STERIK IN MEHREREN REZEPTORZUSTAENDEN\n")

# Zielpatch je Lauf - dieselben Dateien, die check_full.py benutzt.
PATCHES = {"B":  "structures/EGFR_d3_trim.pdb",
           "B2": "structures/EGFR_d3_trim.pdb",
           "C":  "EGFR_6aru_patchC.pdb",
           "F":  "EGFR_6aru_patchF.pdb"}

def lauf_von(design):
    t = str(design).split("_")
    return t[1] if len(t) > 1 else ""

patch_cache = {}
def patch_laden(pfad):
    if pfad not in patch_cache:
        patch_cache[pfad] = lies_pdb(pfad) if os.path.exists(pfad) else None
    return patch_cache[pfad]

def sequon_ids(resdict):
    """N-X-S/T ueber tatsaechlich benachbarte Reste, Luecken werden uebersprungen."""
    ids = sorted(resdict)
    s1 = {i: AA3.get(resdict[i]["name"], "X") for i in ids}
    out = []
    for k in range(len(ids) - 2):
        a, b, c = ids[k], ids[k+1], ids[k+2]
        if b == a + 1 and c == a + 2 and s1[a] == "N" and s1[b] != "P" and s1[c] in "ST":
            out.append(a)
    return out

def clash_bericht(bt, ref, aus_kette, aus_reste):
    """bt = Binderkoordinaten im Rahmen von ref."""
    pts, meta = [], []
    for (ch, rn), v in ref.items():
        if ch == aus_kette and rn in aus_reste:
            continue
        for an, xyz in v["atome"]:
            pts.append(xyz); meta.append((ch, rn))
    if not pts:
        return None
    fa = np.array(pts)
    dm = np.linalg.norm(fa[:, None, :] - bt[None, :, :], axis=-1).min(axis=1)
    nah = np.where(dm < NEAR_CUT)[0]
    proK = {}
    for k in nah:
        proK[meta[k][0]] = proK.get(meta[k][0], 0) + 1
    return dict(clash=int((dm < CLASH_CUT).sum()), near=int(len(nah)),
                mind=float(dm.min()), ketten=proK)

zusammen = []
for _, r in df.iterrows():
    des = lies_pdb(r["PDB"])
    tgt, bnd = kette(des, "A"), kette(des, "B")
    print(f"--- {r['Design']}")
    if not tgt or not bnd:
        print("    Ketten A/B nicht gefunden\n"); continue

    struct_seq, csv_seq = seq_von(bnd), str(r["Sequence"])
    if struct_seq == csv_seq:
        print(f"    Sequenzabgleich: stimmt mit der Struktur ueberein ({len(csv_seq)} aa)")
    else:
        gl = sum(a == b for a, b in zip(struct_seq, csv_seq))
        print(f"    !! SEQUENZABGLEICH FEHLGESCHLAGEN: Struktur {len(struct_seq)} aa, "
              f"CSV {len(csv_seq)} aa, {gl} Positionen gleich")
        print(f"       Struktur: {struct_seq}")
        print(f"       CSV:      {csv_seq}")

    bxyz = np.array([x for _, _, x in koords(bnd)])
    txyz = np.array([x for _, _, x in koords(tgt)])
    tids = [i for i, _, _ in koords(tgt)]
    dmin = np.linalg.norm(txyz[:, None, :] - bxyz[None, :, :], axis=-1).min(axis=1)
    epitop_des = sorted({tids[k] for k in np.where(dmin <= 4.5)[0]})

    zeile = dict(Design=r["Design"], SeqOK=(struct_seq == csv_seq),
                 Epitop=len(epitop_des))

    lauf = lauf_von(r["Design"])
    pfad = PATCHES.get(lauf)
    patch = patch_laden(pfad) if pfad else None
    if patch is None:
        print(f"    Zielpatch fuer Lauf {lauf} nicht gefunden ({pfad}) - "
              f"keine Ueberlagerung moeglich\n")
        zusammen.append(zeile); continue

    pch = sorted({k[0] for k in patch},
                 key=lambda c: -len(kette(patch, c)))[0]
    pres = kette(patch, pch)
    pids, dids = sorted(pres), sorted(tgt)
    if len(pids) != len(dids):
        print(f"    Zielpatch passt nicht ({len(pids)} vs {len(dids)} Reste)\n")
        zusammen.append(zeile); continue
    versatz = sorted({p - d for p, d in zip(pids, dids)})
    if len(versatz) != 1:
        print(f"    Nummerierung nicht gleichmaessig versetzt: {versatz[:5]}\n")
        zusammen.append(zeile); continue
    off = versatz[0]

    P, Q = [], []
    for d, p in zip(dids, pids):
        a = dict(tgt[d]["atome"]).get("CA"); b = dict(pres[p]["atome"]).get("CA")
        if a and b: P.append(a); Q.append(b)
    R, t = kabsch(np.array(P), np.array(Q))
    rmsd = np.sqrt(((np.array(Q) - (np.array(P) @ R.T + t))**2).sum(1).mean())
    print(f"    Zielpatch {os.path.basename(pfad)}: Versatz {off:+d}, "
          f"Ueberlagerung RMSD {rmsd:.2f} A ueber {len(P)} CA")
    if rmsd > 2.0:
        print("    !! Ueberlagerung misslungen - Ergebnisse unten nicht verwertbar\n")
        zusammen.append(zeile); continue

    epitop = [e + off for e in epitop_des]
    hot = hotspots_map.get(str(r["Design"]), [])
    treffer = [h for h in hot if h in epitop]
    print(f"    Epitop in Rezeptornummerierung ({len(epitop)} Reste): "
          f"{', '.join(str(x) for x in epitop[:18])}"
          f"{' ...' if len(epitop) > 18 else ''}")
    if hot:
        print(f"    Hotspots {'/'.join(map(str, hot))}: "
              f"{', '.join(map(str, treffer)) if treffer else 'keiner getroffen'}")
    zeile["Hotspots"] = f"{len(treffer)}/{len(hot)}" if hot else "-"

    b6 = bxyz @ R.T + t                      # Binder im Rahmen von 6ARU
    aus = set(pids)

    for name, ref in refs.items():
        if name == "6ARU":
            bt, ketten_aus = b6, pch
        else:
            q6 = refs.get("6ARU")
            if q6 is None:
                continue
            beste, bestn = None, 0
            for ch in sorted({k[0] for k in ref}):
                g = len(set(kette(ref, ch)) & set(kette(q6, pch)))
                if g > bestn: beste, bestn = ch, g
            if bestn < 50:
                print(f"    {name}: keine passende Rezeptorkette"); continue
            rc, qc = kette(ref, beste), kette(q6, pch)
            A_, B_ = [], []
            for i in sorted(set(rc) & set(qc)):
                x = dict(qc[i]["atome"]).get("CA"); y = dict(rc[i]["atome"]).get("CA")
                if x and y: A_.append(x); B_.append(y)
            R2, t2 = kabsch(np.array(A_), np.array(B_))
            rm2 = np.sqrt(((np.array(B_) - (np.array(A_) @ R2.T + t2))**2).sum(1).mean())
            if rm2 > 4.0:
                print(f"    {name}: Rezeptor-Ueberlagerung {rm2:.1f} A - "
                      f"andere Konformation, nicht verwertbar"); continue
            bt, ketten_aus = b6 @ R2.T + t2, beste
            print(f"    {name}: Rezeptor ueberlagert, RMSD {rm2:.2f} A ueber {len(A_)} CA")

        rep = clash_bericht(bt, ref, ketten_aus, aus)
        if rep is None:
            continue
        kette_txt = (", ".join(f"{c}:{n}" for c, n in sorted(rep["ketten"].items()))
                     if rep["ketten"] else "frei")
        print(f"       {name}: Clashes {rep['clash']}, Atome<5A {rep['near']}, "
              f"min {rep['mind']:.2f} A   [{kette_txt}]")
        zeile[f"{name}_clash"] = rep["clash"]
        zeile[f"{name}_min"] = round(rep["mind"], 2)

        if name == "6ARU":
            nah = []
            for asn in sequon_ids(pres if False else kette(ref, ketten_aus)):
                ax = np.array([x for _, x in kette(ref, ketten_aus)[asn]["atome"]])
                dd = float(np.linalg.norm(ax[:, None, :] - bt[None, :, :],
                                          axis=-1).min())
                if dd <= GLYC_CUT:
                    nah.append((asn, round(dd, 1)))
            zeile["Glyc"] = len(nah)
            if nah:
                print("       Glykosylierungsstellen innerhalb "
                      f"{GLYC_CUT:.0f} A: "
                      + ", ".join(f"N{a} ({d} A)" for a, d in nah))
    zusammen.append(zeile)
    print()

if zusammen:
    z = pd.DataFrame(zusammen)
    print("=" * 100)
    print("UEBERSICHT\n")
    print(z.to_string(index=False))
    z.to_csv("abschlusspruefung.csv", index=False)
    print("\ngeschrieben: abschlusspruefung.csv")
    print("\nIn 1IVO sind die Ketten C und D das gebundene EGF. Treffer dort sind")
    print("kein Ausschlussgrund - Cetuximab bindet dieselbe Region.")

print("""
NICHT GEPRUEFT, und von Hand nachzuholen:
  - Aehnlichkeit zu natuerlichen Proteinen (BLAST gegen nr). Ein Design, das
    einem bekannten Protein stark aehnelt, wirkt in einer Einreichung schlecht.
  - Faltung der Sequenz ohne das Ziel (braucht AlphaFold2, also GPU).
  - Bindung an HER2/HER3/HER4 (ebenfalls Strukturvorhersage).
""")
