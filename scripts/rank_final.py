#!/usr/bin/env python3
"""Endauswahl: die 20 Kandidaten fuer die Einreichung.

Verarbeitet nur bereits erzeugte Ausgaben weiter - check_full.py und
ph_scan.py bleiben unveraendert und werden vorher separat aufgerufen.

Reihenfolge der Kriterien, wie festgelegt:
  0. harte Vorbedingung - kein Kontakt zu Rezeptorteilen, gegen die nicht
     entworfen wurde (aus dem check_full-Bericht)
  1. Vorauswahl der besten 100 nach den BindCraft-Kennzahlen
  2. pH-Schalter   (aus dem ph_scan-Bericht)
  3. Kreuzreaktivitaet (Ersatzmass, siehe unten)
  4. Affinitaet zu hEGFR

Kreuzreaktivitaet ist ein Ersatzmass, kein Messwert: fuer jeden Binder
werden die tatsaechlich beruehrten EGFR-Reste bestimmt und geprueft, wie
viele davon sich in HER2, HER3 und HER4 unterscheiden. Viele abweichende
Kontaktreste sprechen gegen eine Kreuzbindung. Belastbar waere nur eine
Strukturvorhersage gegen die drei Verwandten.

Aufruf:
  python rank_final.py clash_report.txt ph_report.txt [--offset 309]
         [--lax] [--maxclash N] [--mindist X]
"""
import os, re, sys, glob, json, urllib.request
import numpy as np
import pandas as pd

CLASH_REPORT = sys.argv[1]
PH_REPORT    = sys.argv[2]
OFFSET       = int(sys.argv[sys.argv.index("--offset") + 1]) if "--offset" in sys.argv else 309
LAX          = "--lax" in sys.argv
MAXCLASH     = int(sys.argv[sys.argv.index("--maxclash") + 1]) if "--maxclash" in sys.argv else None
MINDIST      = float(sys.argv[sys.argv.index("--mindist") + 1]) if "--mindist" in sys.argv else None
SIGNAL       = 24          # Signalpeptid P00533: reif i  ->  UniProt i+24
NEAR_CUT     = 5.0         # Angstrom, Kontaktdefinition
N_SHORTLIST  = 100
N_FINAL      = 20
CACHE        = "paralogs.json"

PARALOGS = {"EGFR": "P00533", "HER2": "P04626", "HER3": "P21860", "HER4": "Q15303"}

# ---------------------------------------------------------------- Berichte
def parse_clash(path):
    """Zeilen wie:  0   0   12.51   0.73  design.pdb  OK"""
    rows = {}
    pat = re.compile(r"^\s*(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+(\S+\.pdb)\s*(\S*)")
    for line in open(path, encoding="utf-8", errors="replace"):
        m = pat.match(line)
        if m:
            rows[os.path.basename(m.group(5))] = dict(
                clashes=int(m.group(1)), near=int(m.group(2)),
                mindist=float(m.group(3)), rmsd=float(m.group(4)),
                status=m.group(6))
    return rows

def parse_ph(path):
    """Zeilen der Rangliste:  +0.132 kcal/mol   0.8x   design.pdb"""
    rows = {}
    pat = re.compile(r"^\s*([+-][\d.]+)\s+kcal/mol\s+\S+\s+(\S+\.pdb)")
    for line in open(path, encoding="utf-8", errors="replace"):
        m = pat.match(line)
        if m:
            rows[os.path.basename(m.group(2))] = float(m.group(1))
    return rows

clash = parse_clash(CLASH_REPORT)
ph    = parse_ph(PH_REPORT)
print(f"Clash-Bericht: {len(clash)} Eintraege   ph_scan-Bericht: {len(ph)} Eintraege")

# ---------------------------------------------------------------- Kennzahlen
csvs = sorted(glob.glob("designs*/final_design_stats.csv")) + \
       sorted(glob.glob("designs*/**/final_design_stats.csv", recursive=True))
csvs = list(dict.fromkeys(csvs))
if not csvs:
    sys.exit("Keine final_design_stats.csv gefunden.")
frames = []
for c in csvs:
    f = pd.read_csv(c)
    f["Quelle"] = c.split(os.sep)[0]
    frames.append(f)
    print(f"  {c}: {len(f)} Zeilen")
df = pd.concat(frames, ignore_index=True)
df = df[df["Design"].notna()].drop_duplicates(subset="Design", keep="first")
print(f"{len(df)} Designs insgesamt\n")

# PDB zu jedem Design finden
def find_pdb(design, quelle):
    for pat in (f"{quelle}/Accepted/{design}_model*.pdb",
                f"{quelle}/**/{design}_model*.pdb",
                f"**/{design}_model*.pdb"):
        hit = sorted(glob.glob(pat, recursive=True))
        if hit:
            return hit[0]
    return None

df["PDB"] = [find_pdb(d, q) for d, q in zip(df["Design"], df["Quelle"])]
df["Datei"] = df["PDB"].map(lambda p: os.path.basename(p) if p else None)
print(f"{df['PDB'].notna().sum()} Designs mit zugehoeriger PDB-Datei")

# ---------------------------------------------------------------- Schritt 0
df["Clashes"] = df["Datei"].map(lambda f: clash.get(f, {}).get("clashes"))
df["Nahkontakte"] = df["Datei"].map(lambda f: clash.get(f, {}).get("near"))
df["MinDist"] = df["Datei"].map(lambda f: clash.get(f, {}).get("mindist"))
geprueft = df["Clashes"].notna()
print(f"{geprueft.sum()} davon im Clash-Bericht enthalten")

if MAXCLASH is not None or MINDIST is not None:
    ok = geprueft.copy()
    if MAXCLASH is not None:
        ok &= (df["Clashes"] <= MAXCLASH)
    if MINDIST is not None:
        ok &= (df["MinDist"] >= MINDIST)
    kriterium = (f"hoechstens {MAXCLASH} Clashes" if MAXCLASH is not None else "") + \
                (" und " if MAXCLASH is not None and MINDIST is not None else "") + \
                (f"Mindestabstand {MINDIST} A" if MINDIST is not None else "")
else:
    ok = geprueft & (df["Clashes"] == 0)
    if not LAX:
        ok &= (df["Nahkontakte"] == 0)
    kriterium = "0 Clashes" + ("" if LAX else " und 0 Nahkontakte unter 5 A")
raus = int(geprueft.sum() - ok.sum())
print(f"Schritt 0: Kriterium = {kriterium}")
print(f"           {int(ok.sum())} bestehen, {raus} fallen raus\n")
df = df[ok].copy()
if df.empty:
    sys.exit("Nichts uebrig - mit --lax erneut versuchen.")

# ---------------------------------------------------------------- Schritt 1
M = {"ipTM": "Average_i_pTM", "ipAE": "Average_i_pAE", "dG": "Average_dG",
     "dSASA": "Average_dSASA", "Unsat": "Average_n_InterfaceUnsatHbonds",
     "HotRMSD": "Average_Hotspot_RMSD", "SC": "Average_ShapeComplementarity",
     "Hyd": "Average_Surface_Hydrophobicity"}
for k, c in M.items():
    df[k] = pd.to_numeric(df.get(c), errors="coerce")

def z(s, hoch_ist_gut):
    s = s.astype(float)
    sd = s.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(0.0, index=s.index)
    v = (s - s.mean()) / sd
    return v if hoch_ist_gut else -v

df["Guete"] = (2.0 * z(df["ipTM"], True) + 1.0 * z(df["ipAE"], False) +
               1.0 * z(df["dG"], False) + 0.5 * z(df["dSASA"], True) +
               1.0 * z(df["Unsat"], False) + 0.5 * z(df["HotRMSD"], False))
short = df.sort_values("Guete", ascending=False).head(N_SHORTLIST).copy()
print(f"Schritt 1: Vorauswahl auf die besten {len(short)} nach Kennzahlen")
print(f"   i_pTM {short['ipTM'].min():.2f} bis {short['ipTM'].max():.2f}, "
      f"dG {short['dG'].min():.1f} bis {short['dG'].max():.1f}\n")

# ---------------------------------------------------------------- Schritt 2
short["ddG_pH"] = short["Datei"].map(ph)
def stufe(v):
    if pd.isna(v):      return 4, "ungeprueft"
    if v <= -1.3:       return 0, "stark"
    if v <= -0.5:       return 1, "maessig"
    if v <= -0.15:      return 2, "schwach"
    return 3, "keiner"
short[["Schalterstufe", "Schalter"]] = pd.DataFrame(
    [stufe(v) for v in short["ddG_pH"]], index=short.index)
print("Schritt 2: pH-Schalter")
for s, n in short["Schalter"].value_counts().items():
    print(f"   {s}: {n}")
print("   (gestuft statt sortiert - Unterschiede unter 0.15 kcal/mol sind Rauschen)\n")

# ---------------------------------------------------------------- Schritt 3
def hole_sequenzen():
    if os.path.exists(CACHE):
        return json.load(open(CACHE))
    seqs = {}
    for name, acc in PARALOGS.items():
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        with urllib.request.urlopen(url, timeout=30) as r:
            body = r.read().decode()
        seqs[name] = "".join(body.split("\n")[1:])
    json.dump(seqs, open(CACHE, "w"))
    return seqs

def kontaktreste(pdb, cutoff=NEAR_CUT):
    """EGFR-Reste (PDB-Nummerierung Kette A) mit Kontakt zum Binder."""
    import numpy as np
    A, B = [], []
    for line in open(pdb):
        if not line.startswith("ATOM"):
            continue
        ch, rid = line[21], int(line[22:26])
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        (A if ch == "A" else B).append((rid, xyz))
    if not A or not B:
        return []
    ar = np.array([p[1] for p in A]); br = np.array([p[1] for p in B])
    d = np.linalg.norm(ar[:, None, :] - br[None, :, :], axis=-1).min(axis=1)
    return sorted({A[i][0] for i in np.where(d <= cutoff)[0]})

spezifitaet = {}
try:
    seqs = hole_sequenzen()
    from Bio import Align
    al = Align.PairwiseAligner(mode="global", open_gap_score=-11,
                               extend_gap_score=-1, substitution_matrix=
                               Align.substitution_matrices.load("BLOSUM62"))
    karten = {}
    for name in ("HER2", "HER3", "HER4"):
        a = al.align(seqs["EGFR"], seqs[name])[0]
        m = {}
        for (s1, e1), (s2, e2) in zip(a.aligned[0], a.aligned[1]):
            for k in range(e1 - s1):
                m[s1 + k + 1] = seqs[name][s2 + k]      # 1-basiert, UniProt
        karten[name] = m
    for _, r in short.iterrows():
        if not r["PDB"]:
            continue
        kon = kontaktreste(r["PDB"])
        if not kon:
            continue
        abw, ges = 0, 0
        for rid in kon:
            up = rid + SIGNAL if rid < 1000 else rid     # PDB fuehrt reife Nr.
            wt = seqs["EGFR"][up - 1] if 0 < up <= len(seqs["EGFR"]) else None
            if wt is None:
                continue
            for name in ("HER2", "HER3", "HER4"):
                p = karten[name].get(up)
                ges += 1
                if p is None or p != wt:
                    abw += 1
        if ges:
            spezifitaet[r["Design"]] = dict(anteil=100.0 * abw / ges,
                                            n_kontakte=len(kon))
    print(f"Schritt 3: Kreuzreaktivitaet fuer {len(spezifitaet)} Designs abgeschaetzt")
    print("   Ersatzmass auf Sequenzebene - kein Ersatz fuer eine Vorhersage "
          "gegen HER2/3/4\n")
except Exception as e:
    print(f"Schritt 3: uebersprungen ({type(e).__name__}: {e})")
    print("   ohne Netz oder ohne Biopython bleibt die Spalte leer\n")

short["Abweichend%"] = short["Design"].map(
    lambda d: spezifitaet.get(d, {}).get("anteil"))
short["nKontakte"] = short["Design"].map(
    lambda d: spezifitaet.get(d, {}).get("n_kontakte"))
med = short["Abweichend%"].median()
short["Spezifitaetsstufe"] = short["Abweichend%"].map(
    lambda v: 2 if pd.isna(v) else (0 if v >= med else 1))

# ---------------------------------------------------------------- Schritt 4
short["Seed"] = short["Design"].str.extract(r"_s(\d+)_")
final = short.sort_values(
    ["Schalterstufe", "Spezifitaetsstufe", "dG", "ipTM"],
    ascending=[True, True, True, False])

# hoechstens zwei je Rueckgrat, damit ein Fehlschlag nicht alles mitnimmt
gewaehlt, proseed = [], {}
for _, r in final.iterrows():
    s = r["Seed"]
    if proseed.get(s, 0) >= 2:
        continue
    proseed[s] = proseed.get(s, 0) + 1
    gewaehlt.append(r)
    if len(gewaehlt) >= N_FINAL:
        break
if len(gewaehlt) < N_FINAL:                      # auffuellen
    drin = {r["Design"] for r in gewaehlt}
    for _, r in final.iterrows():
        if r["Design"] not in drin:
            gewaehlt.append(r)
        if len(gewaehlt) >= N_FINAL:
            break
top = pd.DataFrame(gewaehlt)

print("=" * 100)
print(f"EINREICHUNG - {len(top)} Kandidaten\n")
sp = ["Design", "Schalter", "ddG_pH", "Abweichend%", "ipTM", "ipAE",
      "dG", "dSASA", "Unsat", "HotRMSD", "Clashes", "MinDist", "Seed"]
aus = top[sp].copy()
aus.columns = ["Design", "Schalter", "ddG", "Abw%", "i_pTM", "i_pAE",
               "dG", "dSASA", "Unsat", "HotRMSD", "Clash", "MinDist", "Seed"]
print(aus.to_string(index=False, float_format=lambda v: f"{v:.2f}"))

print(f"\nRueckgrate: {top['Seed'].nunique()} verschiedene")
print(f"Laengen: {sorted(set(top['Length'].astype(int)))}")
top[["Design", "Sequence", "Length", "Seed", "ddG_pH", "Abweichend%",
     "ipTM", "ipAE", "dG", "dSASA", "Unsat", "HotRMSD", "Clashes", "MinDist",
     "PDB"]].to_csv(
    "einreichung_top20.csv", index=False)
print("\ngeschrieben: einreichung_top20.csv")
with open("einreichung_top20.fasta", "w") as fh:
    for _, r in top.iterrows():
        fh.write(f">{r['Design']}\n{r['Sequence']}\n")
print("geschrieben: einreichung_top20.fasta")
