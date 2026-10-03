#!/usr/bin/env python3
"""Erzeugt die Einreichungsdateien und prueft sie auf formale Fehler.

Liest einreichung_top20.csv (Ausgabe von rank_final.py) und schreibt:
  submission/designs.csv    - name, sequence, molecule_class  (Pflichtspalten)
  submission/designs.fasta  - dieselben Sequenzen als FASTA
  submission/metrics.csv    - optionale Kennzahlen je Design

Prueft dabei: Anzahl, erlaubte Laenge 10-250, nur die 20 Standardaminosaeuren,
eindeutige Namen, keine doppelten Sequenzen, Uebereinstimmung mit der
zugehoerigen PDB-Datei.

Aufruf:  python scripts/make_submission.py [einreichung_top20.csv]
"""
import os, sys
import pandas as pd

QUELLE = sys.argv[1] if len(sys.argv) > 1 else "einreichung_top20.csv"
OUT = "submission"
AA = set("ACDEFGHIKLMNPQRSTVWY")
AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
       "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
       "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

os.makedirs(OUT, exist_ok=True)
d = pd.read_csv(QUELLE)
print(f"{len(d)} Eintraege aus {QUELLE}\n")

fehler = []

# -- Pflichtfelder
for sp in ("Design", "Sequence"):
    if sp not in d.columns:
        sys.exit(f"Spalte {sp} fehlt in {QUELLE}")

# -- Einzelpruefungen
for i, r in d.iterrows():
    name, seq = str(r["Design"]).strip(), str(r["Sequence"]).strip().upper()
    if not (10 <= len(seq) <= 250):
        fehler.append(f"{name}: Laenge {len(seq)} ausserhalb 10-250")
    ungueltig = sorted(set(seq) - AA)
    if ungueltig:
        fehler.append(f"{name}: unzulaessige Zeichen {ungueltig}")
    if len(name) == 0:
        fehler.append(f"Zeile {i}: leerer Name")
    # Sequenz muss mit Kette B der Struktur uebereinstimmen; eine fehlende
    # PDB-Datei ist ein Fehler und kein Grund, die Pruefung zu ueberspringen.
    p = r.get("PDB")
    if not isinstance(p, str) or not p.strip():
        fehler.append(f"{name}: keine PDB-Datei angegeben, Strukturabgleich unmoeglich")
    elif not os.path.exists(p):
        fehler.append(f"{name}: PDB-Datei {p} nicht gefunden, Strukturabgleich unmoeglich")
    else:
        res, seen = [], set()
        for ln in open(p):
            if ln.startswith("ATOM") and ln[21] == "B" and ln[16] in " A":
                k = int(ln[22:26])
                if k not in seen:
                    seen.add(k); res.append((k, ln[17:20].strip()))
        sseq = "".join(AA3.get(n, "X") for _, n in sorted(res))
        if not sseq:
            fehler.append(f"{name}: keine Kette B in {os.path.basename(p)}")
        elif sseq != seq:
            fehler.append(f"{name}: Sequenz weicht von {os.path.basename(p)} ab")

if d["Design"].duplicated().any():
    fehler.append("doppelte Namen: " +
                  ", ".join(d.loc[d["Design"].duplicated(), "Design"]))
if d["Sequence"].duplicated().any():
    fehler.append("doppelte Sequenzen: " +
                  ", ".join(d.loc[d["Sequence"].duplicated(), "Design"]))
if len(d) > 20:
    fehler.append(f"{len(d)} Designs - laut Challenge-1-Spezifikation sind\n                   hoechstens 20 je Teilnehmer in Track 3 zulaessig")

if fehler:
    print("PRUEFUNG FEHLGESCHLAGEN:")
    for f in fehler:
        print("  -", f)
    sys.exit(1)
print("Pruefung bestanden: Laenge, Alphabet, Eindeutigkeit, Strukturabgleich\n")

# -- Ausgabe
sub = pd.DataFrame({"name": d["Design"].str.strip(),
                    "sequence": d["Sequence"].str.strip().str.upper(),
                    "molecule_class": "protein"})
sub.to_csv(f"{OUT}/designs.csv", index=False)
with open(f"{OUT}/designs.fasta", "w") as fh:
    for _, r in sub.iterrows():
        fh.write(f">{r['name']}\n{r['sequence']}\n")

spalten = [c for c in ["Design","Length","Seed","ipTM","ipAE","dG","dSASA",
                       "Unsat","HotRMSD","ddG_pH","Abweichend%","Clashes",
                       "MinDist","PDB"] if c in d.columns]
m = d[spalten].copy()
if "MinDist" in m.columns:
    # Grenzfall = kein Clash, aber Van-der-Waals-Kontakt zu einer Nachbardomaene
    m["steric_edge_case"] = m["MinDist"] < 5.0
m.to_csv(f"{OUT}/metrics.csv", index=False)

print(f"geschrieben: {OUT}/designs.csv ({len(sub)} Zeilen)")
print(f"geschrieben: {OUT}/designs.fasta")
print(f"geschrieben: {OUT}/metrics.csv")
print("\nInhalt von designs.csv:\n")
print(sub.to_string(index=False))
