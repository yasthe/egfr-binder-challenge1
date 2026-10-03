#!/usr/bin/env python3
"""pH-Schalter quantifizieren.

Fuer jedes Design: PROPKA3 auf den Komplex und auf den freien Binder.
Ein Histidin traegt nur dann zum Schalter bei, wenn sein pKa beim Binden
STEIGT: dann ist es bei pH 6.5 im gebundenen Zustand protoniert, bei 7.4
aber nicht, und das Binden wird sauer beguenstigt.

Die Kopplung ist thermodynamisch (Wyman-Linkage), nicht geraten:
  dG_kopplung(pH) = -RT ln[ (1+10^(pKa_geb - pH)) / (1+10^(pKa_frei - pH)) ]
  ddG_Schalter    = dG(6.5) - dG(7.4)      negativ = bei pH 6.5 fester
"""
import glob, os, sys, math, subprocess, tempfile, shutil

RT   = 0.001987 * 298.15          # kcal/mol
PH_L, PH_H = 6.5, 7.4
DES  = sys.argv[1] if len(sys.argv) > 1 else "designs_C/Accepted/*.pdb"

def run_propka(pdb, workdir):
    dst = os.path.join(workdir, os.path.basename(pdb))
    if os.path.abspath(pdb) != os.path.abspath(dst): shutil.copy(pdb, dst)
    r = subprocess.run(["propka3", os.path.basename(dst)], cwd=workdir,
                       capture_output=True, text=True)
    pka = dst[:-4] + ".pka"
    if not os.path.exists(pka):
        print("   PROPKA fehlgeschlagen:", r.stderr.strip()[:200]); return {}
    out, insum = {}, False
    for ln in open(pka):
        if ln.startswith("SUMMARY"): insum = True; continue
        if insum and ln.startswith("-"): break
        f = ln.split()
        if insum and len(f) >= 4 and f[0] == "HIS":
            out[(f[2], int(f[1]))] = float(f[3])
    return out

def coupling(pb, pf, ph):
    return -RT * math.log((1 + 10 ** (pb - ph)) / (1 + 10 ** (pf - ph)))

def binder_only(pdb, path):
    with open(pdb) as fh, open(path, "w") as o:
        for ln in fh:
            if ln.startswith(("ATOM", "HETATM")) and ln[21] == "B":
                o.write(ln)
        o.write("END\n")

files = sorted(glob.glob(DES))
if not files:
    sys.exit(f"keine Dateien unter {DES}")
print(f"{len(files)} Designs\n")
rank = []
for f in files:
    print("=" * 70); print(os.path.basename(f))
    with tempfile.TemporaryDirectory() as wd:
        cplx = run_propka(f, wd)
        fb = os.path.join(wd, "binder.pdb"); binder_only(f, fb)
        free = run_propka(fb, wd)
    his = sorted(k for k in cplx if k[0] == "B")
    if not his:
        print("  keine Histidine im Binder"); continue
    tot = 0.0
    for key in his:
        pb = cplx[key]; pf = free.get(key)
        if pf is None:
            print(f"  His{key[1]}: im freien Binder nicht gefunden"); continue
        d  = pb - pf
        dd = coupling(pb, pf, PH_L) - coupling(pb, pf, PH_H)
        tot += dd
        mark = "  <== Schalter" if d > 0.4 and dd < -0.05 else ""
        print(f"  His{key[1]:<4d} pKa frei {pf:5.2f} -> gebunden {pb:5.2f} "
              f"(D {d:+5.2f})   ddG(6.5-7.4) {dd:+6.3f} kcal/mol{mark}")
    fold = math.exp(-tot / RT) if abs(tot) < 20 else float("inf")
    print(f"  SUMME ddG {tot:+.3f} kcal/mol  ->  {fold:.1f}x fester bei pH 6.5")
    rank.append((tot, fold, os.path.basename(f)))

print("\n" + "=" * 70); print("RANGLISTE (negativer = besserer Schalter)")
for t, fo, n in sorted(rank):
    print(f"  {t:+7.3f} kcal/mol   {fo:6.1f}x   {n}")
print("\nRichtwert: ein brauchbarer Schalter braucht mindestens ~1.3 kcal/mol,")
print("das entspricht etwa einem Faktor 10 zwischen pH 6.5 und 7.4.")
