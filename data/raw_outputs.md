# Raw tool outputs (verbatim, for verification)

**Run naming.** The write-up numbers the runs A–D; the directories are named
after the target patch. The mapping is: run A = `designs_B/`,
run B = `designs_B2/`, run C = `designs_C/`, run D = `designs_F/`.

**German terms appearing in tool output.** `sauber` = passes the strict
criterion (no clashes, no atom within 5 A); `randnah` = borderline, at the edge
of the distance cut; `Abweichend%` = percentage of contacted EGFR residues that
differ in HER2, HER3 and HER4; `Pruefung bestanden` (in the script output) = checks passed.

**Chain labels in the contact tables.** The chain composition below counts CA
atoms, i.e. polymer chains only. The contact tables were computed over all ATOM
and HETATM records, so chain identifiers that hold only heteroatoms (water,
glycans, ions) also appear there and are absent from the composition list: 6ARU
chain D and 1YY9 chain B are of that kind. Contacts attributed to them are
contacts with solvent or carbohydrate, not with protein, and carry no weight.

## BindCraft accepted designs, per run

designs_B (3 accepted), TargetSettings EGFR_patchB, default_filters, helicity -0.3, hotspot A431-436:
EGFRpH_B_l62_s216692_mpnn8   len 62 seed 216692 i_pTM 0.82 dG -50.07 dSASA 1747.58
EGFRpH_B_l69_s170491_mpnn4   len 69 seed 170491 i_pTM 0.76 dG -48.75 dSASA 2088.79
EGFRpH_B_l69_s170491_mpnn8   len 69 seed 170491 i_pTM 0.73 dG -46.07 dSASA 2019.57

designs_B2 (18 accepted), EGFR_patchB, relaxed_filters, hotspots A431,A434,A436:
l47_s84372_mpnn12  0.85 -64.06 2256.30 | l47_s84372_mpnn17  0.85 -64.01 2250.99
l51_s210945_mpnn13 0.74 -51.48 1914.91 | l52_s590370_mpnn4  0.86 -67.09 2373.95
l52_s590370_mpnn9  0.86 -64.80 2329.73 | l53_s457194_mpnn2  0.84 -67.72 2120.46
l53_s457194_mpnn3  0.85 -68.22 2096.48 | l55_s930282_mpnn8  0.82 -57.08 1940.17
l55_s930282_mpnn15 0.83 -51.90 1896.13 | l58_s629057_mpnn5  0.80 -55.67 1819.15
l58_s629057_mpnn8  0.76 -58.83 1789.50 | l49_s595942_mpnn2  0.82 -62.66 1963.20
l49_s595942_mpnn3  0.82 -53.99 1931.19 | l53_s725081_mpnn9  0.64 -56.23 1969.07
l53_s725081_mpnn13 0.66 -61.03 2113.30 | l48_s811300_mpnn13 0.75 -57.05 2123.99
l49_s845828_mpnn1  0.80 -45.25 1604.58 | l49_s845828_mpnn2  0.78 -45.38 1592.71

designs_C (21 accepted), EGFR_patchC, relaxed_filters, hotspots A323,A344,A355:
(columns: design, length, seed, i_pTM, dG, dSASA)
l48_s238691_mpnn1   48 238691 0.74 -54.14 2042.08
l48_s238691_mpnn2   48 238691 0.78 -52.82 2065.02
l55_s661569_mpnn4   55 661569 0.80 -41.70 1799.35
l55_s661569_mpnn6   55 661569 0.80 -39.17 1893.21
l55_s205293_mpnn13  55 205293 0.72 -72.13 2606.81
l56_s447873_mpnn17  56 447873 0.65 -51.71 1723.51
l54_s137423_mpnn1   54 137423 0.59 -40.14 1495.02
l54_s137423_mpnn3   54 137423 0.65 -38.30 1575.93
l49_s11396_mpnn2    49  11396 0.85 -56.78 1879.16
l49_s11396_mpnn4    49  11396 0.85 -54.50 2024.99
l46_s573572_mpnn1   46 573572 0.82 -63.86 2029.61
l46_s573572_mpnn2   46 573572 0.83 -59.97 2033.47
l58_s786788_mpnn9   58 786788 0.80 -44.25 1561.68
l58_s786788_mpnn10  58 786788 0.81 -46.23 1583.33
l45_s295112_mpnn1   45 295112 0.66 -38.81 1295.11
l45_s295112_mpnn2   45 295112 0.65 -33.61 1256.37
l58_s400755_mpnn7   58 400755 0.82 -60.19 1813.24
l46_s483187_mpnn1   46 483187 0.79 -50.48 1955.62
l46_s483187_mpnn2   46 483187 0.78 -50.92 1903.96
l54_s965837_mpnn3   54 965837 0.78 -50.81 2177.47
l54_s965837_mpnn10  54 965837 0.77 -58.07 2296.53

## Per-design values for the seven submitted designs (rank_final.py output)

Design                        ddG_pH  paralogDiff%  i_pTM  i_pAE  BindCraft_dG  BindCraft_dSASA  Unsat  Hotspot_RMSD  clashes  minDist
EGFRpH_B2_l53_s725081_mpnn13   -0.45     48.48      0.66   0.28     -61.03         2113.30       3.50     1.58         0       12.73
EGFRpH_C_l55_s661569_mpnn6     +0.08     54.55      0.80   0.18     -39.17         1893.21       5.00     1.44         2        3.31
EGFRpH_C_l45_s295112_mpnn1     +0.02     58.73      0.66   0.27     -38.81         1295.11       1.50     2.33         0       12.94
EGFRpH_C_l54_s137423_mpnn3     +0.13     64.29      0.65   0.30     -38.30         1575.93       3.00     1.41         0        3.69
EGFRpH_C_l45_s295112_mpnn2     +0.02     58.33      0.65   0.28     -33.61         1256.37       1.50     2.54         0       13.31
EGFRpH_B2_l53_s725081_mpnn9    +0.36     46.88      0.64   0.30     -56.23         1969.07       2.00     1.28         0       12.53
EGFRpH_C_l55_s661569_mpnn4     +0.07     53.09      0.80   0.18     -41.70         1799.35       4.50     1.41         1        3.57

(Unsat and dG in these columns are BindCraft's own figures; the recomputed
Rosetta values are in the InterfaceAnalyzer table below. That recomputation
uses the same mover and score function BindCraft uses internally, on an
unrelaxed pose - see METHODS.md section 6. The Unsat columns are not
comparable: BindCraft uses a BuriedUnsatHbonds filter, the recomputation uses
get_interface_delta_hbond_unsat.)

designs_F (5 accepted), EGFR_patchF, relaxed_filters, hotspots A431,A434,A436:
l37_s14679_mpnn2   37  0.77 -58.78 1657.23
l37_s14679_mpnn11  37  0.83 -60.74 1627.07
l32_s730522_mpnn2  32  0.89 -50.09 1420.13
l32_s730522_mpnn3  32  0.89 -51.21 1432.08
l34_s472129_mpnn6  34  0.89 -63.27 1673.44

## Target patches

EGFR_6aru_patchC.pdb       chain A, residues 300-450, 151 residues
EGFR_6aru_patchF.pdb       chain A, residues 365-520, 156 residues
structures/EGFR_d3_trim.pdb chain A, residues 310-455, 146 residues
structures/EGFR_domainIII.pdb chain A, residues 310-480, 171 residues (not used in a run)

## Reference structures, chain composition (CA counts)

6ARU   A 609 (EGFR), B 210 (cetuximab light), C 212 (cetuximab heavy)
       HEADER TRANSFERASE/IMMUNE SYSTEM, 23-AUG-17
       TITLE  STRUCTURE OF CETUXIMAB FAB MUTANT IN COMPLEX WITH EGFR EXTRACELLULAR DOMAIN
1YY9   A 613 (EGFR), C 211 (cetuximab Fab light), D 220 (cetuximab Fab heavy)
1IVO   A 511, B 510 (EGFR), C 47, D 47 (EGF)

## check_full.py, corrected per-run invocation

clash_B.txt  : sauber 0 / 3
clash_B2.txt : sauber 2 / 18
clash_C.txt  : sauber 2 / 21
clash_F.txt  : sauber 0 / 5

(earlier, incorrect single invocation with the run B patch for all 47: "sauber 7 / 47")

Least-clashing designs, sorted by clashes then atoms<5A
(columns: clashes, atoms<5A, min distance, RMSD, file):
       0         0    12.53   0.71  EGFRpH_B2_l53_s725081_mpnn9_model1.pdb  OK
       0         0    12.73   0.73  EGFRpH_B2_l53_s725081_mpnn13_model2.pdb  OK
       0         0    12.94   0.96  EGFRpH_C_l45_s295112_mpnn1_model2.pdb  OK
       0         0    13.31   1.37  EGFRpH_C_l45_s295112_mpnn2_model2.pdb  OK
       0        13     3.69   0.92  EGFRpH_C_l54_s137423_mpnn3_model1.pdb  randnah
       1         8     3.57   1.14  EGFRpH_C_l55_s661569_mpnn4_model2.pdb
       2         7     3.31   1.07  EGFRpH_C_l55_s661569_mpnn6_model2.pdb
       8        12     1.94   1.01  EGFRpH_C_l56_s447873_mpnn17_model1.pdb
      10        16     2.65   0.96  EGFRpH_C_l48_s238691_mpnn2_model2.pdb
      15        28     2.35   1.39  EGFRpH_C_l54_s137423_mpnn1_model1.pdb
      20        21     1.56   0.93  EGFRpH_C_l58_s786788_mpnn9_model2.pdb
      42        32     1.08   0.98  EGFRpH_C_l58_s786788_mpnn10_model2.pdb
      60        49     2.02   0.87  EGFRpH_B_l62_s216692_mpnn8_model1.pdb
    3224       595     0.07  11.88  EGFRpH_C_l46_s573572_mpnn2_model2.pdb   (worst case)

## ph_scan.py over the accepted set (36 structures)

Best (most negative): -0.45 kcal/mol, EGFRpH_B2_l53_s725081_mpnn13
Worst: +0.519 kcal/mol, EGFRpH_C_l58_s400755_mpnn7
Submitted designs: -0.45, +0.02, +0.02, +0.07, +0.08, +0.13, +0.36

Histidine-bias MPNN variants (28 threaded variants), full ranking range:
+0.132 (b00_03) to +0.505 (b15_07) kcal/mol. All positive. Controls at bias 0.0
gave +0.132 to +0.168.

Yield of the 1-3 histidine / no-cysteine filter, per bias level (out of 40
sequences each): bias 0.0 -> 4 passed; bias 0.5 -> 20; bias 1.0 -> 36;
bias 1.5 -> 30. The best eight per level were carried forward, except at bias
0.0 where only four existed. Hence 4 + 8 + 8 + 8 = 28 threaded variants.

Example per-histidine values (b15_07): His14 pKa free 6.41 -> bound 4.05;
His18 6.10 -> 2.51; His47 5.47 -> 4.29.

## final_checks.py on the seven candidates

Sequence liabilities:
Design                        Len   pI   charge7.4  GRAVY  hydrophob Cys Met NGlyc Deamid maxHydRun
EGFRpH_B2_l53_s725081_mpnn13   53  4.65   -5.9      -1.37    28%      0   1    0      0      2
EGFRpH_C_l55_s661569_mpnn6     55  5.43   -2.9      -0.67    35%      0   2    0      1      4
EGFRpH_C_l45_s295112_mpnn1     45  5.12   -2.0      -0.75    33%      0   1    1      0      3
EGFRpH_C_l54_s137423_mpnn3     54  4.88   -3.9      -1.34    26%      0   2    0      0      4
EGFRpH_C_l45_s295112_mpnn2     45  7.65    0.0      -0.96    31%      0   1    1      0      3
EGFRpH_B2_l53_s725081_mpnn9    53  4.96   -3.9      -1.39    30%      0   1    0      0      2
EGFRpH_C_l55_s661569_mpnn4     55  4.87   -4.9      -0.66    35%      0   2    0      1      4

Pairwise identity (%): mpnn13/mpnn9 85, 661569_mpnn6/mpnn4 89, 295112_mpnn1/mpnn2 80.
All other pairs 15-32.

Rosetta InterfaceAnalyzer (recomputation, same mover as BindCraft, unrelaxed pose):
Design                          dG     dSASA   SC   unsat  nRes  BinderScore expHyd
EGFRpH_B2_l53_s725081_mpnn13  -56.4    2139   0.58   6.0    61    -121.9     22.5%
EGFRpH_C_l55_s661569_mpnn6    -30.9    1896   0.56  15.0    58    -150.1     14.1%
EGFRpH_C_l45_s295112_mpnn1    -32.1    1281   0.73   4.0    50    -107.2     18.8%
EGFRpH_C_l54_s137423_mpnn3    -34.5    1496   0.60   9.0    56    -134.8     14.8%
EGFRpH_C_l45_s295112_mpnn2    -27.1    1256   0.63   4.0    48    -115.3     18.4%
EGFRpH_B2_l53_s725081_mpnn9   -52.6    1945   0.63   8.0    66    -123.4     21.7%
EGFRpH_C_l55_s661569_mpnn4    -36.2    1775   0.57  18.0    59    -142.6     14.7%

Epitopes (mature EGFR numbering, contacts within 4.5 A):
mpnn13:          316,317,318,321,322,323,324,325,326,343,344,346,348,349,350,353,354,355,356,357,358,379,380,382,384,405,406,407,408,409,410
661569_mpnn6:    323,324,325,346,348,349,350,353,355,356,357,358,380,382,384,408,409,411,412,415,417,418,419,420,438,440,441,443
295112_mpnn1:    316,317,322,323,324,325,348,349,350,353,354,355,356,357,358,359,384
137423_mpnn3:    316,317,323,324,325,346,348,349,350,353,355,356,357,358,382,384,385,408,409,412,417,418,438,440
295112_mpnn2:    316,317,322,323,324,325,348,349,350,353,354,355,356,357,358,359,384,409
mpnn9:           316,317,318,321,322,323,324,325,326,342,344,346,348,349,350,353,354,355,356,357,358,380,382,384,405,406,407,408,409,410
661569_mpnn4:    325,346,348,349,350,353,355,356,357,358,380,382,384,408,409,411,412,415,417,418,438,440,441,442,443,444

Hotspots specified vs contacted: run B designs 0/3 (hotspots 431,434,436);
run C designs 2/3 except 661569_mpnn4 at 1/3 (hotspots 323,344,355).

Superposition and contacts in the reference structures
(atoms within 5 A, per chain; patch residues of the design's own run excluded):
mpnn13:        6ARU [C:185, D:49]  1YY9 [A:5, B:46, D:197]
661569_mpnn6:  6ARU [A:3, B:23, C:513, D:42]  1YY9 [A:10, B:35, C:23, D:530]
295112_mpnn1:  6ARU [C:232, D:44]  1YY9 [A:2, B:35, D:252]
137423_mpnn3:  6ARU [A:4, B:64, C:420, D:36]  1YY9 [A:8, B:26, C:81, D:418]
295112_mpnn2:  6ARU [C:210, D:53]  1YY9 [A:3, B:42, D:223]
mpnn9:         6ARU [C:199, D:51]  1YY9 [A:5, B:45, D:215]
661569_mpnn4:  6ARU [A:6, B:22, C:514, D:24]  1YY9 [A:10, B:17, C:25, D:537]
1IVO in all cases: receptor superposition RMSD 16.9 A, refused as not comparable.
Patch superposition RMSD: 0.69-0.73 A (run B designs), 0.92-1.37 A (run C designs).
1YY9 receptor superposition onto 6ARU: RMSD 0.68 A over 609 CA.

Glycosylation sites within 10 A of the binder (6ARU frame):
mpnn13 N328 4.1 | mpnn9 N328 4.1 | 295112_mpnn1 N328 8.2 | 295112_mpnn2 N328 8.2
137423_mpnn3 N328 6.9, N420 9.4 | 661569_mpnn6 N328 7.3, N389 9.6, N420 4.7
661569_mpnn4 N328 7.8, N420 7.3

Cross-reactivity proxy (% of contacted EGFR residues differing in HER2/3/4):
46.88 to 64.29 across the seven.

## BLAST

blastp, ClusteredNR, 3 October 2026, all seven queries:
"No significant similarity found."

## epitope_overlap.py, 3 October 2026

Cetuximab contact set in 6ARU: receptor chain A residues within 4.5 A of
chains B and C (light and heavy chain), 24 residues:
349, 350, 353, 382, 384, 408, 409, 411, 412, 415, 417, 418, 438, 440, 441,
443, 465, 466, 467, 468, 469, 471, 472, 473

Design epitopes in the same frame (receptor chain A of 6ARU, 4.5 A) and their
intersection with that set:

Design                         epitope  shared  share
EGFRpH_B2_l53_s725081_mpnn13       40       9     22%
EGFRpH_C_l55_s661569_mpnn6         27      16     59%
EGFRpH_C_l45_s295112_mpnn1         18       4     22%
EGFRpH_C_l54_s137423_mpnn3         25      13     52%
EGFRpH_C_l45_s295112_mpnn2         19       4     21%
EGFRpH_B2_l53_s725081_mpnn9        40       9     22%
EGFRpH_C_l55_s661569_mpnn4         27      16     59%

Note on the differing epitope sizes: the contact lists earlier in this file
were computed against the design's own target chain (the 146-156 residue
patch). The counts here are against the full receptor chain of 6ARU after
superposition, so residues outside the design patch can also fall within
4.5 A. For s725081_mpnn13 this is 40 residues against 31 earlier.
