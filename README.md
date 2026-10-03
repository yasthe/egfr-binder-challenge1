# pH-conditional EGFR binders — Challenge 1 submission

De novo mini-binders against domain III of the human EGFR extracellular region,
designed for the Anthropic × Adaptyv Protein Design Competition, Challenge 1
(Track 3). Seven designs are submitted.

**Author.** Yannick Heimann. Work carried out September–October 2026. This is a
personal entry to a public competition. It was not funded, does not form part of
any research project, and does not represent any institution; no institutional
or rented compute was used.

**Funding and competing interests.** No funding was received. The author
declares no competing interests.

**Headline result: the pH switch was not achieved.** Of three strategies tried,
none produced a usable difference in binding between pH 6.5 and 7.4. The best
value in the set is −0.45 kcal/mol, against roughly −1.3 kcal/mol needed for a
tenfold change, and the whole observed spread lies within the uncertainty of
the method used to compute it. The designs are therefore submitted as EGFR
binders with **no measurable pH dependence**, not as a demonstrated switch. Why
the approach fails is set out in [METHODS.md](METHODS.md) §5; in short,
interface histidine pKa values fall rather than rise on binding, and in several
cases were already below the test pH to begin with.

**Second result: 40 of 47 designs that passed the
BindCraft filters collide with parts of the receptor the design stage never
saw.** The design stage only ever sees the target patch it is given, so it
cannot notice that the rest of the receptor occupies the space the binder
needs. Under a strict criterion — zero clashes and no atom within 5 Å of
anything outside the design patch — only 4 of 47 pass. That criterion is
harsher than the physics requires; a cut-off that tolerates van der Waals
contact (≥3.0 Å minimum distance, ≤2 clashes) leaves seven designs and
discards 40.

## The seven submitted designs

| Design | aa | i_pTM | ΔG | ΔSASA | SC | unsat. | hotspots | min. dist. | glycan | ΔΔG pH |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| EGFRpH_C_l45_s295112_mpnn1 | 45 | 0.66 | −32.1 | 1281 | 0.73 | 4.0 | 2/3 | 12.9 Å | 8.2 Å | +0.02 |
| EGFRpH_C_l45_s295112_mpnn2 | 45 | 0.65 | −27.1 | 1256 | 0.63 | 4.0 | 2/3 | 13.3 Å | 8.2 Å | +0.02 |
| EGFRpH_B2_l53_s725081_mpnn13 | 53 | 0.66 | −56.4 | 2139 | 0.58 | 6.0 | 0/3 | 12.7 Å | 4.1 Å | −0.45 |
| EGFRpH_B2_l53_s725081_mpnn9 | 53 | 0.64 | −52.6 | 1945 | 0.63 | 8.0 | 0/3 | 12.5 Å | 4.1 Å | +0.36 |
| EGFRpH_C_l54_s137423_mpnn3 | 54 | 0.65 | −34.5 | 1496 | 0.60 | 9.0 | 2/3 | 3.7 Å | 6.9 Å | +0.13 |
| EGFRpH_C_l55_s661569_mpnn4 | 55 | 0.80 | −36.2 | 1775 | 0.57 | 18.0 | 1/3 | 3.6 Å | 7.3 Å | +0.07 |
| EGFRpH_C_l55_s661569_mpnn6 | 55 | 0.80 | −30.9 | 1896 | 0.56 | 15.0 | 2/3 | 3.3 Å | 4.7 Å | +0.08 |

i_pTM is BindCraft's own figure. ΔG (Rosetta energy units), ΔSASA (Å²), shape
complementarity and buried unsatisfied hydrogen bonds were recomputed with the
PyRosetta InterfaceAnalyzer. This is the same mover BindCraft uses internally,
so it is a consistency check rather than an independent assessment; for five of
seven the recomputed energy is less favourable, most plausibly because
BindCraft scores a relaxed pose and this recomputation did not. "hotspots" counts how many of
the three specified target residues are actually contacted. "min. dist." is the
closest approach to receptor regions outside the design patch: the first four
designs touch nothing, the last three sit at van der Waals contact with a
neighbouring domain. "glycan" is the distance to the nearest receptor
N-glycosylation site. ΔΔG pH is in kcal/mol, negative meaning tighter binding at
pH 6.5; the whole column lies within the uncertainty of the method.

Three pairs share a backbone (80–89 % sequence identity within a pair, 15–32 %
across), so the seven designs represent four distinct scaffolds, all directed at
the same region of domain III. Superposed into the EGFR–cetuximab complexes
6ARU and 1YY9, all seven overlap the antibody heavy chain substantially, so all
seven would compete with cetuximab for binding. They do not, however, all bind
the cetuximab epitope: of each design's contacts, the s661569 pair shares 59 %
with cetuximab's own contact set and s137423_mpnn3 shares 52 %, while the
s725081 and s295112 pairs share only 21–22 % and bind a largely different face
of domain III.

## Repository layout

```
submission/     the submitted sequences and their metrics
scripts/        the analysis scripts
targets/        the three target patches cut from 6ARU
data/           BindCraft statistics per run, all check reports, raw tool output
structures/     the seven submitted complexes as PDB
METHODS.md      full methods, results and limitations
```

`6ARU.pdb` is not redistributed here; it can be downloaded from
https://files.rcsb.org/download/6ARU.pdb into the repository root, which is
where the commands below expect it.

## Reproducing the analysis

The design runs need a GPU and a BindCraft installation; the analysis steps do
not. The commands below use the repository layout above: the three target
patches in `targets/` and `6ARU.pdb` in the repository root. In the original
working directory these files sat elsewhere, so the paths recorded in
`data/raw_outputs.md` differ. Each check was called once per design
run with that run's own target patch and numbering offset. A single invocation
across all runs with one patch yields invalid results; this occurred during the
campaign and is described in METHODS.md §4.

```bash
python scripts/check_full.py "designs_C/Accepted/*.pdb"  targets/EGFR_6aru_patchC.pdb 299 > data/clash_C.txt
python scripts/check_full.py "designs_F/Accepted/*.pdb"  targets/EGFR_6aru_patchF.pdb 364 > data/clash_F.txt
python scripts/check_full.py "designs_B2/Accepted/*.pdb" targets/EGFR_d3_trim.pdb      309 > data/clash_B2.txt
python scripts/check_full.py "designs_B/Accepted/*.pdb"  targets/EGFR_d3_trim.pdb      309 > data/clash_B.txt
cat data/clash_B.txt data/clash_B2.txt data/clash_C.txt data/clash_F.txt > data/clash_report.txt

python scripts/ph_scan.py "designs_*/Accepted/*.pdb" > data/ph_report.txt
python scripts/rank_final.py data/clash_report.txt data/ph_report.txt --maxclash 2 --mindist 3.0
python scripts/final_checks.py einreichung_top20.csv
python scripts/make_submission.py einreichung_top20.csv
python scripts/epitope_overlap.py einreichung_top20.csv
```

`rank_final.py` writes `einreichung_top20.csv`, which the three following steps
read. `make_submission.py` validates the sequences against their structures and
refuses to write the submission files if anything disagrees.

## Author and AI use

This work was carried out by Yannick Heimann. All design runs, analyses,
interpretations and decisions are the author's own, and the author takes full
responsibility for the entire content of this repository, including every line
of code and text.

An AI assistant (Anthropic Claude, used September–October 2026) provided
substantial assistance with writing and debugging the analysis scripts in
`scripts/`, with structuring and drafting this documentation, and with
exploratory analysis of the numerical results. Every script was read, run and
checked by the author against the raw outputs in `data/raw_outputs.md`, and
every scientific claim, threshold and conclusion was verified by the author.
The AI assistant is a tool and is not an author; responsibility for any error
rests with the author alone.

This statement concerns AI used as a writing and programming aid. The
deep-learning methods used as scientific instruments — BindCraft, AlphaFold2,
AlphaFold2-multimer and ProteinMPNN — are described in METHODS.md §3 and cited
in §9.

## Software and data

Tools, versions where recorded, and full citations are in METHODS.md §8 and §9.
Some version numbers were not captured at the time of the runs and are marked
as such. Structures used: PDB 6ARU, 1YY9, 1IVO. Sequence reference:
UniProt P00533.

