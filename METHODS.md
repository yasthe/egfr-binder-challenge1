# Methods, results and limitations

Anthropic × Adaptyv Protein Design Competition, Challenge 1 (pH-conditional
EGFR binder), Track 3. All work was carried out on a single desktop machine
with one GPU; no rented or institutional compute was used.

## 1. Task

Design a binder to human EGFR that binds at pH 6.5 and shows no detectable
binding at pH 7.4. The assay construct is the EGFR extracellular region,
residues 25–645 in UniProt P00533-1 numbering. Domain III is the recommended
epitope and PDB 6ARU the suggested reference structure. Accepted formats are
single-chain proteins, nanobodies, scFv and Fab, 10–250 residues, and Tracks 2
and 3 may submit up to 20 designs; these figures are taken from the Challenge 1
specification on the competition page, accessed 3 October 2026.

**Numbering convention.** Residue numbers in this document refer to the mature
EGFR chain, i.e. after removal of the 1–24 signal peptide. UniProt number =
mature number + 24. BindCraft numbers the target chain of its output complexes
from 1, giving a fixed per-run offset to mature numbering: +309 for runs A and
B, +299 for run C, +364 for run D.

## 2. Target patches

Three target patches were cut from chain A of 6ARU:

| File | Residues (mature) | Count | Used in |
| --- | --- | --- | --- |
| `EGFR_d3_trim.pdb` | 310–455 | 146 | runs A, B |
| `EGFR_6aru_patchC.pdb` | 300–450 | 151 | run C |
| `EGFR_6aru_patchF.pdb` | 365–520 | 156 | run D |

Domain III ends near residue 480, so the run D patch extends into domain IV.
This turned out to matter (section 4).

## 3. Design runs

BindCraft with the `default_4stage_multimer` protocol and a helicity bias of
−0.3 throughout. BindCraft generates binder backbones by optimising through
AlphaFold2-multimer and assigns sequences with ProteinMPNN.

| Run | Directory | Patch | Hotspots | Filters | Accepted | Length | i_pTM |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A | `designs_B/` | 310–455 | A431–A436 | default | 3 | 62–69 | 0.73–0.82 |
| B | `designs_B2/` | 310–455 | A431, A434, A436 | relaxed | 18 | 47–58 | 0.64–0.86 |
| C | `designs_C/` | 300–450 | A323, A344, A355 | relaxed | 21 | 45–58 | 0.59–0.85 |
| D | `designs_F/` | 365–520 | A431, A434, A436 | relaxed | 5 | 32–37 | 0.77–0.89 |

Run A used BindCraft's default filters and yielded only three designs; the
remaining runs used the relaxed filter set; designs are therefore not directly comparable across
runs. 47 designs were accepted in total. The number of
trajectories attempted per run is recorded in each run's `trajectory_stats.csv`
and `failure_csv.csv` and is not reproduced here; acceptance *rates* are
therefore not quoted anywhere in this document.

BindCraft writes more than one AlphaFold model per accepted design. The
structure carried forward for each design is the one named in the `PDB` column
of `submission/metrics.csv`. The seven retained designs come variously from
`model1` and `model2`; no selection on model quality was performed, and the
rule by which one model rather than the other entered the check is not
recorded.

## 4. Steric check against the intact receptor

During design, BindCraft only ever sees the 146–156 residue target patch. It
therefore cannot detect that a binder occupies space taken by the rest of the
receptor. Each design was checked afterwards (`check_full.py`) against all
atoms of 6ARU outside its own target patch, counting clashes, atoms within 5 Å
and the minimum heavy-atom distance.

**Definition.** "Clash" and "atoms within 5 Å" are counted by `check_full.py`.
A clash is a heavy-atom contact below a fixed distance threshold set in that
script; *the numerical value of that threshold was not recorded in this
write-up*, which is a reproducibility gap, since the second headline result
rests on it. The two quantities are not interchangeable: among the seven
retained designs, one has zero clashes and 13 atoms within 5 Å of a
neighbouring domain.

**Result.** Under the strict criterion — zero clashes *and* no atom within 5 Å —
**4 of 47 designs pass and 43 do not.** Per run: 0/3 in run A, 2/18 in run B,
2/21 in run C, 0/5 in run D. The worst case reaches 3224 clashes at a minimum
heavy-atom distance of 0.07 Å. All five run D designs fail, which is consistent
with their patch extending into domain IV.

**Selection cut.** The strict criterion is harsher than the physics requires: a
binder that touches a neighbouring domain at van der Waals distance is not
thereby excluded from binding. Sorted by minimum distance, the designs run
13.31, 12.94, 12.73, 12.53 Å — four designs that touch nothing — then 3.69,
3.57, 3.31 Å with zero to two clashes, then 2.65 Å with ten clashes, and
downwards from there. The cut was placed at **minimum distance ≥3.0 Å and ≤2
clashes**, in the gap between 3.31 and 2.65 Å, leaving **seven designs on four
backbones** and discarding 40.

Three caveats on this threshold, none of which were resolved:

- 3.0 Å is below a normal carbon–carbon van der Waals contact (≈3.4–3.7 Å), so
  the three designs retained at 3.31–3.69 Å are at the edge of steric strain,
  not comfortably clear of it. They are flagged in `submission/metrics.csv`.
- The threshold discriminates at a finer scale than the inputs justify. The
  receptor is a single rigid crystal structure, the binders are AlphaFold
  predictions, and the superposition of each design onto its patch carries an
  RMSD of 0.69–1.37 Å.
- The criterion is not strictly comparable across runs, because each run
  exempts a different region of the receptor (146, 151 or 156 residues at
  different positions). The per-run pass rates above are therefore not directly
  comparable.

In practice the clash count does most of the work: the next design below the
cut, at 2.65 Å, already carries ten clashes and is excluded on that ground
alone.

**A methodological error worth recording.** The check was first run in a single
invocation across all 47 designs using the run B patch for all of them. Runs C
and D were designed against different patches whose residues fall outside the
exempted 310–455 range, so the script counted the very surface those designs
were built against as forbidden territory. This produced seven apparently clean
designs, including run D designs with the highest i_pTM of the whole campaign.
After correction — one invocation per run with its own patch, as the script's
own header comments prescribe — only four designs pass the strict criterion and
all run D designs fail. The apparently best candidates were an artefact of the
wrong invocation.

## 5. pH dependence

Histidine pKa values were computed with PROPKA3 in the free and bound states
and converted into a binding free-energy difference between pH 6.5 and 7.4
(`ph_scan.py`). Negative values mean tighter binding at low pH. A usable switch
was defined as ≤ −1.3 kcal/mol; RT·ln(10) is 1.36 kcal/mol at 298 K and
1.42 kcal/mol at 310 K, so this threshold corresponds to roughly a nine- to
tenfold change in affinity.

Three strategies were applied in sequence:

1. **Histidine point mutations** at interface positions of a run B design, each
   position mutated individually with local repacking.
2. **Increased histidine content**, several interface histidines at once.
3. **Surface redesign with a histidine bias** (`mpnn_his.py`): backbone and
   target held fixed, ProteinMPNN redesigning the binder sequence with the
   histidine logit raised by 0.0, 0.5, 1.0 and 1.5. 40 sequences per level at
   temperature 0.2, filtered to one to three histidines and no cysteine. Of 40
   sequences per level, 4, 20, 36 and 30 respectively passed that filter; the
   best eight per level were threaded onto the original backbone and repacked,
   except at bias 0.0 where only four existed, giving 28 variants in total.
   Level 0.0 served as a control.

The filter yields are themselves informative: without a bias, only 4 of 40
ProteinMPNN sequences pass a filter requiring one to three histidines and no
cysteine. The filter rejects on both grounds at once, so the low yield cannot be
attributed to histidine content alone.

**All three strategies failed.** The 28 redesigned variants all fall on the
wrong side, between +0.13 and +0.51 kcal/mol, including the unbiased controls.
Across the accepted structures that were scanned (36 of the 47 accepted designs;
the scan globbed the `Accepted/` directories and the remaining 11 produced no
entry in the report, which was not investigated), values range from
−0.45 to +0.52 kcal/mol. The claim that −0.45 is the best value therefore holds
for the 36 structures scanned, not demonstrably for all 47.

**Why, as far as the data show.** For the variant examined residue by residue,
interface histidine pKa values drop on binding, in one case from 6.41 to 4.05
and in another from 6.10 to 2.51. A histidine whose pKa falls on binding is
predominantly neutral in the complex and therefore binds *worse* at low pH, not
better. The desired direction requires a free pKa near 6.5–7 that *rises* on
binding, for instance a histidine facing an aspartate or glutamate carboxylate
of the target so that the protonated form is stabilised in the complex.

Two qualifications. First, this residue-level diagnosis rests on one variant;
the distribution of pKa shifts across all variants was not tabulated, so
"consistently" would overstate what was measured. Second, two of the three
histidines examined had a free pKa of 6.10 and 5.47, already below the test pH
of 6.5. For those, the switch was lost before binding was considered at all:
the design process never placed a histidine with a suitable free pKa. That is
arguably the more fundamental diagnosis, and it is what the data most directly
support. Neither BindCraft nor ProteinMPNN represents protonation states or
optimises for this quantity. A targeted approach pairing histidines with target
carboxylates was not testable in the time remaining.

## 6. Final checks on the seven candidates

Carried out with `final_checks.py`, covering what the earlier steps did not.

**Sequence against structure.** All seven submitted sequences match chain B of
their checked PDB file character for character. The same comparison is repeated
by `make_submission.py`, which refuses to write the submission files if any
sequence disagrees with its structure.

**Novelty.** blastp against ClusteredNR on 3 October 2026 returns "No
significant similarity found" for all seven queries, with the web interface's
automatic parameter adjustment for short sequences enabled. For 45–55 residue
queries this result is weak evidence: short sequences rarely reach significance
against a database of that size, so the absence of a hit rules out close
relatives of known proteins but does not establish novelty in a stronger sense.

**Sequence liabilities.** No cysteine in any design. Isoelectric point
4.65–7.65, net charge at pH 7.4 between −5.9 and 0, GRAVY −1.39 to −0.66,
hydrophobic fraction 26–35 %, longest contiguous hydrophobic stretch two to
four residues, zero to one deamidation motif (NG or NS) per sequence, one to two
methionines per design as a potential oxidation site. **Two designs carry an
N-glycosylation sequon of their own:** `EGFRpH_C_l45_s295112_mpnn1` and
`_mpnn2`. This is inconsequential for cell-free expression but would matter in
a eukaryotic host.

**Interface recomputation — a consistency check, not an independent one.**
The interface was recomputed with the PyRosetta `InterfaceAnalyzerMover`,
giving ΔG −27.1 to −56.4 REU, buried surface 1256–2139 Å², shape complementarity
0.56–0.73, and 4.0–18.0 buried unsatisfied hydrogen bonds.

This is **not an independent method**. BindCraft computes its own `dG`, `dSASA`
and shape complementarity with the same mover and the same default full-atom
score function (`functions/pyrosetta_utils.py`, `score_interface`), so the two
sets of numbers can differ only in how the pose was prepared. They did differ:
BindCraft calls `score_interface` on the relaxed pose, after `pr_relax`,
whereas the recomputation here was run on the deposited model without
relaxation. (Source read at commit `7713aa0`, 21 September 2026.) For the two run B designs the
figures broadly agree (−56.4 against −61.0); for the other five the recomputed
energy is noticeably less favourable, for example −30.9 against −39.2. **The
most likely explanation is the difference in relaxation state rather than any
disagreement between methods, and this was not resolved.** No genuinely
independent structural rescoring — a different score function or a different
package — was performed.

The unsatisfied hydrogen bond counts are not comparable at all: BindCraft
derives them from a `BuriedUnsatHbonds` XML filter with specific probe radius
and burial settings, while the recomputation used
`get_interface_delta_hbond_unsat()`. The two numbers measure different things
and the earlier draft of this document wrongly compared them. Within the
recomputation alone, the s661569 pair reaches 15 and 18 against 4 for
s295112_mpnn1, and that internal comparison does stand.

**Free binder.** Rosetta score of the isolated binder −107 to −150 REU, exposed
hydrophobic fraction 14.1–22.5 %. No aggregation signal.

**Relation to the cetuximab epitope.** 6ARU (a cetuximab Fab *mutant*) and 1YY9
(wild-type cetuximab Fab) are both EGFR–cetuximab complexes. Superposed into them, all seven designs overlap the cetuximab heavy
chain substantially: 185–514 atoms within 5 Å in 6ARU and 197–537 in 1YY9,
with the s725081 pair at the low end (185 and 197) and the s661569 pair at the
high end (514 and 537). At
the same time they contact the receptor outside their own patch with only 0–10
atoms in those superposed frames. That establishes that the designs and cetuximab
cannot bind simultaneously. Whether they contact the same receptor residues is
a separate question, and `scripts/epitope_overlap.py` answers it: the cetuximab
contact set in 6ARU comprises 24 receptor residues within 4.5 Å of the two Fab
chains (349, 350, 353, 382, 384, 408–418, 438–443 and 465–473), and the share
of each design's epitope that falls inside it is

| Design | epitope | shared with cetuximab | share |
| --- | --- | --- | --- |
| EGFRpH_C_l55_s661569_mpnn4 | 27 | 16 | 59 % |
| EGFRpH_C_l55_s661569_mpnn6 | 27 | 16 | 59 % |
| EGFRpH_C_l54_s137423_mpnn3 | 25 | 13 | 52 % |
| EGFRpH_B2_l53_s725081_mpnn13 | 40 | 9 | 22 % |
| EGFRpH_B2_l53_s725081_mpnn9 | 40 | 9 | 22 % |
| EGFRpH_C_l45_s295112_mpnn1 | 18 | 4 | 22 % |
| EGFRpH_C_l45_s295112_mpnn2 | 19 | 4 | 21 % |

**So the designs do not all bind the cetuximab epitope.** Three of the seven —
the s661569 pair and s137423_mpnn3 — share a majority of their contacts with
it. The other four overlap it only marginally and bind a largely different face
of domain III; they block the Fab sterically through its bulk rather than by
occupying its footprint. An earlier draft of this document asserted a shared
epitope for all seven on the strength of the steric overlap alone, which the
measurement does not support.

The epitope sizes in that table are larger than the per-design contact lists
given in `data/raw_outputs.md` (40 against 31 for s725081_mpnn13) because they
are computed against the full receptor chain of 6ARU after superposition,
whereas the earlier lists were computed against the design's own target patch.

Across the set the design epitopes span mature residues 316–326, 342–359,
379–385, 405–420 and 438–444, with the s295112 pair at 316–359 plus 384 (and
409 for mpnn2) and the s661569 and s137423 designs reaching 408–444.

Superposition onto 1IVO fails as expected (16.9 Å RMSD over the receptor chain)
because 1IVO is the extended, EGF-bound conformation while the other two are
tethered. The accessibility of these epitopes in the extended conformation was
therefore not tested, although the assay construct is the full ectodomain and
samples both states.

**Hotspots.** The epitope of the two run B designs ends at residue 410 and does
not include their specified hotspots 431, 434 and 436. For those designs the
hotspot specification did not take effect, and no BindCraft filter captures
that.
This is consistent with how BindCraft defines its metrics: `Hotspot_RMSD`
(1.58 Å for `s725081_mpnn13`, 1.28 Å for `s725081_mpnn9`) is assigned from
`unaligned_rmsd(trajectory_pdb, mpnn_design_pdb, binder_chain, binder_chain)`,
i.e. the RMSD of the binder chain between the trajectory pose and the
re-predicted ProteinMPNN design. As computed at that call site it does not
evaluate the specified hotspot residues, so for this campaign it carried no
information about whether hotspot targeting succeeded; the separate
`Target_RMSD` column measures target fidelity. Read from the BindCraft source
at commit `7713aa0`, 21 September 2026. Four of the five run C designs contact two of
their three specified hotspots; `s661569_mpnn4` contacts one.

**Glycosylation of the receptor.** N328, an EGFR N-glycosylation site, lies
4.1 Å from the interface of both run B designs; N420 lies 4.7 Å from
`s661569_mpnn6`. If mammalian-expressed receptor is used for the measurement,
this is a risk. Glycans are flexible, so it is not a definitive exclusion.

**Redundancy.** Three pairs share a backbone at 80–89 % identity; across pairs,
identity is 15–32 %. The submission therefore contains four distinct molecules
on seven slots, all directed at the same region of domain III.

## 7. Limitations

- No sequence was folded without its target; monomer stability is unverified.
  This would have been inexpensive to check and was omitted for lack of GPU time
  on the final day, not for any principled reason.
- No structure predictions against HER2, HER3 or HER4 were run, so no
  cross-reactivity assessment exists. A sequence-level comparison shows that
  36–53 % of the contacted EGFR residues are **conserved** in the paralogues,
  which is a risk rather than a clearance. The comparison is unweighted by
  buried surface or per-residue energy, is pooled across the three paralogues,
  Per design the proportion of differing contacted
  residues ranges from 46.9 % to 64.3 %.
- The steric check uses one crystal structure in one conformation and treats
  the receptor as rigid; see also the three caveats in section 4.
- PROPKA is an empirical method applied to a fixed backbone; it does not model
  conformational change on binding or explicit solvent. More importantly, the
  **entire observed spread of ΔΔG values, about 1.0 kcal/mol, is comparable to
  the method's own uncertainty.** No design is demonstrably better than another
  on this axis, and −0.45 kcal/mol may not be distinguishable from zero. No
  uncertainty estimate, replicate run or second pKa predictor was used to test
  this. The negative result — that no design reaches −1.3 kcal/mol — is robust
  to that uncertainty; the ranking among designs is not.
- i_pTM, pLDDT and the BindCraft energies originate from the pipeline that
  produced the designs and are not an independent assessment. **No independent
  structural rescoring was performed either:** the Rosetta recomputation in
  section 6 uses the same mover and score function BindCraft uses internally,
  so the discrepancy for five of seven designs is most plausibly an artefact of
  pose preparation. This remains unresolved.
- Only 7 of the 20 permitted slots are used, and all seven target the same
  surface. Epitope diversification was not attempted.
- Nothing is reported about expression construct, tags, linkers or purification
  feasibility.

## 8. Software versions

| Component | Version |
| --- | --- |
| PyRosetta | PyRosetta4.Release.python310.ubuntu, 2026.29+release.80a0635 |
| Python | 3.10 (conda environment `BindCraft`) |
| Operating system | Windows with WSL 2, Ubuntu |
| BindCraft | *not recorded* |
| AlphaFold2 weights | *not recorded* |
| ProteinMPNN weights | `v_48_020`, soluble weight set |
| PROPKA | *not recorded* |
| NCBI BLAST | web interface, blastp against ClusteredNR, 3 October 2026 |

Entries marked *not recorded* were not captured at the time of the runs.

## 9. Software and structures cited

**Design**

- Pacesa, M., Nickel, L., Schellhaas, C. et al. One-shot design of functional protein binders with BindCraft. *Nature* **646**, 483–492 (2025). doi:10.1038/s41586-025-09429-6
- Dauparas, J., Anishchenko, I., Bennett, N. et al. Robust deep learning-based protein sequence design using ProteinMPNN. *Science* **378**, 49–56 (2022). doi:10.1126/science.add2187
- Goverde, C. A., Pacesa, M., Goldbach, N. et al. Computational design of soluble and functional membrane protein analogues. *Nature* **631**, 449–458 (2024). doi:10.1038/s41586-024-07601-y — source of the solubility-optimised ProteinMPNN weights
- Jumper, J., Evans, R., Pritzel, A. et al. Highly accurate protein structure prediction with AlphaFold. *Nature* **596**, 583–589 (2021). doi:10.1038/s41586-021-03819-2
- Evans, R., O'Neill, M., Pritzel, A. et al. Protein complex prediction with AlphaFold-Multimer. *bioRxiv* 2021.10.04.463034 (v1 2021, v2 2022). doi:10.1101/2021.10.04.463034 — preprint, no journal version; the version corresponding to the weights used is not recorded
- Ovchinnikov, S. et al. ColabDesign. https://github.com/sokrypton/ColabDesign, commit and access date *not recorded* — no publication exists; used for direct access to ProteinMPNN

**Evaluation**

- Alford, R. F., Leaver-Fay, A., Jeliazkov, J. R. et al. The Rosetta All-Atom Energy Function for Macromolecular Modeling and Design. *J. Chem. Theory Comput.* **13**, 3031–3048 (2017). doi:10.1021/acs.jctc.7b00125
- Chaudhury, S., Lyskov, S., Gray, J. J. PyRosetta: a script-based interface for implementing molecular modeling algorithms using Rosetta. *Bioinformatics* **26**, 689–691 (2010). doi:10.1093/bioinformatics/btq007
- Stranges, P. B., Kuhlman, B. A comparison of successful and failed protein interface designs highlights the challenges of designing buried hydrogen bonds. *Protein Sci.* **22**, 74–82 (2013). doi:10.1002/pro.2187 — the reference given by the Rosetta documentation for InterfaceAnalyzer
- Lawrence, M. C., Colman, P. M. Shape complementarity at protein/protein interfaces. *J. Mol. Biol.* **234**, 946–950 (1993). doi:10.1006/jmbi.1993.1648
- Olsson, M. H. M., Søndergaard, C. R., Rostkowski, M., Jensen, J. H. PROPKA3: Consistent Treatment of Internal and Surface Residues in Empirical pKa Predictions. *J. Chem. Theory Comput.* **7**, 525–537 (2011). doi:10.1021/ct100578z
- Søndergaard, C. R., Olsson, M. H. M., Rostkowski, M., Jensen, J. H. Improved Treatment of Ligands and Coupling Effects in Empirical Calculation and Rationalization of pKa Values. *J. Chem. Theory Comput.* **7**, 2284–2295 (2011). doi:10.1021/ct200133y

**Sequence analysis and libraries**

- Altschul, S. F., Gish, W., Miller, W., Myers, E. W., Lipman, D. J. Basic local alignment search tool. *J. Mol. Biol.* **215**, 403–410 (1990). doi:10.1016/S0022-2836(05)80360-2
- Altschul, S. F., Madden, T. L., Schäffer, A. A. et al. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. *Nucleic Acids Res.* **25**, 3389–3402 (1997). doi:10.1093/nar/25.17.3389 — the algorithm behind the `blastp` service actually used here
- NCBI BLAST web service, https://blast.ncbi.nlm.nih.gov/, queried 3 October 2026. The standalone BLAST+ package (Camacho et al. 2009) was **not** used and is therefore not cited.
- ClusteredNR, NCBI database resource. No standalone publication; described in NCBI's annual database-resources paper. Queried through the BLAST web interface on 3 October 2026. *The database release version was not recorded. A negative similarity result is interpretable only against a known database snapshot, which limits the strength of the novelty statement.*
- Kyte, J., Doolittle, R. F. A simple method for displaying the hydropathic character of a protein. *J. Mol. Biol.* **157**, 105–132 (1982). doi:10.1016/0022-2836(82)90515-0 — the hydropathy scale behind the GRAVY values
- Kabsch, W. A solution for the best rotation to relate two sets of vectors. *Acta Crystallogr. A* **32**, 922–923 (1976). doi:10.1107/S0567739476001873 — the superposition algorithm used for every structural superposition reported here, implemented directly in `final_checks.py` and `epitope_overlap.py`
- Berman, H. M., Westbrook, J., Feng, Z. et al. The Protein Data Bank. *Nucleic Acids Res.* **28**, 235–242 (2000). doi:10.1093/nar/28.1.235
- The UniProt Consortium. UniProt: the Universal Protein Knowledgebase. *Nucleic Acids Res.*, current annual database issue. Entry P00533 accessed 3 October 2026; *release version not recorded*.
- Cock, P. J. A., Antao, T., Chang, J. T. et al. Biopython. *Bioinformatics* **25**, 1422–1423 (2009). doi:10.1093/bioinformatics/btp163
- Harris, C. R., Millman, K. J., van der Walt, S. J. et al. Array programming with NumPy. *Nature* **585**, 357–362 (2020). doi:10.1038/s41586-020-2649-2
- McKinney, W. Data Structures for Statistical Computing in Python. *Proc. 9th Python in Science Conf.*, 56–61 (2010). doi:10.25080/Majora-92bf1922-00a

**Structures**

- **6ARU** — Christie, M., Christ, D. Structure of Cetuximab Fab mutant in complex with EGFR extracellular domain. PDB entry, deposited 2017. RCSB DOI 10.2210/pdb6ARU/pdb. *No publication exists for this entry*; its primary reference is listed as "to be published", so it is cited here as a PDB deposition rather than as a journal article.
- **1YY9** — Li, S., Schmitz, K. R., Jeffrey, P. D. et al. Structural basis for inhibition of the epidermal growth factor receptor by cetuximab. *Cancer Cell* **7**, 301–311 (2005). doi:10.1016/j.ccr.2005.03.003. PDB entry deposited 2005, RCSB DOI 10.2210/pdb1YY9/pdb
- **1IVO** — Ogiso, H., Ishitani, R., Nureki, O. et al. Crystal structure of the complex of human epidermal growth factor and receptor extracellular domains. *Cell* **110**, 775–787 (2002). doi:10.1016/S0092-8674(02)00963-7. PDB entry deposited 2002, RCSB DOI 10.2210/pdb1IVO/pdb

EGFR sequence reference: UniProt P00533 (EGFR_HUMAN), isoform P00533-1.
Isoelectric points and GRAVY values were computed with Biopython's `ProtParam`
module (Bjellqvist pKa scales, Kyte–Doolittle hydropathy).

The analysis scripts in `scripts/` and the text of this document were produced
with substantial assistance from an AI assistant (Anthropic Claude,
September–October 2026); see the AI use statement in `README.md`. All outputs
were verified by the author against the raw tool output in
`data/raw_outputs.md`, and the author is solely responsible for the content.
