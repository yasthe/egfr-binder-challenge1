# Methods and design rationale

*Text for the "methodology" field of the Proteinbase submission.*

## Approach

Seven de novo mini-binders (45–55 aa) targeting domain III of the human EGFR
extracellular region. Backbones were generated with BindCraft, which optimises
binder geometry through AlphaFold2-multimer and assigns sequences with
ProteinMPNN. Four runs were carried out against three target patches cut from
PDB 6ARU (mature EGFR residues 300–450, 310–455 and 365–520), with a helicity
bias of −0.3 and the default four-stage multimer protocol. 47 designs passed
the BindCraft filters.

## Selection

All 47 designs were re-evaluated against the intact ectodomain in 6ARU rather
than against their design patch alone, because the design stage never sees the
rest of the receptor. Only 4 of 47 touch nothing outside their own patch;
applying a cut-off that tolerates van der Waals contact (minimum heavy-atom
distance ≥3.0 Å, ≤2 clashes) retains seven and discards 40. Four of the seven
stand 12.5–13.3 Å clear of everything outside their patch; the other three make
van der Waals contact with a neighbouring domain at 3.3–3.7 Å.

Interface quality was then recomputed with the PyRosetta InterfaceAnalyzer:
ΔG −27 to −56 REU, buried surface 1256–2139 Å², shape complementarity
0.56–0.73, 4–18 buried unsatisfied hydrogen bonds. This is the same mover
BindCraft uses internally, so it is a consistency check rather than an
independent assessment. For five of the seven the recomputed energy is less
favourable than BindCraft's, most plausibly because BindCraft scores a relaxed
pose and this recomputation did not. The discrepancy is reported here; it is
not claimed as independent confirmation or refutation. No design
contains cysteine; the exposed hydrophobic fraction of the free binder is
14–23 %. Two designs carry an N-glycosylation sequon of their own, which is
inconsequential if the protein is produced cell-free but not in a eukaryotic
host. blastp against ClusteredNR (3 October 2026) returns no
significant similarity for any of the seven, although for 45–55 residue queries
that is weak evidence of novelty.

## Epitope

Superposed into the EGFR–cetuximab complexes 6ARU and 1YY9, all seven designs
overlap the antibody heavy chain substantially (185–514 and 197–537 atoms within
5 Å respectively) while contacting the receptor outside their own patch with
only 0–10 atoms in those superposed frames. They would therefore compete with cetuximab for binding.
Cetuximab's own contact set was not computed, so no claim of an identical
epitope is made. Design contacts span mature EGFR residues 316–444 (UniProt P00533-1:
340–468), but are not uniform across the set. One caveat: N328, a receptor
N-glycosylation site, lies 4.1 Å from the interface of the two designs on the
s725081 backbone.

## pH dependence: a negative result

The challenge asks for pH-conditional binding, and this was not achieved.
Histidine pKa shifts on binding were estimated with PROPKA3 for 36 of the 47
accepted structures and converted into a binding free-energy difference between
pH 6.5 and 7.4. Values range from −0.45 to +0.52 kcal/mol, against roughly
−1.3 kcal/mol needed for a tenfold change in affinity. The six submitted designs
other than the best lie between +0.02 and +0.36 kcal/mol.

Three strategies were applied and all failed: histidine point mutations at
interface positions, increased histidine content, and full redesign of the
binder surface with a histidine bias in ProteinMPNN (bias levels 0.0–1.5, 160
sequences, 28 threaded and repacked; all 28 between +0.13 and +0.51 kcal/mol,
including the unbiased controls).

Where the approach breaks down is visible at residue level: interface histidine
pKa values fall on binding — in the variant examined, 6.41 to 4.05 and 6.10 to
2.51 — which favours the neutral form in the complex and therefore weakens
rather than strengthens binding at low pH. Two of those three histidines also
started below pH 6.5 in the free state, so no switch was available even before
binding was considered. The desired direction needs a free pKa near 6.5–7 that
rises on binding, for example a histidine positioned against a target
carboxylate. Neither BindCraft nor ProteinMPNN represents protonation states, so
this geometry does not arise by chance.

These designs are submitted as EGFR binders with no measurable pH dependence —
the whole observed spread lies within the uncertainty of the method — and not
as a demonstrated switch.

## Limitations

None of the seven sequences was folded without its target, so monomer stability
is unverified. No structure predictions against HER2, HER3 or HER4 were run; a
sequence-level comparison shows 36–53 % of contacted EGFR residues are conserved
in the paralogues, which is a cross-reactivity risk rather than a clearance. The
steric check rests on a single crystal structure in one conformation and treats
the receptor as rigid, and the epitope's accessibility in the extended, EGF-bound
conformation was not tested. PROPKA estimates are empirical and use a fixed
backbone; the entire observed spread of pH values is comparable to the method's
own uncertainty, so the negative result is robust but the ranking among designs
is not. No genuinely independent structural rescoring was performed. Seven of 20
permitted slots are used, on four distinct scaffolds, all directed at the same
region of domain III.

**AI use.** All design, analysis and interpretation are my own, and I take full
responsibility for this submission. An AI assistant (Anthropic Claude,
September–October 2026) helped write and debug the analysis scripts, draft the
documentation, and explore the numerical results; all scripts and all scientific
claims were checked by me against the raw outputs. The AI assistant is not an
author. The deep-learning design methods used (BindCraft, AlphaFold2,
ProteinMPNN) are described and cited in the methods.
