This entry was produced with an AI assistant (Anthropic Claude) writing all
code and analysis, while I set the objective, ran every computation and made
the decisions; scripts and raw outputs are public at
github.com/yasthe/egfr-binder-challenge1.

The hypothesis was the obvious one for this challenge: place histidines at the
interface so that the protonated form is favoured in the bound state, giving
tighter binding at pH 6.5 than at 7.4. I generated 47 mini-binders of 45-69 aa
against domain III of the EGFR ectodomain with BindCraft, using three target
patches cut from 6ARU, and then tried three ways of installing a switch:
histidine point mutations at interface positions, higher histidine content, and
full redesign of the binder surface with a histidine bias in ProteinMPNN.

None of it worked, and the reason looks specific rather than incidental.
PROPKA3 estimates span -0.45 to +0.52 kcal/mol across the 47, against roughly
-1.3 needed for a tenfold change. In the cases examined residue by residue, the
pKa of interface histidines falls on binding instead of rising, which favours
the neutral form in the complex, and several histidines sit below pH 6.5 even
in the free state. Neither BindCraft nor ProteinMPNN represents protonation, so
the required geometry does not arise by chance. 11 of the 47 designs contain no
histidine at all.

I am therefore submitting three EGFR binders with no pH dependence worth the
name (+0.02 to +0.08 kcal/mol). What they do have is a second filter. Every
design was re-checked against the intact ectodomain in 6ARU, not only against
its own design patch, and 40 of 47 collide with receptor regions that the
design stage never sees. Seven survived that check; four of those fall below
the novelty threshold here, on the structural axis rather than the sequence
axis, and were withdrawn. The three remaining span two backbones and two
binding modes: two share 59 percent of their contacts with the cetuximab
epitope, one binds a different face of domain III at 22 percent shared.
Recomputed interface metrics are dG -31 to -36 REU, buried surface 1281-1896
A^2, shape complementarity 0.56-0.73, and none contains cysteine.

Known gaps: monomer stability is unverified, as no sequence was folded without
its target; no predictions were run against HER2, HER3 or HER4, where 36-53
percent of the contacted residues are conserved; the steric check uses a single
rigid crystal structure in one conformation; and the whole pH spread lies
within PROPKA's own uncertainty.
