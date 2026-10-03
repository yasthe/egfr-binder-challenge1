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
PROPKA3 estimates span -0.45 to +0.52 kcal/mol across the set, against roughly
-1.3 needed for a tenfold change. In the cases examined residue by residue, the
pKa of interface histidines falls on binding instead of rising, which favours
the neutral form in the complex, and several histidines sit below pH 6.5 even
in the free state. Neither BindCraft nor ProteinMPNN represents protonation, so
the required geometry does not arise by chance. 11 of the 47 designs contain no
histidine at all.

What I am submitting is therefore seven EGFR binders without a working pH
switch. What may still make them worth testing is the second filter they
passed. Every design was re-checked against the intact ectodomain in 6ARU, not
only against its own design patch: 40 of 47 collide with receptor regions the
design stage never sees, and only these seven survive. They span four
independent backbones and two distinct binding modes - three share 52-59
percent of their contacts with the cetuximab epitope, the other four bind a
largely different face of domain III at 21-22 percent shared. Recomputed
interface metrics: dG -27 to -56 REU, buried surface 1256-2139 A^2, shape
complementarity 0.56-0.73. No cysteines, no significant blastp hit.

Known gaps, in case they bear on the decision: no target-free folding, so
monomer stability is unverified; no predictions against HER2, HER3 or HER4,
where 36-53 percent of the contacted residues are conserved; the steric check
uses one rigid crystal structure in one conformation; and the whole pH spread
lies within PROPKA's own uncertainty, so the negative result holds but the
ranking among the seven does not.
