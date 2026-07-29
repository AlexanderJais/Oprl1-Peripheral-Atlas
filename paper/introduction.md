# *Oprl1* is the dominant opioid receptor transcript of mouse peripheral neurons, from neurogenesis onward

## Introduction

The nociceptin receptor is the most abundant opioid receptor transcript in mouse peripheral
neurons, in every division of the peripheral nervous system and from the earliest stage at which
those neurons can be identified.

Peripheral opioid pharmacology has been built on the mu receptor. Loperamide and eluxadoline act
on enteric MOR, as do the peripherally acting antagonists methylnaltrexone, naloxegol, naldemedine
and alvimopan. Difelikefalin, a peripherally restricted kappa agonist, was approved in 2021 for
pruritus in adults on haemodialysis. Peripheral restriction is the design principle behind all of
them, and it turns receptor identity into a practical constraint: a restricted agonist is worth
building for the receptor the target neurons carry.

NOP, encoded by *OPRL1* at 20q13.33, shares roughly 47% amino acid identity with MOR, DOR and KOR
overall and 64 to 67% within the transmembrane domains (Mollereau et al., *FEBS Lett* 1994;
Meunier et al., *Nature* 1995). The divergence is concentrated where it determines pharmacology. A
glutamine at position 280 in TM6 replaces the histidine conserved in the classical receptors, and
acidic residues in ECL2 (Glu194, Glu199) anchor the basic core of N/OFQ (Thompson et al., *Nature*
2012). The endogenous ligand begins Phe-Gly-Gly-Phe against the Tyr-Gly-Gly-Phe of the enkephalins
and dynorphins. NOP is therefore insensitive to naloxone and to the morphinans, and the classical
opioid peptides are weak at it. Downstream the four receptors converge: NOP couples to Gi/Go,
inhibits Cav2.2, activates GIRK, and suppresses TRPV1 responses in cultured human DRG neurons at
picomolar concentrations of N/OFQ (Anand et al., *PAIN* 2016). A NOP agonist can therefore inhibit
peripheral neuronal output without engaging the receptor that carries the reward and respiratory
liabilities. Cebranopadol, a dual NOP/MOP agonist, met its primary endpoint in two Phase 3 acute
pain trials in 2025; sunobinop, an oral selective NOP partial agonist, is in Phase 2 for insomnia
and Phase 1b for interstitial cystitis and overactive bladder.

The evidence that NOP is abundant in peripheral neurons is protein-level and sparse. Anand and
colleagues reported NOP immunoreactivity in 75 to 80% of small and medium neurons in human lumbar
and sacral DRG, with a several-fold increase in NOP-positive suburothelial fibres in detrusor
overactivity and painful bladder syndrome. NOP-eGFP reporter mice place the receptor in laminae I
to III and in small-diameter DRG somata alongside CGRP and MOR. Potency for NOP-mediated
inhibition of N-type calcium current spans two orders of magnitude across peripheral ganglia:
approximately 0.5 nM in cervical sympathetic neurons, 26 nM in vestibular afferents and 100 nM in
DRG. A range that wide implies a corresponding difference in receptor expression between ganglia,
and no systematic survey of *Oprl1* across the peripheral nervous system exists to test it. The
immunohistochemistry carries a separate problem. Opioid receptor antibodies frequently label
knockout tissue, and protein localisation claims resting on them stay provisional until reproduced
with knock-in reporters, validated in situ hybridisation, autoradiography or PET.

Single-cell transcriptomics can supply the missing survey. Library preparation constrains how its
output is read. *Oprm1* spans roughly 250 kb of mouse genome and *Oprl1* spans 6 kb, and
single-nucleus libraries retain unspliced pre-mRNA, so genes gain signal in proportion to the
intronic sequence they carry. In mouse nodose ganglion, the one tissue for which whole-cell and
nuclear data exist side by side, nuclear preparation multiplies *Oprm1* by 29 and *Oprl1* by 0.47,
which moves *Oprm1* from second place to first. Human peripheral ganglia are recovered post mortem
or surgically and frozen before dissociation, so nearly every human single-cell dataset from them
is single-nucleus and subject to this effect.

We recomputed the four opioid receptors from source matrices across 19 peripheral neuronal
populations: eight sensory populations spanning all three developmental origins of the peripheral
sensory system, six sympathetic populations drawn from five ganglia, the pelvic ganglion, two
parasympathetic ganglia, and enteric submucosal neurons at two postnatal ages. One QC policy, one
neuron definition, one marker gate and one ambient-RNA check apply to all of them. Confidence
intervals come from resampling the datasets rather than the cells, because the claim under test is
reproducibility across the peripheral nervous system. The comparison is made in mouse tissue,
where whole-cell data exist for every ganglion examined.
