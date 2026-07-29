# *Oprl1* is the dominant opioid receptor transcript of mouse peripheral neurons, from neurogenesis onward

## Introduction

The nociceptin receptor is the most abundant opioid receptor transcript in mouse peripheral
neurons, in every division of the peripheral nervous system and from the earliest stage at which
those neurons can be identified. Work on peripheral opioid receptors has been organised around the
mu receptor for four decades, and the therapeutics that organisation produced are real:
loperamide, eluxadoline, and the peripherally acting mu antagonists methylnaltrexone, naloxegol,
naldemedine and alvimopan all act on enteric MOR, and difelikefalin, a peripherally restricted
kappa agonist, was approved in 2021 for pruritus in adults on haemodialysis. Peripheral
restriction is the design principle behind each of them, and it makes the identity of the receptor
these neurons carry a question with a direct consequence for which restricted agonist is worth
building.

NOP, encoded by *OPRL1* at 20q13.33, shares roughly 47% amino acid identity with MOR, DOR and KOR
overall and 64 to 67% within the transmembrane domains (Mollereau et al., *FEBS Lett* 1994;
Meunier et al., *Nature* 1995). The divergence is concentrated where it matters for pharmacology.
A glutamine at position 280 in TM6 replaces the histidine conserved in the classical receptors,
and acidic residues in ECL2 (Glu194, Glu199) anchor the basic core of N/OFQ (Thompson et al.,
*Nature* 2012). The endogenous ligand begins Phe-Gly-Gly-Phe against the Tyr-Gly-Gly-Phe of the
enkephalins and dynorphins. NOP is therefore insensitive to naloxone and to the morphinans, and
the classical opioid peptides are weak at it. Downstream the four receptors converge: NOP couples
to Gi/Go, inhibits Cav2.2, activates GIRK, and suppresses TRPV1 responses in cultured human DRG
neurons at picomolar concentrations of N/OFQ (Anand et al., *PAIN* 2016). A ligand acting there
produces opioid-like inhibition of neuronal output while leaving the receptor responsible for
reward and respiratory depression untouched. Cebranopadol, a dual NOP/MOP agonist, met its primary
endpoint in two Phase 3 acute pain trials in 2025; sunobinop, an oral selective NOP partial
agonist, is in Phase 2 for insomnia and Phase 1b for interstitial cystitis and overactive bladder.

The evidence that NOP is abundant in peripheral neurons is protein-level and sparse. Anand and
colleagues reported NOP immunoreactivity in 75 to 80% of small and medium neurons in human lumbar
and sacral DRG, with a several-fold increase in NOP-positive suburothelial fibres in detrusor
overactivity and painful bladder syndrome. NOP-eGFP reporter mice place the receptor in laminae I
to III and in small-diameter DRG somata alongside CGRP and MOR. Functional potency for NOP-
mediated inhibition of N-type calcium current spans two orders of magnitude across peripheral
ganglia: approximately 0.5 nM in cervical sympathetic neurons, 26 nM in vestibular afferents and
100 nM in DRG. That range implies a corresponding difference in receptor expression across
ganglia, and it has never been measured. The immunohistochemistry carries a separate problem.
Opioid receptor antibodies frequently label knockout tissue, and protein localisation claims that
rest on them are provisional until reproduced with knock-in reporters, validated in situ
hybridisation, autoradiography or PET.

Single-cell transcriptomics can supply the missing survey, and a preparation artefact governs how
its output should be read. *Oprm1* spans roughly 250 kb of mouse genome; *Oprl1* spans 6 kb.
Single-nucleus libraries retain unspliced pre-mRNA, so genes with long introns gain signal in
proportion to the intronic sequence they carry. In mouse nodose ganglion, the one tissue for which
whole-cell and nuclear data exist side by side, nuclear preparation multiplies *Oprm1* by 29 and
*Oprl1* by 0.47, which moves *Oprm1* from second place to first. The same reversal occurs in the
dorsal root ganglion: whole-cell data place *Oprl1* first at 9.31 CPM against *Oprm1* 8.23, and
nuclei from the same tissue in a different laboratory place *Oprm1* first at 94.92 against *Oprl1*
7.60. Human peripheral ganglia are recovered post mortem or surgically and are frozen before
dissociation, so essentially every human single-cell dataset from them is single-nucleus and
subject to this effect.

We recomputed the four opioid receptors from source matrices across 19 peripheral neuronal
populations, holding QC, neuron definition, marker gate and ambient-RNA check constant across all
of them. The populations comprise eight sensory populations spanning all three developmental
origins of the peripheral sensory system, six sympathetic populations drawn from five ganglia, the
pelvic ganglion, two parasympathetic ganglia, and enteric submucosal neurons at two postnatal
ages. Intervals come from resampling the datasets rather than the cells, since the question is
whether a result reproduces across the peripheral nervous system rather than whether it reaches
significance in one ganglion.

*Oprl1* is the highest-expressed opioid receptor in 18 of the 19 peripheral populations. The
exception is enteric submucosal neurons at postnatal day 7, where *Oprk1* leads by 4.59-fold and
falls to a thirtieth of that level by day 24. The ordering is present at embryonic day 9.5 in the
earliest otic neuroblasts and holds at every subsequent age to adulthood, while *Oprm1* and
*Oprk1* are transiently expressed in the embryonic ganglion and extinguished before maturity.
Against the wider denominator of every G protein-coupled receptor measured in the same cells,
*Oprl1* holds a median rank of 31st of about 280 detected, the 90th percentile, while *Oprm1* and
*Oprd1* sit at the 48th and 47th percentiles. *Oprl1* occupies the top decile of the receptor
repertoire of these neurons and stands outside their twenty most abundant receptors, and it is the
only opioid receptor above the median. Genes co-expressed with *Oprl1* across thirteen populations
form a synaptic and axonal adhesion programme, and the association with *Scn1a* reported in the
vagus and the dorsal root ganglion holds between populations that differ in sodium channel class
and disappears within a population that does not.

The human question stays open, and the preparation is what leaves it open. In a cross-species
atlas of dorsal root ganglion nuclei covering human, macaque, mouse and guinea pig on one
protocol, the mouse arm reverses in the direction the artefact predicts, so no human single-
nucleus dataset can establish a species difference. The only non-nuclear human peripheral ganglion
data located for this work, bulk RNA-seq of six superior cervical ganglia, places *OPRL1* last of
the four at 0.45 TPM, in tissue where mouse gives *Oprl1* 18.68 CPM against *Oprm1* 0.27. Human
NOP protein is reported in three-quarters of small and medium DRG neurons by immunolabelling, by a
reagent class that requires knockout validation it has not received. Resolving the discrepancy
requires in situ hybridisation or a validated protein readout in human ganglion tissue, and the
results below are mouse results until that exists.
