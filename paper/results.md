# Results

## *Oprl1* is the most broadly expressed opioid receptor gene across the peripheral nervous system

Every peripherally restricted opioid drug in clinical use acts at MOR or KOR. The four opioid
receptors have rarely been measured against one another in the same peripheral neurons. We ranked
*Oprl1*, *Oprm1*, *Oprd1*, and *Oprk1* within each of 19 neuronal populations from 16 mouse ganglia
and plexuses, covering the sensory, sympathetic, parasympathetic, and enteric divisions. The data
come from 12 published single-cell studies deposited between 2016 and 2026; we generated no new
sequencing (Table S1). We took the whole-cell fraction of each deposit and applied one analysis to
all of them. *Oprl1* carries the highest population mean of the four in 18 of the 19 populations.

The cell bootstrap reported throughout describes how far a population's ordering depends on which
cells were captured. It is not a test across animals, and four of the nineteen populations come
from a single one, so the inference is made at the unit the claim is about. Taking the population
as that unit, *Oprl1* leads 18 of 19 (sign test against a coin-flip null, p = 3.8 x 10^-5). Taking
the deposit, so that the five autonomic ganglia of one experiment count once and a deposit leads
only if *Oprl1* is first in every population it contributes, 11 of 12 (p = 0.0032). The one deposit
that does not lead is GSE263422, on the strength of its postnatal day 7 enteric sample.

The ranking holds in SMART-seq full-length libraries sequenced to a median of 2.7 million reads per
neuron and in droplet libraries at 1,570 UMI per cell. Separation from the other three receptors is
widest in the spiral ganglion, where *Oprl1* reaches 41.94 CPM in 74.8% of neurons and *Oprm1*,
*Oprd1*, and *Oprk1* are absent from all 226 cells, and in the sphenopalatine ganglion, where
*Oprl1* reaches 19.50 CPM against *Oprm1* at 0.29, *Oprd1* at 0.04, and *Oprk1* at 0.01 CPM
(Figures 1A and 1P). The two largest droplet datasets give the narrowest separations: *Oprl1*
exceeds *Oprm1* by 1.21-fold in 26,047 nodose neurons and by 1.13-fold in 31,802 dorsal root
ganglion neurons, each holding in every bootstrap replicate with a 95% confidence interval
excluding 1 (Figures 1F and 1H).

The geniculate ganglion was sequenced full-length by Dvoryanchikov et al. (2017) and by droplet
capture by Zhang et al. (2019), and the stellate ganglion appears in two unrelated studies. *Oprl1*
ranks first in all four (Figures 1B, 1D, 1N, and 1O).

Five populations were dissected and sequenced together in one experiment (Sivori et al., 2024),
which places them on a common scale. Across four sympathetic ganglia and the sphenopalatine
ganglion, *Oprl1* varies 1.7-fold, from 16.33 CPM in the pelvic ganglion to 27.55 CPM in the
stellate ganglion. The other three receptors vary between 6-fold and 302-fold over the same five
populations. Across the full panel, *Oprm1*, *Oprk1*, and *Oprd1* each rank second in some
populations, and none of them ranks second consistently. NOP is available to autonomic neurons at a
similar level regardless of ganglion, and MOR, DOR, and KOR are distributed by division and cell
type.

The enteric submucosal plexus before weaning ranks *Oprk1* first. At postnatal day 7, *Oprk1*
reaches 122.51 CPM against 26.70 CPM for *Oprl1* in both samples at that age (Figure 1S). By
postnatal day 24, *Oprk1* has fallen to 3.77 CPM at 1.1% detection, *Oprl1* stands at 35.57 CPM,
and *Oprl1* ranks first in all three samples (Figure 1R). Both ages come from one deposit, one
tissue, and one platform (Li et al., 2025), and the 32-fold fall in *Oprk1* is the largest change
of the four receptors between them: *Oprd1* rises 7-fold, *Oprm1* falls 4-fold, and *Oprl1* moves
1.3-fold. High submucosal *Oprk1* is a property of the neonatal plexus, and the mature plexus
matches the rest of the panel.

The difference between *Oprl1* and the other three is principally one of prevalence. In the 15
populations that record detection for all four receptors, *Oprl1* is the most widely detected in
14, reaching 7.2% to 91.7% of neurons against 0.2% to 12.0% for *Oprm1* (Figure S1A). Among the
neurons that do express them the receptors sit much closer together: the median margin falls from
3.97-fold on the population mean to 1.17-fold on the mean among positive cells, and in 2 of the 8
populations where a second receptor is detected in at least 100 neurons the ordering changes, the
vestibular ganglion and the jugular (Figures S1B and S1C). A mean over positive cells is floored by
the one-count detection limit and therefore favours the rarely detected receptors, so that
comparison is conservative for prevalence and unreliable where the positive-cell count is small.
*Oprl1* is available to more peripheral neurons than the other three opioid receptors rather than
present at a higher level in the neurons that carry them.

## Nuclear libraries scale with gene length and whole-cell libraries do not

Single-nucleus libraries from vagal and dorsal root ganglia place *Oprm1* first, at 244.49 CPM
against 5.22 for *Oprl1* and at 94.92 against 7.60, and the whole-cell fractions of those same
atlases place *Oprl1* first (Figures 2A and 2B). The two preparations sample different RNA. A
nucleus holds transcript that is still being made, unspliced and not yet exported, so most of a
nuclear library is intronic; single-nucleus quantification counts reads across the entire gene body
because counting exons alone would discard that majority. The number of reads assigned to a gene
then scales with the length of its transcription unit as well as with the number of transcripts the
neuron has made. Whole-cell libraries sample the cytoplasmic pool, which is the spliced mRNA
available for translation, and 3' counting of that pool carries no length term. Across 14,876
protein-coding genes quantified in both preparations, nuclear expression rises with genomic span
(r = +0.403) and whole-cell expression stays flat over a 68-fold range of length (r = +0.022;
Figure 2E). *Oprm1* spans 279.7 kb and *Oprl1* spans 7.1 kb, giving nuclear-to-whole-cell ratios
of 25.4 and 0.45 (Figure 2F), and *Oprm1* is detected in 71.9% of nuclei against 5.5% to 32.4% of
cells (Figures 2C and 2D). The nuclear values report a combination of transcript abundance and
gene length.

Intronic sequence is the component that carries the length effect. Genomic span stands in for the
pre-mRNA a nucleus holds, and the two can be separated: exonic length is the union of a gene's exons
over all its transcripts, and intronic length is the remainder of the span. Median whole-cell
expression holds between 14.1 and 19.6 CPM across deciles of intronic length whose medians run from
0.93 to 229.2 kb, a 246-fold range, while nuclear expression rises from 8.0 to 91.9 CPM over the
same deciles (Figure 2G). Over the individual genes, intronic length correlates with the
nuclear-to-whole-cell ratio at r = +0.474 and exonic length at r = +0.377; holding the other fixed,
intronic length gives +0.362 and exonic length +0.197. Neither reaches the +0.592 of the whole span,
and length in any form accounts for 35% of the variance at most, so the two components do not
partition the effect cleanly and neither is the whole of it. The union of exons overstates the
mature transcript of any one isoform and understates introns by the same amount, which works against
the intronic account rather than for it.

The four receptors separate the two components sharply. *Oprm1* carries 271.7 kb of intron and
*Oprd1* 29.8 kb, against 8.6 kb in *Oprk1* and 3.4 kb in *Oprl1* (Figure 2H). Ordered by intronic
length the four run *Oprm1*, *Oprd1*, *Oprk1*, *Oprl1*, which is their order by the
nuclear-to-whole-cell ratio exactly; ordered by exonic length they do not, because *Oprk1* carries
the most exonic sequence of the four at 9.2 kb and gains least of the four in nuclei, at a ratio of
1.14. *Oprl1* and *Oprk1* are intron-poor, at fractions of 0.47 and 0.48 against a genome-wide
median of 0.84, while *Oprm1* and *Oprd1* are intron-rich at 0.97 and 0.88. That is why *Oprl1* is
depleted in nuclei rather than merely unchanged: it has little unspliced sequence to contribute, and
the transcripts a nucleus does hold are counted against a gene body that is mostly exon.

Dissociation is a second preparation effect, and it is separable from this one. Whole-cell
libraries are made from live enzymatic dissociation, which axotomises every neuron, and they carry
that injury: *Fos* reads 448.82 CPM in whole cells against 0.87 in nuclei, and *Atf3*, the axotomy
marker, 151.72 against 6.15. Scored against genes of the same length, the twelve dissociation-
induced genes sit at a median 1.9th percentile of their own length class, and 11 of the 12 fall
below the lowest of the four receptors. The receptors sit at the 31.8th, 60.1st, 93.6th and 99.2nd
percentiles (Figure S2). No claim here may assume the whole-cell libraries are unperturbed, and the
length relationship is not a restatement of the injury response.

Several of the largest sensory ganglion atlases are single-nucleus, which is consistent with, and
may have reinforced, the prevailing view that *Oprm1* is the dominant opioid receptor of peripheral
neurons. That view predates those atlases and rests on pharmacology, conditional knockouts and
reporter lines rather than on transcript counts.
