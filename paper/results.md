# Results

## *Oprl1* is the highest-expressed opioid receptor gene across the peripheral nervous system

We measured the four opioid receptor transcripts in 19 neuronal populations from 16 mouse
peripheral ganglia and plexuses, using whole-cell data throughout and a single analysis applied to
every dataset (Figure 1; Table S1). *Oprl1* is the highest-expressed of the four in 18 of the 19
populations. The exception is enteric submucosal neurons at postnatal day 7, treated below.

The margin over the second receptor varies with the ganglion and with sequencing depth. In the
spiral ganglion, sequenced to a median of 2.7 million reads per neuron, *Oprl1* reaches 41.94 CPM
in 74.8% of neurons and *Oprm1*, *Oprd1*, and *Oprk1* are detected in none of the 226 cells
(Figure 1A). All four genes are quantified in that annotation, so the three zeros are measured
absences. In the geniculate ganglion, *Oprl1* exceeds *Oprm1* 31-fold in full-length data (Figure
1B) and reproduces on a second platform from a second laboratory (Figure 1D, bootstrap support
0.974). The narrowest margins occur in the two largest droplet datasets: 1.21-fold in 26,047
nodose neurons and 1.13-fold in 31,802 dorsal root ganglion neurons, both with support 1.000 and
95% intervals excluding 1 (Figures 1F and 1H).

Sympathetic neurons carry the pattern at larger margins than most sensory ganglia. *Oprl1* exceeds
the second receptor 44.41-fold in celiac neurons, 32.42-fold in the pelvic ganglion, 25.72-fold in
the superior cervical ganglion, and 14.78- and 14.69-fold in the lumbar and thoracic chains, and
exceeds *Oprm1* in all seven, from 4.9-fold in one stellate dataset to 69-fold in the superior
cervical ganglion (Figures 1I to 1O). Six sympathetic ganglia from four deposits on two platform
classes agree, the stellate ganglion is measured twice from two deposits (Figures 1N and 1O), and
every sample that could be ordered on its own cells places *Oprl1* first.

The parasympathetic division gives the largest margin measured. In 2,014 sphenopalatine neurons
*Oprl1* reaches 19.50 CPM against *Oprm1* 0.29, *Oprd1* 0.04, and *Oprk1* 0.02, a 66.86-fold lead
over the second receptor (Figure 1P). These neurons are cholinergic and not noradrenergic
(*Slc18a3* 252 CPM, *Chat* 26, *Th* 2, *Dbh* 42), against *Th* 390 to 840 CPM in the sympathetic
ganglia dissected and sequenced alongside them in the same experiment. Intrinsic cardiac neurons
give a 3.65-fold margin, the narrowest in this division and the only population in which *Oprm1*
is the second receptor rather than a trace transcript (Figure 1Q).

Five panels of Figure 1 come from one experiment in which six autonomic ganglia were dissected and
sequenced together (Figures 1I, 1J, 1L, 1N, and 1P), which removes laboratory, platform, and
sequencing depth from the comparison between divisions. *Oprl1* is first in all five and its level
varies 1.7-fold across them, from 16.33 CPM in the pelvic ganglion to 27.55 CPM in the stellate
ganglion, while the margin over the second receptor varies 11-fold, from 6.08-fold to 66.86-fold.
The margin therefore reports how far the other three receptors fall in a given ganglion rather
than how high *Oprl1* rises. The sphenopalatine ganglion carries the largest margin of the five,
with the lowest *Oprm1* and the lowest *Oprk1* among them.

Enteric submucosal neurons place *Oprk1* first at postnatal day 7, at 122.51 CPM against *Oprl1*
26.70, a 4.59-fold lead in the opposite direction, in both samples at that age (Figure 1S). By
postnatal day 24 *Oprk1* has fallen to 3.77 CPM at 1.1% detection while *Oprl1* holds at 35.57
CPM, and *Oprl1* leads 2.20-fold in all three samples (Figure 1R). The two ages come from one
laboratory, one platform, and one tissue, and the change between them is specific to *Oprk1*.
*Oprm1* reads 0.18 CPM at 0.2% detection at postnatal day 24, the lowest value for that gene in
any population measured here.

Depth constrains the margin without determining it. The two largest sensory margins come from
full-length libraries and the three smallest from droplet libraries, where the NodoMap median is
1,570 UMI per cell against 2.7 million reads in the spiral ganglion. Depth does not account for
the ordering within the droplet data: the superior cervical ganglion and the vestibular ganglion
sit at comparable depth and give margins of 25.72-fold and 4.29-fold, and the five ganglia of one
experiment span an 11-fold range in margin at one depth and one chemistry.

## Nuclear libraries scale with gene length and whole-cell libraries do not

Whole-cell and nuclear libraries of the same vagal tissue disagree about which opioid receptor is
highest. The four whole-cell deposits of the NodoMap atlas place *Oprl1* first at 11.13 CPM against
*Oprm1* 8.42; the 765 nuclei in the same atlas, from the same tissue and the same integration,
place *Oprm1* first at 244.49 against *Oprl1* 5.22 (Figure S1A). The dorsal root ganglion gives the
same disagreement in a second laboratory: 9.31 CPM against 8.23 in 31,802 whole cells, and 94.92
against 7.60 in 13,243 nuclei from six samples, with *Oprm1* first in each of the six (Figure S1B).
One of the two preparations misreports the abundance of these transcripts. Which one can be
decided without appeal to either result.

Genomic span decides it. The length of a transcription unit is a property of the locus and carries
no information about how much mature message a neuron holds, so a measurement that reports
abundance should not track it. We computed both preparations for every protein-coding gene in the
atlas and kept the 12,557 that reach 1 CPM in both, filtering on expression and never on length.
Whole-cell levels do not move with span: median expression runs between 19.4 and 29.6 CPM across
span deciles whose medians run from 4.0 kb to 253.0 kb, a 60-fold range of length, and the
correlation is r = -0.065. Nuclear levels rise monotonically over those same deciles, from 14.0 to
114.4 CPM, r = +0.43 (Figure S1E). The ratio between the two preparations therefore climbs with
length, from 0.43 in the shortest decile to 4.75 in the longest, crossing 1 at 21 kb (Figure S1F).
The nuclear measurement is the one a nuisance variable predicts.

The composition of nuclear RNA accounts for this. A nucleus holds nascent transcript that has not
been spliced or exported, and single-nucleus quantification counts reads across the whole gene
body, introns included, since restricting to exons discards most of a nuclear library. Signal then
accrues with the length of the transcription unit rather than with the number of finished
transcripts. Whole-cell libraries sample the cytoplasmic pool, which is the mature mRNA a neuron
translates, and 3' counting of that pool has no length term. The flat curve in Figure S1E is that
absence measured rather than assumed.

The four opioid receptors span 7 to 280 kb, which is why the distortion reaches them. *Oprl1*, at
7.1 kb, has a nuclear-to-whole-cell ratio of 0.45, the 35th percentile among the 1,982 genes within
a factor of 1.6 of its length. *Oprk1*, at 17.7 kb, sits at 1.14 and the 65th percentile of its own
length class. Both behave as ordinary genes of their size, and the ordering between them survives
the preparation. *Oprm1*, at 279.7 kb, has a ratio of 25.4 against a median of 4.45 for genes of
comparable length, which places it at the 96th percentile of that class. Length accounts for most
of its gain and not all of it, so the nuclear *Oprm1* value is inflated twice over and is not an
estimate of transcript abundance. *Oprd1* stays out of this comparison at 0.22 CPM in whole cells.
Detection follows the levels: *Oprm1* is called in 71.9% of the nuclei against 5.5% to 32.4% of
cells across the four whole-cell deposits, while *Oprl1* is called in 8.5% of the same nuclei,
inside the 5.3% to 27.5% whole-cell range (Figures S1C and S1D). Both genes were counted in the
same libraries, so sequencing depth does not separate them.

Every value in Figure 1 is therefore read from whole cells, and the ordering it reports is an
ordering of mature transcript. Mixing preparations moves the result part of the way toward the
nuclear answer: the trigeminal margin is 2.08-fold over whole cells and 1.18-fold once the nuclear
barcodes of the same atlas are pooled in, and in the NodoMap vagal ganglia 507 nuclei among 26,554
nodose barcodes take *Oprm1* from 9.86 to 14.63 CPM and put it first. Both iPain atlases are
majority single-nucleus, and every value drawn from them here comes from their whole-cell fraction
alone.

Preparation is confounded with laboratory in the two paired comparisons, since no laboratory here
ran both preparations on one tissue. The genome-wide result does not rest on that contrast: it is
internal to one atlas and one chemistry, and it identifies the biased measurement from the
behaviour of 12,557 genes rather than from any receptor. The same reasoning carries to human
tissue, where post-mortem and surgical ganglia are frozen and essentially all single-cell data from
them is single-nucleus. No such dataset can establish a species difference in the receptor ordering
while the mouse arm of the one cross-species nuclear experiment reverses in the same direction.
