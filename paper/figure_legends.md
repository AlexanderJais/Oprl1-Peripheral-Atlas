# Figure legends

## Figure 1. *Oprl1* is the highest-expressed opioid receptor gene in peripheral neurons

Mean expression of the four opioid receptor genes in each of 19 peripheral neuronal populations,
grouped by division of the peripheral nervous system. Roman numerals give the cranial nerve of
ganglia that belong to one. Bars are pseudobulk means over all neurons of the population, in the
unit each dataset was measured in, ordered within each panel from highest to lowest. *Oprl1* is
blue and the other three receptors are gray. Sample size is the number of neurons contributing to
the panel.

(A) Spiral ganglion (VIII), GSE114997, full-length SMART-seq, median 2,678,701 reads per neuron.
*Oprl1* 41.94 CPM in 74.8% of neurons; *Oprm1*, *Oprd1*, and *Oprk1* detected in no cell. All four
genes are present in the annotation, so the three zeros are measured absences.

(B and D) Geniculate ganglion (VII) on two platforms from two laboratories. (B) GSE102443, full-
length SMART-seq, FPKM, 31.30-fold over the second receptor, bootstrap support 1.000. (D)
GSE135801, 3' droplet, 4.02-fold, support 0.974 (95% interval 0.99 to 36.70).

(C) Vestibular ganglion (VIII), GSE309608, four mice, 4.29-fold (3.85 to 4.81), support 1.000,
same ordering in each mouse.

(E) Trigeminal ganglion (V), iPain atlas, whole-cell neurons only, 2.08-fold (1.41 to 3.26).

(F and G) Inferior (nodose) and superior (jugular) ganglia of the vagus (X), NodoMap atlas, whole-
cell datasets only. Nodose 1.21-fold (1.15 to 1.28), support 1.000, first in all four contributing
datasets; jugular 1.16-fold (0.92 to 1.46), support 0.898, first in two of three.

(H) Dorsal root ganglion, a spinal ganglion with no cranial nerve, iPain atlas, whole-cell neurons
only, 1.13-fold (1.08 to 1.19), the narrowest margin measured and the second largest sample.

(I to O) Sympathetic ganglia, ordered by margin. Celiac 44.41-fold (20.42 to 197.63), pelvic
32.42-fold (22.00 to 54.20), superior cervical 25.72-fold (15.93 to 52.14), lumbar chain
14.78-fold (12.00 to 18.83), thoracic chain 14.69-fold (9.35 to 27.27), stellate 6.08-fold (5.16
to 7.27), stellate 4.68-fold (3.02 to 8.08). Support 1.000 throughout. (I), (J), (L), and (N) are
from GSE232789; (K) from the two untreated animals of GSE231766; (M) from GSE78845; (O) from
GSE231924, which profiles cardiac-projecting neurons and places *Oprl1* first in each of eight
mice. The stellate ganglion appears twice, in (N) and (O), from two deposits, two laboratories,
and two neuron-calling routes. The pelvic ganglion carries both sympathetic and parasympathetic
neurons and is grouped here with the sympathetic block.

(P and Q) Parasympathetic ganglia. (P) Sphenopalatine ganglion (VII), also called the
pterygopalatine ganglion, GSE232789, 66.86-fold (40.54 to 132.54), the largest margin in the
study. These neurons are cholinergic and not noradrenergic (*Slc18a3* 252 CPM, *Chat* 26, *Th* 2,
*Dbh* 42) against *Th* 390 to 840 CPM in the sympathetic ganglia of the same experiment. (Q)
Intrinsic cardiac nervous system, GSE330884, three mice, 3.65-fold (3.29 to 4.08), the only
population in this division in which *Oprm1* is the second receptor. These ganglia are intrinsic
to the heart and receive vagal preganglionic input, so they carry no cranial nerve numeral.

(R and S) Enteric submucosal neurons of the small intestine at two ages, GSE263422, one laboratory
and one platform. (R) Postnatal day 24, *Oprl1* first at 2.20-fold (2.05 to 2.37) in all three
samples. (S) Postnatal day 7, *Oprk1* first at 122.51 CPM against *Oprl1* 26.70, 4.59-fold (4.04
to 5.24) in both samples. *Oprk1* falls to 3.77 CPM at 1.1% detection by day 24 while *Oprl1*
holds between 26 and 36 CPM.

All panels use whole-cell data. Both iPain atlases are majority single-nucleus and are restricted
here to `suspension_type == "cell"`; pooling preparations compresses the trigeminal margin from
2.08-fold to 1.18-fold, for the reason given in Figure S1. Neurons were called on raw counts with
glial and immune barcodes excluded, and every population passed the same marker gate and ambient-
RNA check before any receptor value was read from it. Margins are the ratio of the highest to the
second-highest receptor; intervals in parentheses are 95% bootstrap intervals over 10,000
resamples of the cells, and support is the fraction of those resamples retaining the observed top
receptor. Panels (I), (J), (L), (N), and (P) come from a single experiment in which six autonomic
ganglia were dissected and sequenced together, so laboratory, platform, and sequencing depth are
constant across those five.

See also Figure S1, Figure S3, and Table S1.

## Figure S1. A nuclear preparation reverses the receptor ordering, related to Figure 1

Whole-cell and single-nucleus measurements of the same tissue, and the gene property that separates
them. Solid bars are whole-cell data and hatched bars are nuclear data throughout. *Oprl1* is blue
where the panel distinguishes genes by colour, as in Figure 1.

Axes in (A) and (B) are linear, so the size of the nuclear gain reads directly off them. At that
scale the whole-cell ordering the panels are compared against, *Oprl1* 11.13 CPM over *Oprm1* 8.42
in the vagal ganglia and 9.31 over 8.23 in the dorsal root ganglion, is smaller than the line
weight; the values are given below and the ordering itself is Figure 1.

(A) Vagal ganglia (X), NodoMap. The whole-cell value is the unweighted mean over the atlas's four
whole-cell deposits, 30,640 neurons; the nuclear value is the atlas's own 765 nuclei, the same
tissue and the same integration. Whole cell: *Oprl1* 11.13 CPM, *Oprm1* 8.42, *Oprk1* 4.85,
*Oprd1* 0.22. Nuclear: *Oprm1* 244.49, *Oprd1* 5.38, *Oprl1* 5.22, *Oprk1* 4.63. *Oprm1* rises
29-fold and *Oprl1* falls by half, which moves *Oprm1* from second to first.

(B) Dorsal root ganglion, a second tissue and a second laboratory. Whole cell, iPain atlas, 31,802
neurons: *Oprl1* 9.31 CPM, *Oprm1* 8.23, *Oprk1* 2.78, *Oprd1* 0.20. Nuclear, GSE201654 mouse arm,
six samples and 13,243 nuclei: *Oprm1* 94.92, *Oprl1* 7.60, *Oprd1* 2.33, *Oprk1* 2.17. *Oprm1* is
first in each of the six nuclear samples.

(C) Detection rate of *Oprm1* divided by detection rate of *Oprl1*, within each of the five deposits
that make up the NodoMap integration, so per-dataset sequencing depth cancels. The four whole-cell
deposits give 1.18 (Bai), 1.04 (Buchanan), 0.54 (Kupari) and 0.92 (Zhao); the nuclear deposit gives
8.46. Deposits are named for the authors of the data. Bai, Buchanan, Kupari and Zhao are published
cell suspensions that the NodoMap authors re-analysed; Cheng is the deposit the atlas labels
in-house, the 765 nuclei the NodoMap authors generated themselves, and is the nuclear arm of (A) and
(D) as well. *Oprl1* is detected in 8.5% of those nuclei against 7.7% of Zhao's cells, so the
nuclear libraries are not simply deeper.

(D) Ratio of nuclear to whole-cell level in (A) against the genomic span of the gene, both axes
logarithmic, for all eight opioid genes. The dashed line at 1 is no change, and the axis is
logarithmic so that a fall and a rise of the same factor sit the same distance from it. *Oprm1* spans 280 kb and gains 29.05-fold, *Oprd1* 34 kb
and 24.44-fold, *Pnoc* 25 kb and 1.52-fold, *Oprk1* 18 kb and 0.95-fold, *Pdyn* 14 kb and
0.99-fold, *Penk* 9 kb and 0.39-fold, *Oprl1* 7 kb and 0.47-fold, *Pomc* 6 kb and 0.36-fold.
Log-log Pearson r = 0.88, Spearman rho = 0.95, p = 0.0004 over the eight genes. Dropping *Oprm1*
and *Oprd1*, the two genes that gain most, leaves r = 0.95. Dashed line marks no change. Spans are
Ensembl GRCm39 gene loci (`data/raw/ensembl_gene_spans.csv`); taking them instead from a deposit's
own coordinate columns, as an earlier version of this analysis did, omits *Pnoc*, which GSE102443
does not quantify.

Nuclear libraries retain unspliced pre-mRNA, so a gene's signal scales with how much intron it
carries. This is a bias in a known direction rather than a failure of the assay, and it is fatal to
this comparison because the four receptors span 7 to 280 kb. It is not confined to the receptors:
the ordering by span holds across all eight opioid genes, including three that are not receptors.

Preparation is confounded with laboratory in both (A) and (B), since no laboratory here has run
both preparations on one tissue, and an intron-inclusive alignment would produce the same signature
as pre-mRNA retention. What separates those readings from a laboratory effect is that the shift is
ordered by gene length (D) and reproduces in two tissues across two independent pairs of
laboratories (A and B).

Every panel of Figure 1 therefore uses whole-cell data, and no human single-nucleus dataset can
establish a species difference in the ordering: the mouse arm of the one cross-species nuclear
experiment, shown in (B), reverses in the same direction.
