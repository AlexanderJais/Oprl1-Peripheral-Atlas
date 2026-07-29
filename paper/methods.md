# Methods

## Data

No new sequencing was generated. All data are public single-cell or single-nucleus RNA-seq
deposits, listed with their accessions, studies, platforms, sample sizes and quality-control
outcomes in Table S1. Nineteen neuronal populations from 16 mouse ganglia and plexuses were drawn
from 12 deposits published between 2016 and 2026, and were reached through GEO and CZ CELLxGENE.

Every value reported in Figure 1 comes from whole-cell libraries. Both iPain atlases and the
NodoMap atlas contain nuclear as well as whole-cell fractions; those were restricted to
`suspension_type == "cell"` before any receptor value was read, for the reason established in
Figure 2. The nuclear fractions are used only in Figure 2, where they are the subject.

## Quantification

Expression is per-cell counts per million, computed as the count of a gene in a cell divided by
that cell's total counts and multiplied by 10^6, and a population's level is the mean of that
quantity over all its neurons. GSE102443 was deposited as FPKM and is analysed in that unit;
levels are compared only within a population, never across units or platforms. GSE231924 was
deposited as Seurat log-normalised values at a scale of 10,000, from which CPM is recovered as
expm1(x) x 100, with the scale verified from the row sums rather than assumed.

A gene absent from a dataset's annotation is recorded as a null and never as a measured zero. Every
level table carries an `in_matrix` field derived from the matrix itself, so that a reference gap
and a measured absence remain distinguishable. The three zeros in the spiral ganglion are measured
absences by that test: all four receptors are present in the annotation and three are detected in
none of the 226 neurons.

## Cell filtering and neuron calling

Barcodes were retained at 2,000 or more UMI and 1,000 or more detected genes. Neurons were then
called on raw counts, so that the threshold does not move with library size, and glial and immune
barcodes were excluded on CPM, so that those thresholds do not move with sequencing run.

The neuronal criterion is *Snap25* at 5 counts or more together with a second marker chosen for the
tissue: *Tubb3* for the autonomic and vestibular ganglia, *Th* for the superior cervical ganglion,
*Phox2b* for the intrinsic cardiac nervous system, and *Elavl4* for the enteric plexus. Barcodes
carrying a glial, myelinating or immune signature were excluded whatever their *Snap25* count, at
*Sox10* below 200 CPM, *Plp1* below 500 CPM and *Ptprc* below 50 CPM. *Snap25* alone admits
ambient-dominated barcodes, which is what these exclusions remove. Populations deposited as sorted
or author-filtered neurons were taken as supplied.

## Marker gate

Before any receptor number was read from a population, its identity was confirmed against
*Snap25*, *Phox2b*, *Slc17a6*, *Tac1*, *Calca* and *Actb*. *Snap25* and *Actb* are enforced: they
must be present in the matrix and non-zero, and the pipeline raises rather than warns if they are
not. The remaining markers are tissue-specific, and are recorded without being enforced. All 21
populations in the quality panel pass.

## Ambient RNA

No ambient correction is applied anywhere in this work, so a transcript read in neurons carries
whatever the same transcript contributes to the surrounding soup. That contribution is measured
instead of removed: within a single dissociation, each receptor's mean CPM in the called neurons is
compared with its mean CPM in the non-neuronal barcodes captured beside them, as a log2 ratio with
a pseudocount of 0.01 CPM. *Plp1* and *Ptprc* serve as negative controls.

*Oprl1*'s enrichment on its own conflates how neuronal the transcript is with how cleanly the two
compartments separated in that dissociation. *Snap25* is neuronal by definition and calibrates the
second, so the quantity reported is *Oprl1*'s enrichment minus *Snap25*'s. Across the eleven
peripheral populations where both are available, *Oprl1* sits within 0.67 log2 of *Snap25*, median
-0.14. The check ran on 14 of 21 populations; a deposit of sorted or author-filtered neurons has no
non-neuronal compartment to score against, and the reason is recorded per population rather than
left as a gap.

## Receptor ordering, bootstrap and sample agreement

The four receptors are ranked within a population, on the same cells and the same chemistry, which
is what makes the ordering comparable between populations measured on different platforms. A rank
is reported as determinate only when the top receptor is strictly above the runner-up. The margin
is the ratio of the highest receptor to the second-highest.

Uncertainty in that ordering is estimated by resampling the cells of a population with replacement,
10,000 times, from a fixed seed. The 95% interval on the margin is the 2.5th and 97.5th percentile
of the resampled margins, and support is the fraction of resamples in which the observed top
receptor is still top. A pseudobulk ordering with a 1.8% margin and one with a 27% margin are
otherwise reported identically, which is what this number prevents.

Where a deposit resolves into biological samples, each sample was also ordered on its own cells,
and the count of samples placing the top receptor first is reported alongside the bootstrap.

## Whole-cell against nuclear preparation

Figure 2 compares the two preparations within the NodoMap integration: the same tissue, the same
chemistry and the same annotation, so the contrast is between preparations rather than between
studies. Preparation is nonetheless confounded with laboratory in the two paired tissue
comparisons, since no laboratory here ran both preparations on one tissue, and the genome-wide
analysis is what does not rest on that contrast.

Genomic spans are Ensembl GRCm39 gene loci, fetched once by symbol lookup so that every gene is
measured on one annotation. Per-dataset feature tables carry coordinates for the features they
quantify, but not for every gene and not on one assembly, which had previously dropped genes from
the comparison and misplaced others.

Genes were kept if they are protein-coding by Ensembl biotype and reach 0.1 CPM in both
preparations, which leaves 14,876. The floor exists only to keep ratios off a near-zero
denominator; selection is on expression and never on length, and it is low enough to admit all four
opioid receptors, which must not be filtered out of an analysis about them. Correlations are
Pearson on log10 span against log10 expression, with Spearman on the untransformed values reported
alongside. The floor was varied from none at all to 1 CPM: the whole-cell correlation runs from
+0.117 to -0.065 and the nuclear one from +0.399 to +0.433, so the choice of floor sharpens the
separation rather than creating it. Span deciles are ten equal-count bins of genes.

## Figures and tables

Figures are built to Cell Press artwork specification at the journal's fixed column widths, in
Helvetica or a metric-compatible substitute, and are written at the size they will print rather
than rescaled at submission. The build refuses to write a figure that exceeds the journal's width
or height or that places any element outside the canvas. Table S1 is assembled from the analysis
tables under `results/`; the build recomputes every margin from the levels Figure 1 is drawn from
and refuses to write the table if a recorded margin disagrees with it.

## Software and availability

Analysis in Python 3.11 with pandas, numpy, scipy, anndata, h5py and matplotlib. All analysis
code, the intermediate tables, and the scripts that produce every figure and table in this
manuscript are in the repository accompanying this paper; the accession for every dataset is in
Table S1.
