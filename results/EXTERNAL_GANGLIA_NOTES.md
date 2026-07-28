# Ganglia added from external atlases

Exploratory results not yet wired into the pipeline. Source: CZ CELLxGENE collection
`03608e22-227a-4492-910b-3cb3f16f952e` (iPain Atlas), mouse.

Both atlases are majority single-nucleus (TG 70,772 of 84,658; DRG 123,645 of 191,798), so
every number below is restricted to `suspension_type == "cell"` for the reason in figure S1.
Pooling all preparations compresses the trigeminal margin from 2.08x to 1.18x.

## Receptor levels, whole-cell neurons, mean CPM

| ganglion | n | Oprl1 | Oprm1 | Oprd1 | Oprk1 | margin | 95% CI | support |
|---|---|---|---|---|---|---|---|---|
| trigeminal | 2,773 | 13.24 | 6.36 | 1.46 | 5.12 | 2.08x | 1.41-3.26 | 1.000 |
| dorsal root | 31,802 | 9.31 | 8.23 | 0.20 | 2.78 | 1.13x | 1.08-1.19 | 1.000 |

## Spiral ganglion (VIII), GSE114997

Shrestha et al. 2018, 226 spiral ganglion neurons, full-length SMART-seq, median library
2,678,701 reads. All four opioid receptors are present in the annotation.

| gene | mean CPM | % of cells detected |
|---|---|---|
| Oprl1 | 41.94 | 74.8 |
| Oprm1 | 0.00 | 0.0 |
| Oprd1 | 0.00 | 0.0 |
| Oprk1 | 0.00 | 0.0 |

Oprl1 is the only opioid receptor detected in any of the 226 neurons. The other three are
quantified and read exactly zero, so this is a measured absence rather than a reference gap.
Library depth is 730-fold above the nodose droplet median, which is what makes a zero
interpretable here.

The same cells carry Scn1a at 310.2 CPM, Pvalb at 882.6 and Scn10a at 0.00: the Nav1.1-positive,
Nav1.8-negative profile that section 3 predicts should be Oprl1-high.

## Vestibular ganglion (VIII), GSE309608

Four mouse vestibular ganglia, 10x raw feature matrices. Barcodes were filtered at 2,000 UMI and
1,000 detected genes, leaving 25,431 cells at a median library of 6,000 to 10,300. Neurons were
called on raw counts (Snap25 >= 5 and Tubb3 >= 5) with glial and immune barcodes excluded
(Sox10 < 200, Plp1 < 500, Ptprc < 50 CPM), giving 6,596 neurons, 25.9% of QC-passing cells, at
Snap25 3,139 CPM and Plp1 163 CPM.

| gene | mean CPM | % of neurons detected |
|---|---|---|
| Oprl1 | 19.88 | 57.6 |
| Oprk1 | 4.63 | 9.7 |
| Oprm1 | 1.78 | 6.7 |
| Oprd1 | 0.09 | 0.5 |

Margin 4.29x over Oprk1 and 11.2x over Oprm1, bootstrap support 1.000. The same ordering appears
in each of the four samples independently. Scn1a reads 398.0 CPM, Scn10a 0.00 and Pvalb 4,101,
the Nav1.1-positive Nav1.8-negative profile also seen in the spiral ganglion.

A first pass using a 500 UMI cutoff and a Snap25 CPM threshold called 88.6% of barcodes neurons
and gave Oprl1 10.30 CPM. That filter was admitting ambient-dominated barcodes; the numbers above
use the stricter one.

## Sympathetic ganglion, GSE78845 and GSE231766: the specificity control

Two datasets, chosen to share nothing but the cell type: different ganglion, different laboratory,
different platform, different chemistry.

## GSE78845, thoracic chain, full-length

Furlan et al. 2016, 298 mouse thoracic sympathetic neurons, full-length, median library 33,099.
Identity confirmed by Th 2,316 CPM, Dbh 1,822, Prph 1,222, Snap25 1,995.

| gene | mean CPM | % of cells detected |
|---|---|---|
| Oprl1 | 62.30 | 82.2 |
| Oprk1 | 4.24 | 7.0 |
| Oprd1 | 3.90 | 10.4 |
| Oprm1 | 1.48 | 4.4 |

Margin 14.69x over the runner-up (95% interval 9.35 to 27.27), 42x over Oprm1, support 1.0000.
Oprl1 sits at the 81.2nd percentile of expressed genes here, above the nodose (73.6) and
geniculate (63.8) figures.

This control was run to test whether the ordering is specific to sensory neurons. It is not.
Sympathetic neurons are not sensory, are neural-crest-derived, and show the pattern at a larger
margin than any sensory ganglion except the geniculate and the spiral. The restriction to
"sensory" in earlier drafts of the title was an assumption that had never been tested, and the
test does not support it.

The caveat that stood against this result was that it rested on 298 cells in the platform class
that produces the largest margins throughout this project. GSE231766 addresses both halves of it.

## GSE231766, superior cervical ganglion, droplet

Reproduced by `src/external/scg_GSE231766.py`, which fetches the GEO tarball and writes
`results/scg_GSE231766_receptor_levels.csv` and `results/scg_marker_checks.csv`. It is the one
analysis in this file with a checked-in script.

Ziegler et al. 2023 (GSE231766), four 10x samples of the mouse superior cervical ganglion: two
untreated animals and two with heart disease at 5 and 18 days after transverse aortic
constriction. Different ganglion, laboratory and platform from GSE78845, and a cranial rather
than a thoracic sympathetic ganglion.

Barcodes were filtered at 2,000 UMI and 1,000 detected genes, leaving 25,199 cells. Sympathetic
neurons were called on raw counts (Snap25 >= 5 and Th >= 5) with glial and immune barcodes
excluded (Sox10 < 200, Plp1 < 500, Ptprc < 50 CPM), giving 2,115 neurons, 8.4% of QC-passing
cells, at Th 1,010 CPM, Dbh 2,285, Prph 2,424, Snap25 1,106, Plp1 98 and Sox10 7.

The headline row uses the 1,382 neurons from the two untreated animals, at a median library of
9,885 UMI:

| gene | mean CPM | % of neurons detected |
|---|---|---|
| Oprl1 | 18.68 | 21.6 |
| Oprk1 | 0.73 | 1.4 |
| Oprd1 | 0.57 | 1.2 |
| Oprm1 | 0.27 | 0.4 |

Margin 25.72x over the runner-up (95% interval 15.93 to 52.14), 69x over Oprm1, support 1.0000.
Oprl1 sits at the 75.4th percentile of expressed genes. All 2,115 neurons together give
Oprl1 21.86, Oprk1 1.26, Oprd1 1.04, Oprm1 0.31, a margin of 17.35x with support 1.0000, and each
of the four samples places Oprl1 first on its own: 14.18, 34.62, 31.53 and 22.64 CPM against
runners-up of 0.43, 1.79, 2.11 and 2.51. The two disease samples therefore neither produce the
result nor obscure it.

Scn1a reads 16.12 CPM against Scn10a 0.09, and Pvalb is absent. Sympathetic neurons carry the
Nav1.1-positive Nav1.8-negative profile without being proprioceptive, so the section 3 association
does not require the sensory context in which it was found.

The two sympathetic datasets agree on direction and on the size of the lead over Oprm1 (42x and
69x), across a 20-fold difference in library depth. That is the strongest cross-platform agreement
in this project, and it comes from the tissue that was introduced to break the pattern.

## Petrosal ganglion (IX)

No public dataset isolates the petrosal. A GEO search across all assay types and organisms returns
one mouse record containing petrosal tissue, GSE145216, whose samples are titled "NJP ganglia":
nodose, jugular and petrosal dissected and sequenced together. Petrosal neurons cannot be
separated from that pool without a marker that distinguishes them, and none is established. The
petrosal is therefore absent from this survey for a reason that is unlikely to change without a
dedicated dissection.

## The proprioceptor prediction

Stated in advance: if Oprl1 tracks Nav1.1 as a general property of sensory neurons, DRG
proprioceptors (Pvalb+/Runx3+/Ntrk3+, the most Nav1.1-dependent sensory population) should sit
at or near the top of Oprl1 expression.

NF2 is the proprioceptor population in this atlas and ranks 1 of 9 subtypes for Pvalb (1576
CPM), Runx3 (86.9) and Scn1a (226.4). It also ranks 1 of 9 for Oprl1 (29.97 CPM). Oprm1 ranks
7 of 9 in the same population.

Across the 9 whole-cell DRG neuronal subtypes (`drg_oprl1_by_subtype.csv`):

| partner | Spearman rho with Oprl1 | p |
|---|---|---|
| Scn1a | +0.862 | 0.0028 |
| Runx3 | +0.837 | 0.0049 |
| Pvalb | +0.817 | 0.0072 |
| Ntrk3 | +0.550 | 0.125 |
| Scn10a | -0.700 | 0.0358 |
| Oprm1 | -0.183 | 0.637 |

Both directions replicate the nodose result: positive with Scn1a, negative with Scn10a.

Oprm1 is highest in PEP1 (13.93) and SST (12.06), the peptidergic nociceptor populations where
the peripheral analgesia literature places it. The pipeline reproduces that distribution while
showing Oprl1 higher overall.
