# Ganglia added from external atlases

Results outside the numbered pipeline. The sensory sections below are exploratory; the autonomic
and enteric sections are reproduced by `src/external/autonomic_ganglia.py` and
`src/external/scg_GSE231766.py`, which fetch their own source files, apply the same marker gate as
the numbered pipeline and write `results/autonomic_receptor_levels.csv`,
`results/autonomic_marker_checks.csv`, `results/scg_GSE231766_receptor_levels.csv` and
`results/scg_marker_checks.csv`.

Source for the sensory sections: CZ CELLxGENE collection
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

## GSE231924, stellate ganglion, the third sympathetic dataset

Cardiac-projecting neurons of the mouse stellate ganglion, deposited already filtered to neurons
and Seurat log-normalised at a scale of 10,000. CPM is recovered exactly as expm1(x) * 100, and
the script verifies the scale from the row sums rather than assuming it. Identity: Dbh 1,134 CPM,
Prph 2,265, Th 616, Snap25 863.

| gene | mean CPM | % of 1,303 neurons detected |
|---|---|---|
| Oprl1 | 23.93 | 11.5 |
| Oprk1 | 5.12 | 2.3 |
| Oprd1 | 5.04 | 2.5 |
| Oprm1 | 4.84 | 2.7 |

Margin 4.68x over the runner-up (95% interval 3.02 to 8.08), 4.9x over Oprm1, support 1.0000. All
eight mice place Oprl1 first individually, from 15.24 to 55.19 CPM. This is the narrowest
sympathetic margin measured, and the only sympathetic dataset where the other three receptors sit
close together rather than at trace levels. The deposit carries no library sizes, so median UMI is
not recoverable and no depth comparison is made against it.

## GSE232789, six autonomic ganglia in one experiment

One laboratory, one platform, six dissections: stellate, coeliac and lumbar chain (sympathetic),
sphenopalatine (cranial parasympathetic) and two pelvic ganglia. This is the only dataset in the
project where the comparison between autonomic divisions carries no laboratory or chemistry
difference. Neurons were called on raw counts (Snap25 >= 5 and Tubb3 >= 5) with glial and immune
barcodes excluded (Sox10 < 200, Plp1 < 500, Ptprc < 50 CPM).

| ganglion | division | n | Oprl1 | Oprm1 | Oprd1 | Oprk1 | margin | support |
|---|---|---|---|---|---|---|---|---|
| stellate | sympathetic | 2,589 | 27.55 | 0.56 | 0.03 | 4.53 | 6.08x | 1.000 |
| coeliac | sympathetic | 247 | 25.03 | 0.42 | 0.03 | 0.56 | 44.41x | 1.000 |
| lumbar chain | sympathetic | 887 | 24.76 | 1.68 | 0.03 | 0.39 | 14.78x | 1.000 |
| sphenopalatine | parasympathetic | 2,014 | 19.50 | 0.29 | 0.04 | 0.02 | 66.86x | 1.000 |
| pelvic | mixed | 1,444 | 16.33 | 0.50 | 0.43 | 0.46 | 32.42x | 1.000 |

Oprl1 is first in all five, within a 1.7-fold band of each other (16.33 to 27.55 CPM), while the
margin over the runner-up ranges 11-fold. The margin is therefore set by how low the other three
receptors fall in a given ganglion, not by how high Oprl1 rises.

The sphenopalatine is cholinergic and not noradrenergic (Slc18a3 252 CPM, Chat 26, Th 2, Dbh 42)
against the sympathetic ganglia in the same dataset (Th 390 to 840, Dbh 803 to 1,176), so the
division assignment is confirmed within the data rather than taken from the sample titles. It
gives the largest lead over Oprm1 of any ganglion in this atlas, 67x.

The stellate row independently reproduces GSE231924: 27.55 against 23.93 CPM, two laboratories and
two platforms on the same ganglion.

## GSE330884, intrinsic cardiac nervous system

The parasympathetic ganglia of the heart, three samples, median library 40,499 UMI. Neurons called
on Snap25 >= 5 and Phox2b >= 2, giving 4,513 neurons at 68.5% of QC-passing cells, confirmed
cholinergic by Slc5a7 937 CPM and Slc18a3 475.

| gene | mean CPM | % of neurons detected |
|---|---|---|
| Oprl1 | 11.13 | 38.2 |
| Oprm1 | 3.05 | 12.0 |
| Oprd1 | 0.36 | 1.2 |
| Oprk1 | 0.00 | 0.0 |

Margin 3.65x over the runner-up (95% interval 3.29 to 4.08), support 1.0000, and the same ordering
in all three samples. This is the narrowest parasympathetic margin measured and the only
parasympathetic dataset in which Oprm1 is the runner-up rather than a trace transcript. Oprk1 is
detected in none of the 4,513 neurons at a median depth of 52,533 UMI, so that is a measured zero.

## GSE263422, enteric submucosal neurons at P7 and P24

Submucosal neurons of the mouse small intestine at two ages, one laboratory, one platform. Neurons
called on Snap25 >= 5 and Elavl4 >= 3 with the same glial and immune exclusions.

| age | n | median UMI | Oprl1 | Oprm1 | Oprd1 | Oprk1 | top | margin | support |
|---|---|---|---|---|---|---|---|---|---|
| P24 | 7,787 | 12,874 | 35.57 | 0.18 | 16.14 | 3.77 | Oprl1 | 2.20x | 1.000 |
| P7 | 1,648 | 8,974 | 26.70 | 0.74 | 2.29 | 122.51 | Oprk1 | 4.59x | 1.000 |

P7 is the only population in this atlas where Oprl1 is not the highest opioid receptor. Both P7
samples agree (Oprk1 119.86 and 125.36 CPM) and all three P24 samples agree (Oprl1 33.16 to
36.82). Oprk1 falls from 122.51 CPM at 32.7% detection to 3.77 CPM at 1.1% between the two ages
while Oprl1 stays between 26 and 36 CPM, so the change is specific to Oprk1 rather than a shift in
capture.

Oprm1 reads 0.18 CPM at 0.2% detection at P24. The literature on opioid-induced constipation
concerns Oprm1 in the myenteric plexus, and this dataset is submucosal, so the two do not sample
the same neurons and this is not a contradiction of it. No myenteric dataset located for this
survey passed neuronal QC: in GSE292613, the closest candidate, Snap25 is detected in 0.4% of
plexus barcodes and Ptprc in 46.5%, so the preparation is glia and immune cells.

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
