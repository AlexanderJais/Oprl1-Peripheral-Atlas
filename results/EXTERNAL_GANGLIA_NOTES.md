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
