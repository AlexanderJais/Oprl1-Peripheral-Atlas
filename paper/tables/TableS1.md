# Table S1

*Related to Figure 1. The legend is in [`paper/table_legends.md`](../table_legends.md); the table the journal receives is [`TableS1.csv`](TableS1.csv), of which this is a rendering.*

## (a) Provenance and sample

| Panel | Division                    | Population              | Deposit   | Study                      | DOI                          | Platform               | Preparation | Neurons | Samples | Median library |
|-------|-----------------------------|-------------------------|-----------|----------------------------|------------------------------|------------------------|-------------|---------|---------|----------------|
| A     | sensory                     | spiral (VIII)           | GSE114997 | Shrestha et al., 2018      | 10.1016/j.cell.2018.07.007   | SMART-seq, full-length | whole cell  | 226     | --      | 2,678,701      |
| B     | sensory                     | geniculate (VII)        | GSE102443 | Dvoryanchikov et al., 2017 | 10.1038/s41467-017-01095-1   | SMART-seq, full-length | whole cell  | 96      | --      | --             |
| C     | sensory                     | vestibular (VIII)       | GSE309608 | Liu et al., 2026           | 10.1073/pnas.2530677123      | 10x droplet, 3'        | whole cell  | 6,596   | 4       | 55,172         |
| D     | sensory                     | geniculate (VII)        | GSE135801 | Zhang et al., 2019         | 10.1016/j.cell.2019.08.031   | droplet, 3'            | whole cell  | 454     | --      | --             |
| E     | sensory                     | trigeminal (V)          | iPain     | Bhuiyan et al., 2024       | 10.1126/sciadv.adj9173       | 10x droplet, 3'        | whole cell  | 2,773   | --      | --             |
| F     | sensory                     | nodose (X)              | NodoMap   | Cheng et al., 2026         | 10.1016/j.cpblue.2026.100072 | 10x droplet, 3'        | whole cell  | 26,047  | --      | 1,570          |
| G     | sensory                     | jugular (X)             | NodoMap   | Cheng et al., 2026         | 10.1016/j.cpblue.2026.100072 | 10x droplet, 3'        | whole cell  | 4,593   | --      | 1,570          |
| H     | sensory                     | dorsal root             | iPain     | Bhuiyan et al., 2024       | 10.1126/sciadv.adj9173       | 10x droplet, 3'        | whole cell  | 31,802  | --      | --             |
| I     | sympathetic                 | celiac                  | GSE232789 | Sivori et al., 2024        | 10.7554/eLife.91576          | 10x droplet, 3'        | whole cell  | 247     | 1       | 56,023         |
| J     | sympathetic                 | pelvic                  | GSE232789 | Sivori et al., 2024        | 10.7554/eLife.91576          | 10x droplet, 3'        | whole cell  | 1,444   | 2       | 40,986         |
| K     | sympathetic                 | superior cervical       | GSE231766 | Ziegler et al., 2023       | 10.1126/science.abn6366      | 10x droplet, 3'        | whole cell  | 1,382   | 2       | 9,885          |
| L     | sympathetic                 | lumbar chain            | GSE232789 | Sivori et al., 2024        | 10.7554/eLife.91576          | 10x droplet, 3'        | whole cell  | 887     | 1       | 73,326         |
| M     | sympathetic                 | thoracic chain          | GSE78845  | Furlan et al., 2016        | 10.1038/nn.4376              | full-length            | whole cell  | 298     | --      | 33,099         |
| N     | sympathetic                 | stellate                | GSE232789 | Sivori et al., 2024        | 10.7554/eLife.91576          | 10x droplet, 3'        | whole cell  | 2,589   | 1       | 27,931         |
| O     | sympathetic                 | stellate                | GSE231924 | Sharma et al., 2023        | 10.7554/eLife.86295          | 10x droplet, 3'        | whole cell  | 1,303   | --      | --             |
| P     | parasympathetic and enteric | sphenopalatine (VII)    | GSE232789 | Sivori et al., 2024        | 10.7554/eLife.91576          | 10x droplet, 3'        | whole cell  | 2,014   | 1       | 37,416         |
| Q     | parasympathetic and enteric | intrinsic cardiac       | GSE330884 | Xu et al., 2026            | 10.1016/j.cell.2026.06.040   | 10x droplet, 3'        | whole cell  | 4,513   | 3       | 52,533         |
| R     | parasympathetic and enteric | enteric submucosal, P24 | GSE263422 | Li et al., 2025            | 10.1038/s41593-025-01962-x   | 10x droplet, 3'        | whole cell  | 7,787   | 3       | 12,874         |
| S     | parasympathetic and enteric | enteric submucosal, P7  | GSE263422 | Li et al., 2025            | 10.1038/s41593-025-01962-x   | 10x droplet, 3'        | whole cell  | 1,648   | 2       | 8,974          |

## (b) Levels, ordering and quality control

| Panel | Population              | Unit | Oprl1 | Oprm1 | Oprd1 | Oprk1  | Oprl1 detected (%) | Top   | Runner-up     | Margin    | 95% CI          | Support | Samples agreeing | Marker gate | Ambient check |
|-------|-------------------------|------|-------|-------|-------|--------|--------------------|-------|---------------|-----------|-----------------|---------|------------------|-------------|---------------|
| A     | spiral (VIII)           | CPM  | 41.94 | 0.00  | 0.00  | 0.00   | 74.8               | Oprl1 | none detected | undefined | --              | 1.000   | n/a              | pass        | not possible  |
| B     | geniculate (VII)        | FPKM | 5.73  | 0.18  | 0.31  | 0.83   | 91.7               | Oprl1 | Oprk1         | 6.92x     | 4.26 to 13.32   | 1.000   | n/a              | pass        | not possible  |
| C     | vestibular (VIII)       | CPM  | 19.88 | 1.78  | 0.09  | 4.63   | 57.6               | Oprl1 | Oprk1         | 4.29x     | 3.85 to 4.81    | 1.000   | 4/4              | pass        | run           |
| D     | geniculate (VII)        | CPM  | 18.98 | 0.01  | 4.72  | 3.61   | 19.6               | Oprl1 | Oprd1         | 4.02x     | 0.99 to 36.70   | 0.974   | n/a              | pass        | not possible  |
| E     | trigeminal (V)          | CPM  | 13.24 | 6.36  | 1.46  | 5.12   | --                 | Oprl1 | Oprm1         | 2.08x     | 1.41 to 3.26    | 1.000   | n/a              | pass        | not possible  |
| F     | nodose (X)              | CPM  | 11.90 | 9.86  | 0.13  | 4.63   | 9.1                | Oprl1 | Oprm1         | 1.21x     | 1.15 to 1.28    | 1.000   | 4/4              | pass        | run           |
| G     | jugular (X)             | CPM  | 9.62  | 8.33  | 0.92  | 0.90   | 7.2                | Oprl1 | Oprm1         | 1.16x     | 0.92 to 1.46    | 0.898   | 2/3              | pass        | run           |
| H     | dorsal root             | CPM  | 9.31  | 8.23  | 0.20  | 2.78   | --                 | Oprl1 | Oprm1         | 1.13x     | 1.08 to 1.19    | 1.000   | n/a              | pass        | not possible  |
| I     | celiac                  | CPM  | 25.03 | 0.42  | 0.03  | 0.56   | 66.4               | Oprl1 | Oprk1         | 44.41x    | 20.42 to 197.63 | 1.000   | 1/1              | pass        | run           |
| J     | pelvic                  | CPM  | 16.33 | 0.50  | 0.43  | 0.46   | 51.0               | Oprl1 | Oprm1         | 32.42x    | 22.00 to 54.20  | 1.000   | 2/2              | pass        | run           |
| K     | superior cervical       | CPM  | 18.68 | 0.27  | 0.57  | 0.73   | 21.6               | Oprl1 | Oprk1         | 25.72x    | 15.93 to 52.14  | 1.000   | 2/2              | pass        | run           |
| L     | lumbar chain            | CPM  | 24.76 | 1.68  | 0.03  | 0.39   | 76.5               | Oprl1 | Oprm1         | 14.78x    | 12.00 to 18.83  | 1.000   | 1/1              | pass        | run           |
| M     | thoracic chain          | CPM  | 62.30 | 1.48  | 3.90  | 4.24   | 82.2               | Oprl1 | Oprk1         | 14.69x    | 9.35 to 27.27   | 1.000   | n/a              | pass        | not possible  |
| N     | stellate                | CPM  | 27.55 | 0.56  | 0.03  | 4.53   | 52.5               | Oprl1 | Oprk1         | 6.08x     | 5.16 to 7.27    | 1.000   | 1/1              | pass        | run           |
| O     | stellate                | CPM  | 23.93 | 4.84  | 5.04  | 5.12   | 11.5               | Oprl1 | Oprk1         | 4.68x     | 3.02 to 8.08    | 1.000   | 8/8              | pass        | not possible  |
| P     | sphenopalatine (VII)    | CPM  | 19.50 | 0.29  | 0.04  | 0.01   | 49.9               | Oprl1 | Oprm1         | 66.86x    | 40.54 to 132.54 | 1.000   | 1/1              | pass        | run           |
| Q     | intrinsic cardiac       | CPM  | 11.13 | 3.04  | 0.36  | 0.00   | 38.2               | Oprl1 | Oprm1         | 3.65x     | 3.29 to 4.08    | 1.000   | 3/3              | pass        | run           |
| R     | enteric submucosal, P24 | CPM  | 35.57 | 0.18  | 16.14 | 3.77   | 32.4               | Oprl1 | Oprd1         | 2.20x     | 2.05 to 2.37    | 1.000   | 3/3              | pass        | run           |
| S     | enteric submucosal, P7  | CPM  | 26.70 | 0.74  | 2.29  | 122.51 | 23.4               | Oprk1 | Oprl1         | 4.59x     | 4.04 to 5.24    | 1.000   | 2/2              | pass        | run           |

## Why an ambient check could not run

- (A) deposit is 226 sorted neurons, no non-neuronal cells
- (B) deposit is 96 sorted neurons, no non-neuronal cells
- (D) deposit is Phox2b-sorted neurons, no non-neuronal cells
- (E) source h5ad is 0.6 GB and was not retained; not re-derived
- (G) scored on the pooled ganglion, not the jugular subset
- (H) source h5ad is 20.9 GB, above this session's disk allowance
- (M) deposit is 298 sorted neurons, no non-neuronal cells
- (O) deposit is author-filtered neurons, no non-neuronal cells
