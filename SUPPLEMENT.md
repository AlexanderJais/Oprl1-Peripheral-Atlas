# Supplementary material

Data-quality analyses. None of these are findings about `Oprl1` biology; they exist to
establish which datasets can carry which claims, and they are kept out of the main README so
the `Oprl1` result is not buried in method.

---

## Figure S1 — nuclear preparation inverts the opioid receptor ordering

![Nuclear preparation bias](figures/figureS1_nuclear_preparation_bias.png)

Six of the eight datasets in this atlas are whole-cell and two are single-nucleus. The two
nuclear datasets are also the only two that do not rank `Oprl1` first. That is a preparation
artefact, not a tissue difference, and this section is the evidence for treating it as one.

Within the nodose ganglion, where whole-cell and nuclear data exist for the *same tissue*, the
size of the shift tracks genomic span (`results/nuclear_bias_vs_gene_length.csv`):

| gene | genomic span | nuclear / whole-cell level |
|---|---|---|
| `Oprm1` | 250 kb | **29.0×** |
| `Oprd1` | 34 kb | **24.4×** |
| `Oprk1` | 18 kb | 0.95× |
| `Oprl1` | 6 kb | 0.47× |
| `Pomc` | 6 kb | 0.36× |
| `Penk` | 5 kb | 0.39× |
| `Pdyn` | 2 kb | 0.99× |

log-log Pearson *r* = 0.84 (n = 7 genes; Spearman rho = 0.58, *p* = 0.18). Nuclear preparations
retain unspliced pre-mRNA, so long-intron genes gain signal. `Oprm1` spans 250 kb against
`Oprl1`'s 6 kb and gains 29-fold, which is enough to move it from second to first.

Weighting the whole-cell baseline by cell count rather than by dataset changes little
(`Oprm1` 25.4×, `Oprl1` 0.45×).

**How much of that correlation is real.** It is carried by two of the seven genes. Dropping
`Oprm1` and `Oprd1` takes *r* from 0.84 to **0.02** (`results/nuclear_bias_sensitivity.csv`,
panel B) — the remaining five genes span 2–18 kb with ratios 0.36–0.99 and show no trend. Three
further caveats:

- `Oprd1`'s 24.4× is a ratio against a 0.22 CPM floor, where the estimate is unstable.
- Preparation is perfectly confounded with laboratory. There is one nuclear nodose dataset
  (765 neurons, in-house) and no in-house whole-cell dataset, so "nuclear" and "in-house"
  cannot be separated in this design.
- `Oprm1` detection rises from 5–9 % to 71.9 % in that dataset. A jump that large in the
  *fraction of cells with any read* fits an intron-inclusive alignment
  (`cellranger --include-introns`) as well as it fits pre-mRNA retention.

So the direction is clear and reproducible, the mechanism is plausible and consistent with the
two largest points, and the quantitative gene-length relationship across all seven genes is not
established.

**What this means for the main result.** Only whole-cell datasets are used for the receptor
ranking. The two nuclear datasets are reported for completeness and excluded from the claim.
The NTS dataset is nuclear, so this atlas makes no claim about a peripheral-against-central
receptor switch.

---

## Figure S2 — the full opioid panel across all 52 NodoMap clusters

![Opioid panel across NodoMap clusters](figures/figureS2_nodose_opioid_dotplot.png)

The complete eight-gene opioid panel across every cluster in the atlas, neuronal and
non-neuronal, in the same dot-plot idiom as the sibling PNOC-Nodose project. Dot size is the
percentage of cells expressing; colour is log(1 + mean CPM). The dashed line separates the
non-neuronal clusters (left) from the jugular and nodose neurons (right).

This is the context for the main-text statement that `Oprl1` is neuron-enriched: essentially
all of the opioid receptor signal in this ganglion sits to the right of that line. `Penk` is
the exception and is largely a fibroblast transcript (FB1–FB4, MGC5), which is why its neuronal
enrichment is negative.

---

## Marker-gene checks

Every dataset passes `check_markers()` before any `Oprl1` number is read from it. `Snap25` and
`Actb` are enforced — absent or zero raises `SanityCheckError` — and `Phox2b`, `Slc17a6`,
`Tac1` and `Calca` are recorded without being enforced, since they are tissue-specific.

| dataset | `Snap25` | `Phox2b` | `Slc17a6` | `Tac1` | `Calca` | `Actb` |
|---|---|---|---|---|---|---|
| GSE102443 (FPKM) | 2742.4 | 39.8 | 65.0 | 417.0 | 98.0 | 413.4 |
| GSE135801 (CPM) | 1707.6 | 777.0 | 295.6 | 453.2 | 35.1 | 1229.3 |
| GSE166648 neurons (CPM) | 890.2 | 30.9 | 83.2 | 29.1 | 17.1 | 242.2 |

Full tables in `results/*_marker_checks.csv`.

## Detection rate is not comparable across platforms

`Oprl1` is detected in 92 % of geniculate neurons on full-length SMART-seq and in 9 % of nodose
neurons on 10x droplet data. That difference is platform sensitivity, not biology, and no claim
in this repository compares detection rates across assays. Detection is compared only *within*
one dataset — between clusters of the same atlas, sequenced together — where it is a fair
comparison.

This is also why the co-expression analysis in the main text reports a depth-stratified odds
ratio rather than an overlap percentage: raw co-detection between two sparsely detected genes
is dominated by per-cell capture depth.

## Bootstrap support for the receptor ordering

`bootstrap_receptor_support()` resamples cells with replacement 2,000 times and reports how
often the observed top receptor stays top (`results/*_rank_support.csv`). This is in the main
README rather than here, because it qualifies the headline claim: `NodoMap:Buchanan` ranks
`Oprl1` first by 1.8 % and holds that ordering in only 54 % of resamples.

## Genes checked and not found

`Calcr` is not detected in a single nodose neuron in the whole-cell NodoMap data, and `Gfral`
(5 cells) and `Gipr` (19 cells) are effectively absent. They are reported in
`results/vagal_oprl1_coexpression.csv` for completeness and carry no statistic. `Pnoc` is not
present in the GSE102443 annotation at all, so it is recorded as unmeasured rather than as zero
everywhere it appears.
