# `Oprl1` — the dominant opioid receptor of peripheral sensory neurons

The NOP receptor `Oprl1` is the highest-expressed opioid-family receptor in mouse peripheral sensory ganglia,
and within the vagus it marks the myelinated, Nav1.1+ afferents rather than the
unmyelinated nociceptors.

| input | tissue | source |
|---|---|---|
| [PNOC-Nodose](https://github.com/AlexanderJais/PNOC-Nodose) | nodose and jugular ganglia | NodoMap atlas, 106,436 cells in 52 clusters, 5 datasets |
| [Oprl1_Junhe](https://github.com/AlexanderJais/Oprl1_Junhe) | geniculate ganglion | GSE102443 (96 cells), GSE135801 (454 cells) |
| this repo | nucleus of the solitary tract | GSE166648, 49,392 neuronal nuclei in 25 subtypes |

Every number is recomputed from the GEO and CELLxGENE source matrices through one pipeline.
Data-quality analyses live in [`SUPPLEMENT.md`](SUPPLEMENT.md); the full code audit is in
[`AUDIT.md`](AUDIT.md).

---

## 1. `Oprl1` is the highest-expressed opioid receptor, in both ganglia

Mean expression of the four opioid-family receptors, measured in the same cells, in each
dataset's own unit (`results/*_opioid_levels.csv`):

| dataset | tissue | `Oprl1` | `Oprm1` | `Oprd1` | `Oprk1` |
|---|---|---|---|---|---|
| GSE102443 (FPKM) | geniculate | **5.73** | 0.18 | 0.31 | 0.83 |
| GSE135801 (CPM) | geniculate | **18.98** | 0.01 | 4.72 | 3.61 |
| NodoMap, nodose neurons (CPM) | nodose | **11.90** | 9.86 | 0.13 | 4.63 |
| GSE166648, NTS neurons (CPM) | NTS | 15.79 | 118.04 | 9.62 | 11.75 |

**In the geniculate ganglion `Oprl1` is roughly an order of magnitude above every other opioid
receptor** — 5.73 FPKM against `Oprk1` 0.83, `Oprd1` 0.31 and `Oprm1` 0.18 — and the same
pattern appears on an independent platform in an independent laboratory (GSE135801).

![Oprl1 in the geniculate ganglion](figures/figure1_geniculate_oprl1.png)

**In the nodose ganglion `Oprl1` is again the highest, at 11.9 CPM,** ahead of `Oprm1` (9.86),
`Oprk1` (4.63) and `Oprd1` (0.13), and it reproduces in each of the four whole-cell datasets
independently (11.98, 11.84, 11.09, 9.59 CPM).

![Oprl1 across the NodoMap atlas](figures/figure2_nodose_oprl1.png)

The nodose lead over `Oprm1` is real but narrow. Resampling the cells 2,000 times, `Oprl1`
stays highest in 100 %, 97 % and 94 % of resamples in Zhao, Bai and Kupari, but in only 54 % in
Buchanan, where the two differ by 1.8 % (`results/*_rank_support.csv`). Buchanan should be read
as "`Oprl1` and `Oprm1` are level", not as a fourth independent win. The geniculate datasets
carry the claim; the nodose datasets are consistent with it.

The eight-dataset version of this comparison is [`figureS3`](SUPPLEMENT.md).

> The two single-nucleus datasets put `Oprm1` first. That is a preparation artefact — nuclear
> preparations retain unspliced pre-mRNA and `Oprm1` spans 250 kb against `Oprl1`'s 6 kb — and
> is documented in [`SUPPLEMENT.md`](SUPPLEMENT.md). They carry no claim here.

## 2. `Oprl1` is expressed across the whole population, not by a subtype

**Geniculate.** Detected in 92 % of the 96 neurons in the deep dataset, at indistinguishable
levels in both divisions of the ganglion: gustatory (Phox2b+) 6.29 FPKM (n = 61),
somatosensory (Phox2b−) 4.77 FPKM (n = 35) — figure 1b. It is pan-geniculate, not a subtype
marker.

**Nodose and jugular.** Present in every one of the 26 neuronal clusters, from 27.7 % of cells
in NGN14 down to 0.8 % in NGN5, and near-absent from the non-neuronal compartment (log2
neuronal enrichment **+3.53**, comparable to `Oprk1` +3.49 and `Oprm1` +3.24, far above `Oprd1`
+1.62). The top clusters are NGN14 (43.5 CPM), NGN17 (31.7), NGN6 (30.5) and NGN19 (30.3) —
figure 2b.

## 3. `Oprl1` marks the myelinated, Nav1.1+ vagal afferents

The NodoMap authors annotated every nodose neuron by organ projection, fibre type, sensor type
and sodium-channel class. Each of those annotations is a property of the **cluster**, not of the
cell, so the honest unit of analysis is the 21 nodose clusters — testing across 26,047 cells
would treat correlated observations as independent and inflate every *p* value.

Tested that way (`results/nodose_oprl1_annotation_tests.csv`, Kruskal–Wallis over cluster means):

| annotation | mean CPM, high vs low | *p* (n = 21 clusters) |
|---|---|---|
| **sodium channel** | Nav1.1 **27.0** vs Nav1.8 **6.8** | **0.0010** |
| **fibre type** | Myelinated **23.1** vs Unmyelinated **4.5** | **0.021** |
| sensor type | Mechanosensor 18.7 vs Nocisensor 9.3 | 0.109 |
| organ projection | Pancreas 25.2 vs Jejunum/Ileum 5.3 | 0.336 |

**Two of the four hold.** `Oprl1` is a receptor of the myelinated, Nav1.1-expressing A-fibre
vagal afferents rather than of the unmyelinated, Nav1.8-expressing C-fibres. The sensor-type
and organ-projection splits point the same way and are the larger effects in raw CPM, but they
do not survive testing at cluster resolution — organ projection in particular has only one
cluster in four of its six levels, so it is descriptive only.

![Where Oprl1 sits](figures/figure3_oprl1_localisation.png)

**The same answer falls out of an unbiased transcriptome-wide scan.** Pseudobulk mean CPM was
computed for all 54,640 genes in each of the 21 clusters and every gene correlated with `Oprl1`
across them; 2,149 of 16,380 expressed genes reach FDR < 5 %
(`results/nodose_oprl1_gene_correlations.csv`). The top correlates are the transcriptional
signature of myelinated afferents — `Adgrg6`/*Gpr126*, the receptor required for Schwann-cell
myelination (rho = 0.88, *q* = 0.0003), with `Chgb` (0.90), `Ptn` (0.89), `Rph3a` (0.88),
`Cacng5` (0.86), `Atp1b1` (0.86) and `Lsamp` (0.84).

## 4. The satiation receptors are not among them

The hypothesis this project was built to test was that `Oprl1` sits on the `Glp1r` and `Cckar`
afferents, which would place a Gi-coupled brake on the first synapse of the gut–brain axis.
**In this atlas it does not, at the population level.** In the same unbiased ranking:

| gene | Spearman rho with `Oprl1` | *q* | rank among 16,380 genes |
|---|---|---|---|
| `Glp1r` | +0.13 | 0.74 | 8,461 |
| `Cckbr` | +0.11 | 0.79 | 8,828 |
| `Cckar` | −0.06 | 0.90 | 11,270 |

All three sit in the middle of the distribution. The `Glp1r`-high clusters NGN20 and NGN21 are
mid-range for `Oprl1`; the `Oprl1`-high cluster NGN14 is `Glp1r`-low.

What does survive is a weaker, cell-level association. Holding cluster identity and capture
depth fixed, a nodose neuron expressing `Glp1r` is 1.46× more likely to also express `Oprl1`
than its neighbours (Mantel–Haenszel, *p* = 6×10⁻⁵; `Cckar` 1.37×, *p* = 1×10⁻⁷). By raw
co-detection, 15.8 % of `Glp1r`+ and 12.8 % of `Cckar`+ nodose neurons carry `Oprl1` — floors,
given droplet dropout.

![Oprl1 and the satiation receptors](figures/figure4_vagal_oprl1_satiation.png)

**So the honest reading is that the brake is on the mechanosensory arm, not the peptide-receptor
arm.** A minority of GLP-1R and CCK-A afferents do carry NOP, but `Oprl1` is not a marker of
those populations; it is a marker of the myelinated mechanosensors that signal distension. If
N/OFQ acts on vagal afferents at all, the prediction from this data is that it damps
mechanosensory satiation signalling, and that is a testable and more specific claim than the one
the project started with.

## 5. Is the substrate for a mechanosensory brake actually there?

Section 3 yields a prediction: if N/OFQ acts on the vagus, it damps **mechanosensory** signalling
rather than peptide-receptor signalling. That prediction needs three things to be true of the
same neurons, and each is checkable here without a new experiment — the mechanotransducer
(`Piezo2`), the Gi effector machinery a NOP receptor would actually work through (GIRK channels
`Kcnj3/6/9`, `Gnai`/`Gnao`, `Cacna1b`), and the *absence* of the nociceptor programme (`Trpv1`,
`Trpa1`, `Scn10a`), which serves as the internal negative control.

![The mechanosensory brake substrate](figures/figure6_mechanosensory_brake.png)

**Supported.** `Piezo2`+ nodose neurons carry `Oprl1` far more often than their neighbours:
22.2 % against 6.9 %, and holding cluster identity and capture depth fixed the odds ratio is
**2.00** (*p* = 6×10⁻²²; `results/vagal_piezo2_oprl1_coexpression.csv`). That is the strongest
co-expression result in this project — a third again the `Glp1r` association. The nociceptor
programme runs the other way at the population level, exactly as a negative control should:
`Trpa1` rho = −0.68, `Trpv1` −0.66, `Scn10a` −0.63 (all *q* < 0.03).

**Not supported, and this matters.** Two parts of the prediction fail:

- `Piezo2` itself is **not** significantly correlated with `Oprl1` across clusters (rho = 0.34,
  *q* = 0.31). The per-cell association above is real but it is a within-cluster effect, not a
  population one.
- The Gi effector module **splits**. `Kcnj9` tracks `Oprl1` (rho = 0.68, *q* = 0.015) but the
  G-protein subunits run against it — `Gnai2` −0.75 (*q* = 0.006), `Gnao1` −0.65 (*q* = 0.022),
  `Cacna1b` −0.52. The claim "the Gi machinery co-localises with the receptor" is not supported
  by this data. `Gnai2` and `Gnao1` are abundant (346 and 520 CPM) and near-ubiquitous, so their
  cluster-level variation may track cell type rather than availability — but that is a
  hypothesis about why the test failed, not a rescue of it.

**And one result that inverts under control.** Crude `Trpv1`/`Oprl1` co-detection is strongly
*negative* (OR 0.37), but holding cluster and depth fixed it flips to slightly positive (OR
1.34, *p* = 5×10⁻⁶). The exclusion of `Oprl1` from nociceptors is a property of which clusters
are which, not of individual cells within a cluster.

**Where that leaves the prediction.** The receptor is on `Piezo2`-expressing, myelinated,
Nav1.1+ afferents and is largely absent from the nociceptor clusters, so the anatomy is
consistent with a brake on mechanosensory transmission. The effector arm is not demonstrated.
The experiment this points to is direct: apply N/OFQ to vagal afferents and measure the
mechanically evoked response — gastric or intestinal distension with recording from nodose
neurons, comparing wild-type with NOP knockout, or *ex vivo* with SB-612111. That is a specific,
falsifiable test, which is what this analysis was for.

## 6. The central relay

`Oprl1` is expressed across all 25 NTS neuronal subtypes, highest in the glutamatergic Glu9
(27.8 CPM, n = 1,155), Glu13 (22.6) and Glu7 (21.9). GSE166648 is a nuclear preparation, so its
receptor *ordering* is not usable (see the supplement) and no peripheral-against-central
comparison is made; the per-subtype `Oprl1` distribution is unaffected by that caveat.

![Oprl1 across NTS neuronal subtypes](figures/figure5_nts_oprl1.png)

---

## Methods

Levels are pseudobulk means per dataset: FPKM for GSE102443, mean per-cell CPM elsewhere.
Transcript rows are summed to gene level, and a gene symbol on more than one annotation row has
its rows summed rather than truncated to the first. Detection percentages are compared only
within one dataset, never across assays.

Statistics, in `src/atlas_common.py`:

- `receptor_rank()` — the ordering of the four receptors in one sample. Reports
  `determinate = False` when every level is zero or the top two are tied.
- `bootstrap_receptor_support()` — the fraction of 2,000 cell resamples in which the observed
  top receptor stays top, plus a 95 % interval on the margin.
- `stratified_odds_ratio()` — Mantel-Haenszel odds ratio holding capture depth, and separately
  depth and cluster identity, fixed. Co-detection in droplet data is confounded by depth: a
  cell detecting any gene tends to detect more genes overall, so a raw overlap percentage shows
  an association whether or not one exists.
- `transcriptome_percentile()`, `leave_one_out_pearson()` — supporting statistics.

Each dataset passes `check_markers()` before any `Oprl1` number is read from it: `Snap25` and
`Actb` are enforced, the tissue-specific markers are recorded. Preparation type is read from
the NodoMap atlas's own `suspension_type` field and checked against the registry in
`atlas_common.DATASETS`.

Figures follow the idiom of the sibling PNOC-Nodose project (`src/atlas_style.py`).

## Reproducing

```bash
pip install -r requirements.txt
export NODOSE_ROOT=/path/to/PNOC-Nodose      # defaults to /home/user/PNOC-Nodose
bash src/00_download_data.sh                 # GEO + the NodoMap atlas, checksummed
python3 src/01_geniculate.py                 # GSE102443 + GSE135801
python3 src/02_nodose.py                     # NodoMap atlas
python3 src/03_nts.py                        # GSE166648, streamed and cached
python3 src/04_synthesis.py                  # cross-tissue ranking + Figure S1
python3 src/05_vagal_coexpression.py         # Oprl1 x Glp1r/Cckar
python3 src/06_oprl1_localisation.py         # where Oprl1 sits: unbiased
python3 src/07_mechanosensory_brake.py       # is the brake substrate present
python3 -m pytest tests -q                   # 33 unit tests
```

| figure | contents |
|---|---|
| `figure1_geniculate_oprl1` | `Oprl1` per neuron by division, detection, opioid panel |
| `figure2_nodose_oprl1` | `Oprl1` on the published UMAP, across 26 neuronal clusters, by dataset |
| `figure3_oprl1_localisation` | `Oprl1` by organ, fibre type, sensor type, Nav class; top correlates |
| `figure4_vagal_oprl1_satiation` | `Oprl1` against `Glp1r`/`Cckar`: UMAP, per cluster, odds ratios |
| `figure5_nts_oprl1` | `Oprl1` in NTS neurons and across the 25 subtypes |
| `figure6_mechanosensory_brake` | `Piezo2`, the Gi effector module, and the nociceptor control |
| `figureS1`–`figureS3` | supplementary — see [`SUPPLEMENT.md`](SUPPLEMENT.md) |

| table | contents |
|---|---|
| `nodose_oprl1_by_cluster_annotated.csv` | per-cluster `Oprl1` with the atlas's annotations |
| `nodose_oprl1_annotation_tests.csv` | Kruskal-Wallis over cluster means, per annotation |
| `vagal_brake_module_correlations.csv` | `Oprl1` against the transduction, Gi and nociceptor modules |
| `vagal_piezo2_oprl1_coexpression.csv` | `Piezo2`/`Oprl1` co-detection, stratified |
| `nodose_oprl1_by_annotation.csv` | `Oprl1` by organ projection, fibre type, sensor type, Nav class |
| `nodose_oprl1_gene_correlations.csv` | every expressed gene correlated with `Oprl1` across clusters |
| `vagal_oprl1_coexpression.csv` | co-detection and stratified odds ratios per partner gene |
| `vagal_oprl1_satiation_by_cluster.csv` | `Oprl1` and the satiation panel per nodose cluster |
| `vagal_cluster_level_correlation.csv` | cluster-level Spearman against each partner |
| `oprl1_across_datasets.csv` | `Oprl1`'s rank among the four receptors, per dataset |
| `top_receptor_by_dataset.csv` | top receptor per dataset, with bootstrap support |
| `nodose_oprl1_by_cluster.csv` | `Oprl1` level and detection in each of the 52 clusters |
| `nts_by_subtype.csv` | `Oprl1` across the 25 NTS neuronal subtypes |
| `geniculate_per_cell_GSE102443.csv` | per-cell `Oprl1` FPKM, split gustatory/somatosensory |
| `*_opioid_levels.csv`, `*_receptor_rank.csv`, `*_rank_support.csv` | per-tissue levels, ordering, support |

## Provenance

- GSE102443 Dvoryanchikov et al. 2017, *Nat Commun*, 96 geniculate neurons, SMART-seq
- GSE135801 Zhang et al. 2019, *Cell* (Zuker lab), 454 Phox2b+ geniculate neurons
- GSE166648 Ludwig et al., dorsal vagal complex snRNA-seq, 72,128 nuclei
- NodoMap Cheng et al. 2026, *Cell Press Blue* 1:100072, doi:10.1016/j.cpblue.2026.100072,
  via CZ CELLxGENE collection `982f9f44-031c-4c8c-91ee-dcaa53b10151`
