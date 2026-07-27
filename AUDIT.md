# Code audit — Oprl1-Peripheral-Atlas

Full read of `src/*.py`, `src/00_download_data.sh`, `README.md`, and all 16 tables in
`results/`, at commit `cbc374f`. Every numeric claim in the README was recomputed from the
committed CSVs; every finding below was reproduced rather than inferred.

No analysis code was changed. Several findings would move published numbers, so the fixes
are described but left for the author to approve.

**Verified correct.** The README's arithmetic is sound. All seven `nuclear_over_whole_cell`
ratios, all six `Oprl1`:`Pnoc` ratios, the log-log Pearson *r* (0.8356 → "0.84"), the Spearman
rho and *p* (0.577 / 0.175 → "0.58 / 0.18"), the 4.4× NTS:nodose `Pnoc` ratio, the 96/454 cell
counts, the 91.67 % and 9.1 % detection rates, and the 49,392 neuronal nuclei all reproduce
exactly from the committed tables. The overall design — comparing only within-sample ranks,
ratios and percentiles across assays — is the right call for this data, and the README is
unusually candid about its own weak claims.

---

## Severity 1 — broken

### 1.1 The GSE166648 cache path always raises `TypeError`

`src/03_nts.py:39`

```python
if set(TARGETS).issubset(set(cached) | z["absent"].tolist()):
```

`set | list` is not a supported operation. Confirmed:

```
TypeError: unsupported operand type(s) for |: 'set' and 'list'
```

The cache is written on the first run, and every subsequent run crashes on load — before the
fallback re-stream is ever reached. The README advertises this cache as a headline feature
("streamed once … then cached", "cached so later runs are cheap"); in practice `03_nts.py`
runs exactly once and then cannot be run again without deleting `data/gse166648_targets.npz`.

Fix: `set(cached) | set(z["absent"].tolist())`. Same line, `allow_pickle=True` at `03_nts.py:36`
is unnecessary — every stored array is a plain string or float array and `np.load` reads the
file without it. Dropping it removes an arbitrary-code-execution surface on a cache file.

### 1.2 A gene missing from the NodoMap annotation is recorded as a measured zero

`src/02_nodose.py:95` and `src/02_nodose.py:128`

`load_counts()` skips genes absent from the atlas with a `[warn]` print (`02_nodose.py:38-39`).
Downstream, `m.get(g, 0.0)` returns `0.0` for those genes — and the row is then stamped:

```python
"in_matrix": True,          # 02_nodose.py:95, hardcoded
"Pnoc_in_matrix": True,     # 02_nodose.py:128, hardcoded
```

So a reference gap is written to `nodose_opioid_levels.csv` as `mean_level = 0.0,
in_matrix = True` — indistinguishable from a gene measured at zero.

This is precisely the error the project exists to correct. `atlas_common.ligand_receptor_ratio`
documents it explicitly:

> Absence from the quantified annotation and a measured zero are treated the same way here …
> They are not the same evidence, so callers report `<ligand>_in_matrix` alongside this value.

`01_geniculate.py:79` computes `in_matrix` honestly (`bool(g in mat.index)`), which is how
`Pnoc`'s absence from GSE102443 was caught in the first place — the finding the README's
"Corrections" section is built on. `02_nodose.py` hardcodes the same flag to `True`, so the
identical bug in the nodose data would be invisible. It currently reports no warnings, so no
published number is wrong today; the guard simply is not there.

Fix: have `load_counts` return the set of found symbols and derive both flags from it.

---

## Severity 2 — results do not match the code that claims to produce them

### 2.1 `Oprl1_over_Pnoc` is `NaN` in the committed CSVs, but the code returns `inf`

`results/geniculate_ligand_receptor.csv` and `results/ligand_receptor_by_tissue.csv` both carry
an **empty** field for GSE102443. `atlas_common.ligand_receptor_ratio` (`atlas_common.py:104`)
cannot produce that value: with `Oprl1 = 5.7329` finite and `Pnoc` absent (`levels.get("Pnoc",
0.0)` → `0.0`, verified), it returns `np.inf`, and pandas serialises `np.inf` as the literal
string `inf`, not as blank. Both verified against pandas 3.0.5.

The committed results were therefore produced by a different version of this function than the
one in the tree. Whichever behaviour is intended, the repo currently cannot reproduce its own
outputs — which for a project whose opening claim is "Every number here is recomputed from the
GEO source matrices through one pipeline" is the most consequential finding in this audit.

Re-running the pipeline end to end and committing the regenerated tables would settle it. `inf`
is arguably the better value (it is what the docstring promises, and `04_synthesis.py:166`
already filters on `np.isfinite`), but `NaN` avoids `inf` leaking into a log-scale axis.

### 2.2 The README's sanity-gene gate is not implemented

`README.md:141` states:

> Each dataset is checked against `Snap25`, `Phox2b`, `Slc17a6`, `Tac1`, `Calca` and `Actb`
> before any opioid number is read from it.

In fact:

- `01_geniculate.py:63-66` **prints** the six levels. It does not assert or gate anything; the
  script proceeds identically whichever values come back.
- `02_nodose.py:59` loads `SANITY_GENES` into `counts` and then never reads them again. No
  print, no check.
- `03_nts.py:31` puts them in `TARGETS`; only `Slc17a6` and `Gad1` are used, for the subtype
  table. `Snap25`, `Phox2b`, `Tac1`, `Calca`, `Actb`, `Slc17a7`, `Glp1r`, `Calcr` are streamed,
  cached, and discarded.

A stated QC control that exists in one of three pipelines, as a print rather than a check,
should either be implemented or removed from Methods.

### 2.3 Three README numbers are not produced by any script

- **`README.md:103`, "97.6th percentile of expressed genes"** for `Penk` in GSE102443.
  `transcriptome_percentile()` is called exactly once in the codebase
  (`01_geniculate.py:120-121`) and only for `Oprl1`. The `Penk` percentile appears in no table.
  The two percentiles that *are* computed — `Oprl1` at 63.8 (GSE102443) and 65.7 (GSE135801) —
  are written to CSV and never mentioned in the README.
- **`README.md:117`, "GSE102443 quantifies 17,111 genes"** — printed at `01_geniculate.py:40`,
  never persisted.
- **`README.md:181`, "Cheng et al. 2026, *Cell Press Blue* 1:100072"** — "Cell Press Blue" is
  not a journal title. Needs checking against the actual citation.

---

## Severity 3 — claims that overstate what the data supports

### 3.1 The gene-length mechanism rests on two of seven points

This is the paper's mechanistic claim, the bold sentence at `README.md:61`, and panel B's title
("The inversion **is** a gene-length effect"). Recomputed from
`results/nuclear_bias_vs_gene_length.csv`:

| points | *n* | log-log Pearson *r* |
|---|---|---|
| all | 7 | **0.836** |
| drop `Oprm1` | 6 | 0.695 |
| **drop `Oprm1` + `Oprd1`** | **5** | **0.016** |

Remove the two genes that gain signal and the relationship is *gone* — not weakened, absent
(*r* = 0.016). The remaining five genes span 2–18 kb with ratios 0.36–0.99 and show no trend.
The correlation is two points at the top right and a flat cloud at the bottom left.

The README discloses the weak Spearman (rho = 0.58, *p* = 0.18) and calls it "what seven points
can support" — genuinely good practice — but then states the conclusion in bold as settled and
titles the figure panel with it. A leave-one-out or leave-two-out line belongs next to the
Spearman.

Compounding this:

- **`Oprd1`'s 24.4× is a ratio against a floor value.** Whole-cell `Oprd1` is 0.22 CPM at
  0.22–1.26 % detection. A 24× ratio off that denominator is not a stable estimate, and it is
  one of the two points carrying the entire correlation.
- **`Pomc` is almost certainly ambient.** It contributes 11.5 CPM in nodose neurons at 11.5 %
  detection, with `log2_enrichment` = 0.305 over non-neurons (`nodose_neuron_enrichment.csv`) —
  i.e. essentially no neuronal enrichment, the signature of ambient RNA. `Pomc` is a
  hypothalamic/pituitary transcript. It is used as one of the seven regression points.

### 3.2 Preparation type is perfectly confounded with laboratory

The entire nuclear-artefact argument compares four whole-cell NodoMap datasets against
**one** nuclear dataset, `NodoMap:inhouse`, **765 cells** (2.4 % of the 31,405 neurons). "Nuclear"
and "in-house lab" are the same variable in this design — there is no second nuclear nodose
dataset and no in-house whole-cell dataset to separate them. Batch, dissociation protocol,
reference build and alignment settings are all equally consistent with the observation.

A specific alternative the README does not address: `Oprm1` detection jumps from 5–9 % in the
whole-cell datasets to **71.9 %** in `inhouse`, alongside the 244 CPM level. A >7× jump in the
*fraction of cells with any read* is more consistent with an intron-inclusive alignment
(`cellranger --include-introns`) or a different reference than with pre-mRNA retention alone.
That would produce the same gene-length correlation through a different mechanism, and it is
checkable from the NodoMap metadata.

`README.md:44` frames this as "Within the nodose ganglion, where whole-cell and nuclear data
exist for the *same tissue*" — true, but the same-tissue control does not address the
same-lab confound.

### 3.3 No uncertainty quantification anywhere in the project

The headline is "6/6 whole-cell datasets rank `Oprl1` first". One of those six is:

```
NodoMap:Buchanan   Oprl1 9.5932   vs   Oprm1 9.4257     margin: 1.8%
```

A 1.8 % gap between pseudobulk means over 4,599 cells is counted as a win identically to Bai's
27 % gap. There is no bootstrap, no confidence interval, no permutation test, and no
per-dataset dispersion in any of the 16 output tables. A cell-level bootstrap of the rank
statistic is roughly ten lines and would tell the reader which of the six results are secure.

### 3.4 The `Penk` divergence is computed on prep-mixed pooled data

`README.md:105` cites `log2` enrichment −1.47 for `Penk` in nodose/jugular neurons. That comes
from `nodose_neuron_enrichment.csv`, built at `02_nodose.py:137-138` from `is_neuron` across
**all five** NodoMap datasets — including the nuclear `inhouse` set.

The project's central methodological argument is that whole-cell and nuclear preparations must
not be pooled before a tissue claim is made (`04_synthesis.py:8-9`: "datasets are stratified by
preparation before any tissue claim is made"). This table pools them. `Penk` is short (5 kb) so
the nuclear contribution is unlikely to flip the sign, but the number as published is
internally inconsistent with the stated method and should be recomputed whole-cell-only.

### 3.5 Two README labelling errors

- **`README.md:75-78`** labels the `NodoMap:Bai/Buchanan/Kupari/Zhao` rows as tissue
  **`nodose`**. Those rows are computed over `is_neuron`, which is *nodose plus jugular*
  (`02_nodose.py:79-80`), and `ligand_receptor_by_tissue.csv` correctly labels them
  `nodose+jugular`. The README's own next row separates "jugular" as a distinct tissue, so the
  table as printed implies a nodose-only split that was not computed. Same issue at
  `README.md:91` ("1.44 CPM in nodose neurons" — that is `NodoMap:inhouse`, nodose+jugular).
- **`README.md:93`**, "Glu15 (n = 98, 40.0 CPM, **7.1 % of nuclei**)". 7.14 is `Pnoc_pct` — the
  fraction of *Glu15 nuclei with any detected `Pnoc`*. Glu15 is 98/49,392 = **0.20 %** of nuclei.
  As written the sentence reads as a population size and is off by 35×.

---

## Severity 4 — robustness and correctness risks not currently triggered

### 4.1 `receptor_rank()` crashes on a `NaN` level

`atlas_common.py:79` — `.astype(int)` on a rank containing `NaN`:

```
IntCastingNaNError: Cannot convert non-finite values (NA or NaN) to integer
```

### 4.2 `receptor_rank()` on an all-zero sample produces four rank-1 rows

Verified: `method="min"` gives every gene rank 1 and `fraction_of_top` = `NaN`. That propagates
to `top_receptor_by_dataset` as four rows for one dataset, and `04_synthesis.py:127`
(`tt.loc[ds, "top_receptor"]`) then returns a Series into `colour.get(...)` → `TypeError:
unhashable`. A guard for `top <= 0` returning an explicit "undetermined" would be clearer.

### 4.3 The whole-cell baseline is an unweighted mean over datasets of very unequal size

`04_synthesis.py:82-84`:

| dataset | cells | weight given |
|---|---|---|
| NodoMap:Bai | 1,031 | 25 % |
| NodoMap:Buchanan | 4,599 | 25 % |
| NodoMap:Kupari | 1,993 | 25 % |
| **NodoMap:Zhao** | **23,017** | **25 %** |

Zhao contributes 74 % of the cells and 25 % of the weight. Defensible as one-vote-per-study,
but it is undocumented and it sets the denominator of every ratio in the gene-length analysis.
Worth a one-line comment stating the choice, and a check that the conclusion survives
cell-weighting.

### 4.4 `nuclear_over_whole_cell` has no zero guard

`04_synthesis.py:91` divides without a pseudocount. `Oprd1` already sits at 0.22 CPM; a gene at
0 would yield `inf`, which flows into `np.log10` at line 96 and then into `np.corrcoef`,
silently returning `NaN` for *r*.

### 4.5 Duplicate gene symbols in the NodoMap atlas would be silently undercounted

`02_nodose.py:40`: `int(hit if np.isscalar(hit) else hit.iloc[0])` takes the **first** matching
row and discards the rest. `atlas_common.collapse_isoforms()` exists for exactly this problem
and is not used here. `01_geniculate.py:39` does sum duplicates. The two pipelines disagree.

### 4.6 The Zuker orientation heuristic can silently transpose the matrix

`01_geniculate.py:49`: `if x.shape[0] < x.shape[1]: x = x.T`. Correct for 454 cells × ~20k
genes, but it is an unvalidated guess — for any matrix with more cells than genes it silently
produces a transposed result and CPM normalisation along the wrong axis. Checking for known
gene symbols in the index is unambiguous and costs one line.

Related: the code assumes `GSM4037432_..._expressed454.xlsx` holds raw counts and divides by
column sums. Nothing verifies that. If the file is already per-cell normalised the rank and
ratio claims survive (both are within-sample and scale-invariant), but the absolute CPM values
in `geniculate_opioid_levels.csv` would be wrong.

### 4.7 `02_nodose.py` documents a fallback it does not implement

`02_nodose.py:62-64`:

```python
# nCount_RNA is the library size the authors recorded; fall back to the sum
# over the genes present only if it is missing.
total = obs["nCount_RNA"].values.astype(float)
```

There is no fallback. A missing column raises `KeyError`.

### 4.8 `03_nts.py` guards three genes and not the fourth

`03_nts.py:143-152` wraps `Pnoc`, `Slc17a6` and `Gad1` in `if … in cpm` guards; `Oprl1`
(lines 147-148) is unguarded and would raise `KeyError`. Same table, same loop.

### 4.9 `00_download_data.sh` cannot detect a failed download

`curl -sS -O` without `-f` writes the HTTP error body to the output file and exits 0, so
`set -euo pipefail` does not catch it. A GEO 404 or a proxy error page lands on disk named
`GSE166648_snRNA_unnormdata.csv.gz`, and the failure surfaces much later as a confusing parse
error. `-L` is also absent, so a redirect is saved as the payload.

Fix: `curl -fsSL -O`. For a reproducibility-first project, checksums on the four downloaded
files would be a further improvement.

---

## Severity 5 — maintainability

- **Dead code in `atlas_common.py`.** `abundance_matched_percentile` (30 lines, `:122-151`),
  `collapse_isoforms` (`:158`), `stream_gene_rows` (`:165`) and `TISSUE_COLORS` (`:47`) have no
  callers. `stream_gene_rows` is a near-duplicate of the streaming loop inlined at
  `03_nts.py:49-65` — two implementations of the same non-trivial parser, only one of which is
  exercised. `abundance_matched_percentile` is referenced obliquely by `02_nodose.py:135-136`
  ("the same abundance-matched statistic"), which then computes a plain log2 ratio instead.
- **`04_synthesis.py:144`**: `v2 = v[v.gene != "Pnoc"]` is a no-op — `Pnoc` has no genomic span
  and was already removed by `dropna` at line 94. Harmless, but it implies a filter that is not
  doing anything.
- **`04_synthesis.py:177`**: panel C is titled "Receptor exceeds its own ligand in every
  **peripheral** dataset" while the plotted selection (lines 162-165) includes
  `NTS / GSE166648`, which is central. The panel also shows 8 bars against the README table's 6
  rows; `NodoMap:inhouse` and NTS appear in the figure and not in the table.
- **`04_synthesis.py:92`**: `bias.reset_index().rename(columns={"index": "gene"})` works only
  because the source Series carry *conflicting* index names (`"Gene"` from the span groupby,
  `"gene"` elsewhere), which makes pandas set the combined index name to `None`. Verified: if
  the names are ever made consistent, `reset_index()` emits a `Gene` column, the rename becomes
  a no-op, and `v.gene` raises `AttributeError`. `.rename_axis("gene")` before `reset_index` is
  explicit.
- **`04_synthesis.py:76-78`** re-reads the GSE102443 FPKM file, re-hardcoding the filename that
  `01_geniculate.py:23` already defines as `FPKM_FILE`. It also means the synthesis step cannot
  run from `results/` alone — and since `data/raw/` is git-ignored, Figure 1 cannot be
  regenerated from a fresh clone without re-downloading ~110 MB.
- **`02_nodose.py:91`** stamps `ASSAY_LABEL["NodoMap"]` = "…(106,436 cells)" onto every row,
  including the 765-cell `inhouse` row whose own `n_cells` column reads 765.
- **`atlas_common.py:36`**: `NODOSE_ROOT` is a hardcoded absolute path
  (`/home/user/PNOC-Nodose`). The README tells the reader to edit the source. An
  `os.environ.get("NODOSE_ROOT", …)` default would keep the same behaviour without a source
  edit.
- **`02_nodose.py:79`**: `f"nodose+jugular"` — f-string with no placeholder.
- **`02_nodose.py:114-119`**: the rank-printing loop relies on `b.tissue.iloc[0]` happening to
  be `"nodose"` after a `str.contains("nodose")` filter that also matches `"nodose+jugular"`.
  Correct today, entirely dependent on row order.
- **Repo hygiene.** No `LICENSE`, no tests, no CI, and `requirements.txt` pins nothing — for a
  project whose thesis is that recomputation changes conclusions, an unpinned pandas/numpy is a
  real reproducibility gap (pandas 3.0 changed several `groupby`/`reset_index` behaviours this
  code depends on).
- **Figure 1 panel B**: the `Penk` and `Pomc` labels overlap in the rendered PNG despite the
  manual offsets at `04_synthesis.py:149`.

---

## Recommended order of work

1. Fix `03_nts.py:39` (§1.1) — one character class, restores the cache.
2. Derive `in_matrix` honestly in `02_nodose.py` (§1.2) — restores the project's own guard
   against its central error.
3. Re-run all four scripts and commit regenerated tables (§2.1) — establishes that the repo
   reproduces its own results.
4. Add leave-one-out *r* to the gene-length section and soften panel B's title (§3.1); name the
   prep/lab confound and the `--include-introns` alternative (§3.2).
5. Bootstrap the rank statistic (§3.3) — the largest single gain in credibility for the effort.
6. Fix the three README labelling errors (§3.5) and either implement or drop the sanity gate
   (§2.2).
7. `curl -fsSL` (§4.9), zero guards (§4.4), delete dead code (§5).
