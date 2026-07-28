# Code audit — Oprl1-Peripheral-Atlas

Full read of `src/*.py`, `src/00_download_data.sh`, `README.md`, and all tables in `results/`,
originally at commit `cbc374f`. Every numeric claim in the README was recomputed from the
committed CSVs; every finding was reproduced rather than inferred.

**Status: all findings resolved.** The whole pipeline was then re-run end to end against freshly
downloaded source data — GSE102443, GSE135801, GSE166648 from GEO and the NodoMap integrated
atlas from CZ CELLxGENE — and every table and figure in the repository was regenerated from it.
A 28-test suite (`tests/`) now covers the shared statistics.

Two things changed the shape of the repo after the audit and are recorded here for context:

- The receptor-to-ligand ratio analysis (`Oprl1`:`Pnoc`) was removed at the author's direction.
  The project reports `Oprl1` expression; `Pnoc` is carried as a descriptive column only. This
  deleted findings 3.5 and part of 2.1 rather than fixing them, and a unit test now fails if a
  `ligand_receptor_ratio` is reintroduced.
- The gene-length / nuclear-preparation analysis was **kept as a full result** at the author's
  direction. The audit's statistical objection to it (3.1) is now stated in the README and
  plotted in the figure rather than removed.

---

## Verified correct, before and after

The README's arithmetic was sound and remains so. All seven `nuclear_over_whole_cell` ratios,
the log-log Pearson *r* (0.8356 → "0.84"), the Spearman rho and *p* (0.577 / 0.175), the 96/454
cell counts, the 91.67 % and 9.1 % detection rates, and the 49,392 neuronal nuclei all reproduce
exactly from a clean re-run.

Two claims the audit flagged as unverifiable were checked against source and are **correct**:

- `Penk` at the **97.6th** percentile of expressed genes in GSE102443. The value was right; it
  simply was not computed by any script. `transcriptome_percentile` now runs for every opioid
  gene and the number is in `geniculate_opioid_levels.csv`.
- The `Cheng et al. 2026, Cell Press Blue 1:100072` citation, which matches the sibling
  PNOC-Nodose repository and the CELLxGENE collection metadata. Left as the author has it, with
  the DOI and collection ID added.

One number was wrong: GSE102443 quantifies **17,110** genes, not 17,111. The README is corrected
and the count is now persisted in `geniculate_opioid_levels.csv`.

The nodose per-cluster results independently reproduce the sibling project's published figure
(`figure2_oprl1_ranked_bars`): NGN14 27.7 %, NGN6 22.2 %, NGN19 22.1 %, NGN16 18.8 %, NGN7
18.6 %. Two pipelines, written separately, agreeing to the decimal.

---

## Severity 1 — broken → fixed

### 1.1 The GSE166648 cache path always raised `TypeError` — FIXED

`03_nts.py:39` evaluated `set(cached) | z["absent"].tolist()`. `set | list` is a `TypeError`, so
the cache was written on the first run and crashed on every run after it, before the fallback
re-stream could be reached.

Fixed by converting to a set. `allow_pickle=True` was also removed from the `np.load` — every
stored array is plain strings or floats, so it was never needed, and it was an
arbitrary-code-execution surface on a cache file. **Verified**: the second run of `03_nts.py`
now prints `[cache] gse166648_targets.npz` and completes in seconds instead of re-streaming.

### 1.2 A gene missing from the NodoMap annotation was recorded as a measured zero — FIXED

`02_nodose.py` hardcoded `in_matrix: True` while `load_counts` silently skipped absent genes,
so a reference gap would have been published as `mean_level = 0.0, in_matrix = True`. This is
the exact error the project's "Corrections" section exists to document.

`load_counts` now returns the set of symbols it found, and both `in_matrix` and the level are
derived from it — an unmeasured gene gets a null level, never `0.0`. The same treatment was
applied to `01_geniculate.py` and `03_nts.py` so all three pipelines agree.

---

## Severity 2 — results not reproducible from the code → fixed

### 2.1 `Oprl1_over_Pnoc` was `NaN` in the CSVs but `inf` in the code — RESOLVED BY REMOVAL

The committed tables carried an empty field for GSE102443, which `ligand_receptor_ratio` could
not produce: it returned `np.inf`, and pandas serialises that as the string `inf`. The
committed results had been produced by a different version of that function.

The receptor-to-ligand analysis has since been removed entirely, so the discrepancy is moot.
The underlying distinction it turned on — a gene absent from the annotation is not a measured
zero — is preserved as the `in_matrix` convention above, which is the part that mattered.

### 2.2 The README's sanity-gene gate was not implemented — FIXED

The README claimed each dataset was "checked against `Snap25`, `Phox2b`, `Slc17a6`, `Tac1`,
`Calca` and `Actb` before any opioid number is read from it". In fact `01_geniculate.py` printed
the levels without gating on them, and `02_nodose.py` and `03_nts.py` loaded the markers and
never looked at them.

`atlas_common.check_markers()` now implements the gate for real. `Snap25` and `Actb` are
enforced — absent or zero raises `SanityCheckError` — while the tissue-specific markers are
recorded without being enforced. All three pipelines call it before reading any opioid number
and write the result to `*_marker_checks.csv`.

### 2.3 Three README numbers were produced by no script — FIXED

`Penk`'s percentile and the gene count are now computed and persisted (see above); the citation
was verified rather than changed.

---

## Severity 3 — overstated claims → stated accurately

### 3.1 The gene-length mechanism rests on two of seven points — NOW STATED, ANALYSIS KEPT

| points | *n* | log-log Pearson *r* |
|---|---|---|
| all | 7 | **0.836** |
| drop `Oprm1` | 6 | 0.695 |
| **drop `Oprm1` + `Oprd1`** | **5** | **0.016** |

Remove the two genes that gain signal and the relationship is absent, not merely weakened.

The author elected to keep this as a full result. It is kept, and the fragility is now part of
it rather than something a reader has to notice: `leave_one_out_pearson()` computes the table,
`nuclear_bias_sensitivity.csv` records it, panel C of the synthesis figure plots it, the README
states the *r* = 0.02 figure in the same paragraph as the *r* = 0.84 figure, and the figure
panel title was changed from "The inversion **is** a gene-length effect" to "tracks gene
length". The `Oprd1`-against-a-0.22-CPM-floor problem and the `Pomc`-is-ambient problem are
named in the README.

### 3.2 Preparation type is perfectly confounded with laboratory — NOW STATED

One nuclear nodose dataset (765 neurons, in-house), no in-house whole-cell dataset: "nuclear"
and "in-house" cannot be separated in this design. The alternative explanation — `Oprm1`
detection rising from 5–9 % to 71.9 %, which fits an intron-inclusive alignment as well as
pre-mRNA retention — is now named in both the README and `atlas_common.DATASETS`.

Preparation type is also no longer hardcoded: `02_nodose.py` reads it from the atlas's own
`suspension_type` field and raises if it disagrees with the registry.

### 3.3 No uncertainty quantification anywhere — FIXED, AND IT MATTERED

`bootstrap_receptor_support()` resamples cells 2,000 times and reports how often the observed
top receptor stays top, with a 95 % interval on the margin. Running it changed what the
headline claim can say:

| dataset | margin | bootstrap support |
|---|---|---|
| NodoMap:Zhao | 1.19× | 1.00 |
| NodoMap:Bai | 1.27× | 0.97 |
| NodoMap:Kupari | 1.33× | 0.94 |
| **NodoMap:Buchanan** | **1.02×** | **0.54** |

Buchanan's "`Oprl1` first" is a coin flip. The README now says so, `04_synthesis.py` prints a
warning for any dataset below 0.95 support, and the synthesis figure prints those in red. The
6/6 count is retained but is no longer presented as six independent wins.

### 3.4 `Penk` enrichment was computed on prep-mixed pooled data — FIXED

`nodose_neuron_enrichment.csv` pooled all five NodoMap datasets, including the nuclear one,
contradicting the stratification the rest of the project rests on. It is now computed from the
whole-cell datasets only. The value barely moved (log2 −1.47 → **−1.46**), which is the right
outcome: the number was defensible, the method was not.

### 3.5 Two README labelling errors — ONE FIXED, ONE REMOVED

The "7.1 % of nuclei" error (a detection rate presented as a population size, off by 35×) was in
the ligand section and went with it. The `nodose` / `nodose+jugular` mislabelling is fixed:
per-dataset rows are labelled for what they are throughout.

---

## Severity 4 — robustness → fixed

| # | finding | fix |
|---|---|---|
| 4.1 | `receptor_rank()` crashed on a `NaN` level (`IntCastingNaNError`) | raises a clear `ValueError` naming the genes; a missing gene must be absent from the index, not `NaN` |
| 4.2 | an all-zero sample produced four rank-1 rows, which reached `colour.get(Series)` as an unhashable `TypeError` | `receptor_rank` returns `determinate = False` for an all-zero or tied sample; `04_synthesis.py` raises if two datasets tie for top |
| 4.3 | the whole-cell baseline was an unweighted mean over datasets of very unequal size (Zhao: 74 % of cells, 25 % of the weight) | both are computed; `nuclear_bias_vs_gene_length.csv` carries a `_cellweighted` column. The conclusion survives (`Oprm1` 29.1× → 25.4×) |
| 4.4 | `nuclear_over_whole_cell` divided with no zero guard, so a zero level would send `inf` into `np.log10` and silently `NaN` the correlation | the denominator is masked to strictly positive values |
| 4.5 | duplicate gene symbols were silently truncated to the first annotation row | `csr_gene_columns()` sums duplicate rows and reports which genes were affected; unit-tested |
| 4.6 | the Zuker orientation heuristic (`shape[0] < shape[1]`) would silently transpose any matrix with more cells than genes | orientation is decided by which axis carries recognisable gene symbols, and raises if neither does. A non-integer matrix now warns that CPM would be meaningless |
| 4.7 | `02_nodose.py` documented an `nCount_RNA` fallback it did not implement | `library_size()` implements it |
| 4.8 | `03_nts.py` guarded three genes in the subtype loop and left `Oprl1` unguarded | `Oprl1`'s absence raises `SanityCheckError` explicitly — there would be no result to report |
| 4.9 | `curl -sS -O` wrote HTTP error bodies to disk and exited 0, invisible to `set -e` | `curl -fsSL` with retries, atomic `.partial` renames, and SHA-256 verification of all five source files |

Also fixed: `stream_gene_rows` now counts malformed rows and raises a named "header offset"
error if every row disagrees with the header, instead of surfacing later as an empty-array
error; and `03_nts.py` raises if the GSE166648 metadata has duplicate barcodes, which would
have silently duplicated nuclei through `meta.loc[common]`.

---

## Severity 5 — maintainability → fixed

- **Dead code removed.** `abundance_matched_percentile` (30 lines), `collapse_isoforms` and the
  stale `TISSUE_COLORS` are gone. `stream_gene_rows` was a near-duplicate of the streaming loop
  inlined in `03_nts.py`; there is now one implementation, in `atlas_common`, used by `03_nts.py`
  and covered by six unit tests.
- **`v2 = v[v.gene != "Pnoc"]`** was a no-op (`Pnoc` has no span and was already dropped by
  `dropna`). Removed.
- **Panel titles now match their contents.** The old panel C was titled "in every *peripheral*
  dataset" while plotting the central NTS dataset, and showed 8 bars against the README's 6
  rows. That panel is gone; the current figures were checked against their own captions.
- **`bias.reset_index().rename(columns={"index": "gene"})`** worked only because the source
  Series carried *conflicting* index names, which made pandas set the combined name to `None`.
  Replaced with an explicit `.rename_axis("gene")`.
- **The duplicated FPKM filename literal** is now `atlas_common.GENICULATE_FPKM`, used by both
  `01_geniculate.py` and `04_synthesis.py`.
- **`ASSAY_LABEL`** no longer stamps "(106,436 cells)" onto a row whose own `n_cells` reads 765;
  cell counts live in the `n_cells` column only.
- **`NODOSE_ROOT`** reads from the environment, defaulting to the previous hardcoded path, so
  relocating the sibling project no longer requires a source edit.
- **`f"nodose+jugular"`** (f-string with no placeholder) and the order-dependent rank-printing
  loop are gone.
- **Repo hygiene.** `requirements.txt` now carries lower bounds and records the exact versions
  the pipeline was last verified against. `tests/` has 28 unit tests. `00_download_data.sh`
  fetches every input including the NodoMap atlas and checksums all five.

## Two tiers of dataset, removed

Datasets entered the project at different times and were held to different standards. The
geniculate and vagal data passed `check_markers()`, and NodoMap alone had an ambient-RNA check;
every ganglion added afterwards lived in a notes file with an ad-hoc loader and neither check. The
split was historical, not methodological.

- **`src/external/` is gone.** `src/10_peripheral_ganglia.py` now handles every ganglion outside
  the geniculate and vagal pipelines with one QC policy (2,000 UMI, 1,000 genes), one neuron
  definition (raw-count positive markers, CPM ceilings on `Sox10`, `Plp1` and `Ptprc`), one marker
  gate and one ambient check. It reproduced every previously published number exactly.
- **`ac.ambient_enrichment()`** is new and covered by four unit tests, including the pseudocount
  that keeps a measured zero in the non-neuronal compartment from returning an infinite ratio.
- **`src/03_nts.py`** now writes the same ambient table, which closed the last population where
  the check was possible but had not been run.
- **`src/11_quality_panel.py`** states, for all 21 populations, which of the four checks ran, what
  each returned, and the reason where one could not. Seven populations are neurons-only deposits
  and one, the iPain dorsal root, has a 20.9 GB source matrix above this session's disk allowance.
  Those reasons are in the table rather than absent from it.

The check changed how one number is read. Raw *Oprl1* neuronal enrichment ranges from +0.75 to
+3.48 log2 across populations, which looks like a large difference in how neuronal the transcript
is. It is not: *Snap25* tracks it, and the difference between them spans only −0.67 to +0.58. The
populations with low *Oprl1* enrichment are the ones where the non-neuronal pool still contains
neurons, which the positive control reveals and the raw number hides.

## Still open

- **No `LICENSE` file.** Choosing one is the author's call, not the auditor's.
- **No CI.** The test suite runs in under a second and would sit naturally in a GitHub Actions
  workflow; not added without direction.
- **`Pomc` in the gene-length regression.** It contributes 11.5 CPM in nodose neurons at 11.5 %
  detection with a log2 neuronal enrichment of 0.31 — the signature of ambient RNA for a
  hypothalamic/pituitary transcript. It is flagged in the README but still one of the seven
  points. Dropping it moves *r* by −0.002, so nothing turns on it either way.
