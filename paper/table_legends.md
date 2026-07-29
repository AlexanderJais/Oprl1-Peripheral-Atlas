# Table legends

## Table S1. The nineteen peripheral neuronal populations of Figure 1, related to Figure 1

One row per panel of Figure 1: the deposit each population was read from, the preparation and
chemistry it was measured on, the neurons and biological samples contributing, the mean level of
each of the four opioid receptor genes, and the statistics behind the ordering. Panel letters are
those of Figure 1. The table the journal receives is
[`paper/tables/TableS1.csv`](tables/TableS1.csv), one sheet of twenty-six columns;
[`paper/tables/TableS1.md`](tables/TableS1.md) renders the same rows for reading, split into (a)
provenance and sample and (b) levels, ordering and quality control. Both are written by
`paper/tables.py` from the tables under `results/`.

Study and DOI name the publication that generated each deposit, taken from the deposit itself
rather than from a search on the tissue: eleven of the twelve GEO records name their publication,
either in the citation field or in the paper's own data availability statement. GSE309608 is the
exception and is matched by its BioProject accession (PRJNA1335542), its six contributors and its
title, all of which the PNAS paper carries. Full citations are in
[`paper/references.md`](references.md).

Divisions are the blocks of Figure 1. The pelvic ganglion carries both sympathetic and
parasympathetic neurons and is grouped with the sympathetic block, as in the figure;
`results/dataset_quality_panel.csv` records it as mixed autonomic. All nineteen populations are
whole-cell, for the reason given in Figure 2; both iPain atlases are majority single-nucleus and
are restricted here to `suspension_type == "cell"`.

Levels are pseudobulk means over all neurons of the population, in the unit each dataset was
measured in: FPKM for the full-length geniculate deposit, CPM elsewhere. Levels are compared only
within a panel, never across units or platforms. Detection is the percentage of that population's
neurons in which *Oprl1* is called, and is reported per dataset for the same reason: it scales with
sequencing depth, which differs by two orders of magnitude across these deposits. Median library is
per neuron, in UMI for droplet deposits and reads for full-length ones.

Margins are the ratio of the highest receptor to the second-highest on the same cells. Intervals
are 95% bootstrap intervals over 10,000 resamples of the cells, and support is the fraction of
those resamples that retain the observed top receptor. Samples agreeing counts the biological
samples that place the top receptor first on their own cells.

Neurons were called on raw counts with glial and immune barcodes excluded, and every population
passed the same marker gate before any receptor value was read from it. The ambient-RNA check is
the neuron-to-non-neuron ratio of each receptor within one dissociation; where it could not run,
the table gives the reason rather than a blank, and those reasons are listed under the rendered
table. A deposit of sorted or author-filtered neurons has no non-neuronal compartment to score
against, which is what most of them are.

Blank cells are quantities the source tables do not record, and are left empty rather than
carried over from a comparable population:

- **Samples** is recorded only for the populations `src/10_peripheral_ganglia.py` derived from
  source. For the rest, sample agreement is in the samples-agreeing column, which is why panel F
  reads 4/4 with no sample count beside it. `n/a` marks a deposit that does not resolve into
  biological samples here.
- **Median library** is not recorded for panels B, D, E, H and O. The two iPain populations were
  not re-derived from their source matrices, and the GSE231924 deposit carries no library sizes.
- **Detection** is not recorded for panels E and H, the two iPain populations, for the same reason.
- **Panel A** has no margin: the spiral ganglion quantifies all four receptors and detects one, so
  there is no second receptor to divide by. The three zeros are measured absences, not a reference
  gap.
