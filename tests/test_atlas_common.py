"""Unit tests for the statistics and loaders every published number passes through.

Run with:  python3 -m pytest tests -q
"""

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import atlas_common as ac  # noqa: E402


# --------------------------------------------------------------- receptor_rank

def test_receptor_rank_orders_and_marks_determinate():
    out = ac.receptor_rank(pd.Series({"Oprl1": 10.0, "Oprm1": 5.0,
                                      "Oprd1": 1.0, "Oprk1": 0.5}))
    assert list(out.gene) == ["Oprl1", "Oprm1", "Oprd1", "Oprk1"]
    assert list(out["rank"]) == [1, 2, 3, 4]
    assert out.determinate.all()
    assert out.fraction_of_top.iloc[1] == pytest.approx(0.5)


def test_receptor_rank_rejects_nan_rather_than_crashing_on_cast():
    # A NaN level means the caller has a data bug: an unmeasured gene must be
    # absent from the index, not present as NaN.
    with pytest.raises(ValueError, match="non-finite"):
        ac.receptor_rank(pd.Series({"Oprl1": 1.0, "Oprm1": np.nan,
                                    "Oprd1": 0.5, "Oprk1": 0.2}))


def test_receptor_rank_all_zero_is_not_a_winner():
    out = ac.receptor_rank(pd.Series({g: 0.0 for g in ac.RECEPTORS}))
    assert not out.determinate.any()
    assert out.fraction_of_top.isna().all()


def test_receptor_rank_exact_tie_is_not_determinate():
    out = ac.receptor_rank(pd.Series({"Oprl1": 5.0, "Oprm1": 5.0,
                                      "Oprd1": 1.0, "Oprk1": 0.5}))
    assert not out.determinate.any()


def test_receptor_rank_ignores_genes_absent_from_the_index():
    out = ac.receptor_rank(pd.Series({"Oprl1": 3.0, "Oprm1": 1.0}))
    assert set(out.gene) == {"Oprl1", "Oprm1"}


# ------------------------------------------------- bootstrap_receptor_support

def test_bootstrap_support_is_high_for_a_clear_winner():
    rng = np.random.default_rng(0)
    per_cell = pd.DataFrame({
        "Oprl1": rng.normal(100, 5, 500), "Oprm1": rng.normal(10, 5, 500),
        "Oprd1": rng.normal(1, 1, 500), "Oprk1": rng.normal(1, 1, 500)})
    out = ac.bootstrap_receptor_support(per_cell, n_boot=200, seed=1)
    assert out["top_gene"] == "Oprl1"
    assert out["support"] == 1.0


def test_bootstrap_support_is_not_confident_when_the_top_two_are_level():
    # Two receptors drawn from the same distribution: one of them wins the
    # sample by noise, and support has to stay well short of certainty. This is
    # the NodoMap:Buchanan case, where Oprl1 leads Oprm1 by 1.8%.
    rng = np.random.default_rng(0)
    per_cell = pd.DataFrame({
        "Oprl1": rng.normal(10, 8, 800), "Oprm1": rng.normal(10, 8, 800),
        "Oprd1": rng.normal(1, 1, 800), "Oprk1": rng.normal(1, 1, 800)})
    out = ac.bootstrap_receptor_support(per_cell, n_boot=400, seed=2)
    assert 0.2 < out["support"] < 0.95
    assert out["margin_lo"] < 1.0 < out["margin_hi"]  # CI spans "no difference"


def test_bootstrap_support_is_reproducible_for_a_fixed_seed():
    per_cell = pd.DataFrame(np.random.default_rng(3).normal(5, 2, (200, 4)),
                            columns=ac.RECEPTORS)
    a = ac.bootstrap_receptor_support(per_cell, n_boot=100, seed=7)
    b = ac.bootstrap_receptor_support(per_cell, n_boot=100, seed=7)
    assert a == b


def test_bootstrap_support_handles_an_empty_signal():
    per_cell = pd.DataFrame(np.zeros((50, 4)), columns=ac.RECEPTORS)
    assert np.isnan(ac.bootstrap_receptor_support(per_cell, n_boot=10)["support"])


# --------------------------------------------------- transcriptome_percentile

def test_transcriptome_percentile_places_a_gene_in_its_own_distribution():
    levels = pd.Series({f"g{i}": float(i) for i in range(1, 201)})
    levels["target"] = 100.5
    pct = ac.transcriptome_percentile(levels, "target")
    assert pct == pytest.approx(100 * 100 / 201, abs=0.5)


def test_transcriptome_percentile_needs_a_real_distribution():
    assert np.isnan(ac.transcriptome_percentile(
        pd.Series({f"g{i}": float(i) for i in range(50)}), "g10"))


def test_transcriptome_percentile_of_an_absent_gene_is_nan():
    assert np.isnan(ac.transcriptome_percentile(pd.Series({"a": 1.0}), "Oprl1"))


# ---------------------------------------------------------------- check_markers

def test_check_markers_passes_a_healthy_matrix():
    levels = pd.Series({g: 100.0 for g in ac.SANITY_GENES})
    tbl = ac.check_markers(levels, "ok")
    assert tbl.in_matrix.all()
    assert set(tbl[tbl.enforced].gene) == set(ac.REQUIRED_MARKERS)


def test_check_markers_raises_when_a_required_marker_is_missing():
    levels = pd.Series({g: 100.0 for g in ac.SANITY_GENES if g != "Snap25"})
    with pytest.raises(ac.SanityCheckError, match="Snap25"):
        ac.check_markers(levels, "no-neurons")


def test_check_markers_raises_when_a_required_marker_is_silent():
    levels = pd.Series({g: 100.0 for g in ac.SANITY_GENES})
    levels["Actb"] = 0.0
    with pytest.raises(ac.SanityCheckError, match="Actb"):
        ac.check_markers(levels, "dead-matrix")


def test_check_markers_tolerates_a_missing_tissue_specific_marker():
    levels = pd.Series({g: 100.0 for g in ac.SANITY_GENES if g != "Phox2b"})
    tbl = ac.check_markers(levels, "no-phox2b")
    assert not tbl.set_index("gene").loc["Phox2b", "in_matrix"]


# ------------------------------------------------------------ stream_gene_rows

def _write(tmp_path, text, gz=False):
    p = tmp_path / ("m.csv.gz" if gz else "m.csv")
    if gz:
        with gzip.open(p, "wt") as fh:
            fh.write(text)
    else:
        p.write_text(text)
    return p


def test_stream_gene_rows_reads_requested_rows(tmp_path):
    p = _write(tmp_path, '"","c1","c2"\n"Oprl1",1,2\n"Oprm1",3,4\n')
    frame, absent, totals = ac.stream_gene_rows(p, ["Oprl1"])
    assert list(frame.index) == ["c1", "c2"]
    assert frame["Oprl1"].tolist() == [1.0, 2.0]
    assert absent == [] and totals is None


def test_stream_gene_rows_reports_absent_genes_rather_than_zeroing_them(tmp_path):
    p = _write(tmp_path, "gene,c1,c2\nOprl1,1,2\n")
    frame, absent, _ = ac.stream_gene_rows(p, ["Oprl1", "Pnoc"], totals=True)
    assert absent == ["Pnoc"]
    assert "Pnoc" not in frame.columns


def test_stream_gene_rows_totals_cover_every_row(tmp_path):
    p = _write(tmp_path, "gene,c1,c2\nOprl1,1,2\nOprm1,3,4\nActb,5,6\n")
    _, _, totals = ac.stream_gene_rows(p, ["Oprl1"], totals=True)
    assert totals.to_dict() == {"Oprl1": 3.0, "Oprm1": 7.0, "Actb": 11.0}


def test_stream_gene_rows_handles_gzip(tmp_path):
    p = _write(tmp_path, "gene,c1,c2\nOprl1,1,2\n", gz=True)
    frame, _, _ = ac.stream_gene_rows(p, ["Oprl1"])
    assert frame["Oprl1"].tolist() == [1.0, 2.0]


def test_stream_gene_rows_names_a_bad_header_offset(tmp_path):
    # Every data row disagrees with the header width: the old code produced an
    # opaque empty-array error much later.
    p = _write(tmp_path, "gene,c1,c2,c3\nOprl1,1,2\nOprm1,3,4\n")
    with pytest.raises(ValueError, match="header offset"):
        ac.stream_gene_rows(p, ["Oprl1"], totals=True)


def test_stream_gene_rows_raises_when_nothing_matches(tmp_path):
    p = _write(tmp_path, "gene,c1\nActb,1\n")
    with pytest.raises(ValueError, match="none of"):
        ac.stream_gene_rows(p, ["Oprl1"])


# ----------------------------------------------------------- csr_gene_columns

class _FakeCSR(dict):
    """Minimal stand-in for an h5py CSR group."""

    def __init__(self, dense):
        import scipy.sparse as sp
        m = sp.csr_matrix(np.asarray(dense, dtype=float))
        super().__init__(data=m.data, indices=m.indices, indptr=m.indptr)
        self.attrs = {"shape": np.array(m.shape)}


def test_csr_gene_columns_extracts_the_requested_columns():
    dense = [[1, 0, 3], [0, 5, 0], [7, 0, 9]]
    out = ac.csr_gene_columns(_FakeCSR(dense), {"a": [0], "c": [2]}, 3)
    assert out.tolist() == [[1, 3], [0, 0], [7, 9]]


def test_csr_gene_columns_sums_duplicate_annotation_rows():
    # A symbol on two annotation rows must be summed, not silently truncated to
    # the first row, which would undercount the gene.
    dense = [[1, 2, 0], [0, 4, 0], [5, 6, 0]]
    out = ac.csr_gene_columns(_FakeCSR(dense), {"dup": [0, 1]}, 3)
    assert out.ravel().tolist() == [3, 4, 11]


# ------------------------------------------------------- leave_one_out_pearson

def test_leave_one_out_pearson_exposes_a_two_point_correlation():
    # Five points with no trend plus two extreme ones: the full r is high and
    # collapses when either extreme point is dropped.
    x = [2, 5, 6, 6, 18, 34, 250]
    y = [0.99, 0.39, 0.47, 0.36, 0.95, 24.4, 29.0]
    out = ac.leave_one_out_pearson(x, y, ["Pdyn", "Penk", "Oprl1", "Pomc",
                                          "Oprk1", "Oprd1", "Oprm1"])
    full = out[out.dropped == "(none)"].pearson_r.iloc[0]
    assert full > 0.8
    assert out[out.dropped == "Oprm1"].pearson_r.iloc[0] < full
    assert len(out) == 8


def test_leave_one_out_pearson_reports_the_full_fit_first():
    out = ac.leave_one_out_pearson([1, 2, 3], [1, 2, 3], ["a", "b", "c"])
    assert out.dropped.iloc[0] == "(none)"
    assert out.delta_vs_full.iloc[0] == 0.0


# ------------------------------------------------------------------- registry

def test_dataset_registry_is_internally_consistent():
    assert set(ac.PREP) == set(ac.TISSUE_OF) == set(ac.DATASETS)
    assert set(ac.PREP.values()) == {"whole cell", "nuclear"}
    assert ac.WHOLE_CELL_NODOSE == ["NodoMap:Bai", "NodoMap:Buchanan",
                                    "NodoMap:Kupari", "NodoMap:Zhao"]


def test_no_ligand_receptor_ratio_survives_in_the_api():
    # This project reports Oprl1 expression. Receptor-to-ligand ratios were
    # removed deliberately; re-adding one should fail a test, not pass review.
    assert not hasattr(ac, "ligand_receptor_ratio")
