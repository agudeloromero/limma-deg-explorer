"""
Tests for the Python Shiny app (run from Shiny_app/python):  pytest -q

The Enrichr web service is mocked, so the tests run offline.
"""

import io
import zipfile

import numpy as np
import openpyxl
import pandas as pd
import pytest

import analysis
from analysis import GO_SETS, filter_degs, read_limma_csv, run_go, validate_limma_table
from report import build_report_zip


def _member(z: zipfile.ZipFile, ext: str) -> bytes:
    return z.read(next(n for n in z.namelist() if n.endswith(ext)))


@pytest.fixture
def limma_table():
    rng = np.random.default_rng(1)
    n = 400
    p = np.sort(
        np.concatenate(
            [
                rng.uniform(1e-12, 1e-3, 120),  # clear DEGs
                rng.uniform(1e-3, 1, n - 120),
            ]
        )
    )
    return pd.DataFrame(
        {
            "gene": [f"Gene{i}" for i in range(n)],
            "logFC": rng.choice([-1, 1], n) * rng.uniform(1.2, 4, n),
            "AveExpr": rng.uniform(0, 10, n),
            "P.Value": p,
            "adj.P.Val": np.minimum(p * 3, 1),
            "B": rng.normal(0, 1, n),
        }
    )


@pytest.fixture
def fake_enrichr(monkeypatch):
    """Replace the Enrichr query with a deterministic offline result."""
    calls = []

    def fake(genes, background, library):
        calls.append((genes, background, library))
        g = list(genes)
        return pd.DataFrame(
            {
                "Gene_set": library,
                "Term": [f"{library} term {i}" for i in range(20)],
                "Overlap": [f"{5 + i}/100" for i in range(20)],
                "P-value": np.linspace(1e-8, 0.2, 20),
                "Adjusted P-value": np.linspace(1e-6, 0.5, 20),
                "Combined Score": np.linspace(100, 1, 20),
                "Genes": [";".join(g[:3])] * 20,
            }
        )

    analysis.clear_go_cache()
    monkeypatch.setattr(analysis, "_enrichr_query", fake)
    return calls


# ── Rule 4: validation ───────────────────────────────────────


def test_valid_table_passes(limma_table):
    assert validate_limma_table(limma_table) == []


def test_minimal_columns_are_enough(limma_table):
    assert validate_limma_table(limma_table[["gene", "logFC", "adj.P.Val"]]) == []


def test_missing_columns_are_named(limma_table):
    msg = validate_limma_table(limma_table.drop(columns=["logFC"]))
    assert len(msg) == 1 and "logFC" in msg[0]


def test_non_numeric_and_out_of_range(limma_table):
    bad = limma_table.assign(logFC="x")
    assert any("logFC" in m for m in validate_limma_table(bad))
    bad = limma_table.assign(**{"adj.P.Val": 2.0})
    assert any("between 0 and 1" in m for m in validate_limma_table(bad))


def test_unreadable_file(tmp_path):
    f = tmp_path / "bad.csv"
    f.write_bytes(b"\x00\x01\x02")
    df, problems = read_limma_csv(str(f))
    assert df is None and problems


def test_read_drops_missing_values(tmp_path, limma_table):
    t = limma_table.copy()
    t.loc[0, "logFC"] = np.nan
    f = tmp_path / "t.csv"
    t.to_csv(f, index=False)
    df, problems = read_limma_csv(str(f))
    assert problems == [] and len(df) == len(t) - 1


# ── DEG filtering ────────────────────────────────────────────


def test_filter_degs(limma_table):
    d = filter_degs(limma_table, 0.05, 1.0)
    assert (d["adj.P.Val"] <= 0.05).all() and (d["logFC"].abs() >= 1).all()
    assert len(filter_degs(limma_table, 1, 0)) == len(limma_table)


# ── GO enrichment ────────────────────────────────────────────


def test_run_go_filters_and_uses_background(limma_table, fake_enrichr):
    d = filter_degs(limma_table, 0.5, 0.5)
    res = run_go(d["gene"], limma_table["gene"], GO_SETS["GO_BP"], 0.1)
    assert res.error is None
    assert (res.table["p.adjust"] <= 0.1).all()
    assert {"Description", "Count", "GeneRatio", "geneID"} <= set(res.table.columns)
    genes, background, _ = fake_enrichr[0]
    assert len(background) == len(limma_table)  # tested genes as universe
    assert all(g.isupper() for g in genes)


def test_run_go_too_few_genes(limma_table):
    res = run_go(["A", "B"], limma_table["gene"], GO_SETS["GO_BP"], 0.05)
    assert res.table is None and "5 DEGs" in res.error


def test_run_go_reports_network_errors(limma_table, monkeypatch):
    def boom(*_):
        raise ConnectionError("offline")

    analysis.clear_go_cache()
    monkeypatch.setattr(analysis, "_enrichr_query", boom)
    res = run_go(limma_table["gene"][:50], limma_table["gene"], GO_SETS["GO_BP"], 0.05)
    assert res.table is None and "Enrichr could not be reached" in res.error


# ── Rule 6: report folder ────────────────────────────────────


def _report(limma_table, go_results):
    d = filter_degs(limma_table, 0.05, 1.0)
    params = {"deg_pval": 0.05, "deg_lfc": 1.0, "enr_pval": 0.1, "input_name": "test.csv"}
    return d, zipfile.ZipFile(io.BytesIO(build_report_zip(limma_table, d, params, go_results)))


def test_report_zip_contents(limma_table, fake_enrichr):
    d = filter_degs(limma_table, 0.05, 1.0)
    go = {k: run_go(d["gene"], limma_table["gene"], lib, 0.1) for k, lib in GO_SETS.items()}
    d, z = _report(limma_table, go)

    names = z.namelist()
    assert len(names) == 3
    html = _member(z, ".html").decode()
    assert "&le; 0.05" in html and "test.csv" in html and html.count("data:image/png") == 5

    csv = pd.read_csv(io.BytesIO(_member(z, ".csv")))
    assert len(csv) == len(d) and csv["adj.P.Val"].is_monotonic_increasing

    wb = openpyxl.load_workbook(io.BytesIO(_member(z, ".xlsx")))
    assert wb.sheetnames == ["GO_BP", "GO_MF", "GO_CC"]


def test_report_without_enrichment_still_builds(limma_table, monkeypatch):
    analysis.clear_go_cache()
    monkeypatch.setattr(
        analysis, "_enrichr_query", lambda *_: (_ for _ in ()).throw(ConnectionError())
    )
    d = filter_degs(limma_table, 0.05, 1.0)
    go = {k: run_go(d["gene"], limma_table["gene"], lib, 0.1) for k, lib in GO_SETS.items()}
    _, z = _report(limma_table, go)
    html = _member(z, ".html").decode()
    assert "Enrichr could not be reached" in html
    wb = openpyxl.load_workbook(io.BytesIO(_member(z, ".xlsx")))
    assert wb.sheetnames == ["No_results"]


def test_report_with_no_degs(limma_table):
    empty = limma_table.iloc[0:0]
    params = {"deg_pval": 0.05, "deg_lfc": 1.0, "enr_pval": 0.1}
    z = zipfile.ZipFile(io.BytesIO(build_report_zip(limma_table, empty, params, {})))
    html = _member(z, ".html").decode()
    assert "No DEGs at the current cutoffs" in html


# ── The app object loads ─────────────────────────────────────


def test_app_imports():
    import app

    assert app.app is not None
