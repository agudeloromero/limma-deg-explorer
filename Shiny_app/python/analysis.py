"""
analysis.py  —  Data logic for the Limma DEG Explorer (Python Shiny)

Rule 3: everything that computes lives here, outside the Shiny server,
so the same functions drive the on-screen plots, the downloadable
report and the unit tests.

Author:      Patricia Agudelo-Romero
Institution: The Kids Research Institute Australia
"""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

# ── Constants ────────────────────────────────────────────────

REQUIRED_COLUMNS = ["gene", "logFC", "adj.P.Val"]
OPTIONAL_COLUMNS = ["AveExpr", "t", "P.Value", "B"]
TABLE_COLUMNS = ["gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"]

COLOURS = {"NS": "#B3B3B3", "Up": "#7B2D8B", "Down": "#2E8B57"}

# Enrichr gene-set libraries (R uses clusterProfiler ont = "BP"/"MF"/"CC")
GO_SETS = {
    "GO_BP": "GO_Biological_Process_2023",
    "GO_MF": "GO_Molecular_Function_2023",
    "GO_CC": "GO_Cellular_Component_2023",
}
GO_LABELS = {
    "GO_BP": "GO Biological Process",
    "GO_MF": "GO Molecular Function",
    "GO_CC": "GO Cellular Component",
}
MIN_DEGS_FOR_ENRICHMENT = 5


# ── Rule 4: validate at the boundary ─────────────────────────


def validate_limma_table(df: pd.DataFrame) -> list[str]:
    """Return a list of human-readable problems (empty list = valid)."""
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        return [
            f"Missing required column(s): {', '.join(missing)}. "
            f"Required: {', '.join(REQUIRED_COLUMNS)}."
        ]
    problems = []
    if df.empty:
        problems.append("The uploaded file contains no rows.")
    for col in ("logFC", "adj.P.Val"):
        if not pd.api.types.is_numeric_dtype(df[col]):
            problems.append(f"Column '{col}' must be numeric.")
    if not problems:
        p = df["adj.P.Val"].dropna()
        if ((p < 0) | (p > 1)).any():
            problems.append("Column 'adj.P.Val' must contain values between 0 and 1.")
    return problems


def read_limma_csv(path: str) -> tuple[pd.DataFrame | None, list[str]]:
    """Read and validate an uploaded limma table: (data, problems)."""
    try:
        df = pd.read_csv(path)
    except Exception as exc:  # unreadable / not a CSV
        return None, [f"Could not read the file as CSV ({exc.__class__.__name__})."]
    problems = validate_limma_table(df)
    if problems:
        return None, problems
    df = df.dropna(subset=REQUIRED_COLUMNS).copy()
    df["gene"] = df["gene"].astype(str)
    return df, []


# ── DEG filtering ────────────────────────────────────────────


def filter_degs(df: pd.DataFrame, pval: float, lfc: float) -> pd.DataFrame:
    """Genes passing adj.P.Val <= pval and |logFC| >= lfc."""
    mask = (df["adj.P.Val"] <= pval) & (df["logFC"].abs() >= lfc)
    return df.loc[mask].copy()


def deg_table(degs: pd.DataFrame) -> pd.DataFrame:
    """DEG table with the columns present, sorted by adj.P.Val."""
    cols = [c for c in TABLE_COLUMNS if c in degs.columns]
    return degs[cols].sort_values("adj.P.Val").reset_index(drop=True)


# ── GO enrichment (Enrichr via gseapy; requires internet) ────


@dataclass(frozen=True)
class GoResult:
    table: pd.DataFrame | None  # filtered, sorted terms (None if none)
    error: str | None = None  # message for the user if the query failed


def _enrichr_query(genes: tuple, background: tuple, library: str) -> pd.DataFrame:
    """One Enrichr request (unfiltered results). Raises on network/API errors."""
    import gseapy as gp

    enr = gp.enrichr(
        gene_list=list(genes),
        gene_sets=library,
        organism="mouse",
        background=list(background) if background else None,
        outdir=None,
        no_plot=True,
    )
    return enr.results.copy()


# Results are cached per (genes, background, library), so moving the
# enrichment P-value slider re-filters instantly without a new query.
_GO_CACHE: dict = {}


def go_cache_key(genes, background, library: str) -> tuple:
    # Enrichr libraries use upper-case symbols; mouse symbols are matched
    # case-insensitively this way (e.g. Krt5 -> KRT5).
    g = tuple(sorted({str(x).upper() for x in genes}))
    b = tuple(sorted({str(x).upper() for x in background})) if background is not None else ()
    return g, b, library


def is_go_cached(genes, background, library: str) -> bool:
    return go_cache_key(genes, background, library) in _GO_CACHE


def clear_go_cache() -> None:
    _GO_CACHE.clear()


def run_go(genes, background, library: str, cutoff: float) -> GoResult:
    """
    Over-representation analysis for one GO library.

    genes       DEG symbols
    background  all tested genes (the uploaded table), so the universe
                matches the experiment rather than the whole genome
    cutoff      enrichment adj. P-value threshold (the only filter applied)
    """
    if len(genes) < MIN_DEGS_FOR_ENRICHMENT:
        return GoResult(
            None, f"At least {MIN_DEGS_FOR_ENRICHMENT} DEGs are needed for enrichment analysis."
        )
    key = go_cache_key(genes, background, library)
    if key not in _GO_CACHE:
        try:
            _GO_CACHE[key] = normalise_enrichr(_enrichr_query(*key))
        except Exception as exc:
            return GoResult(
                None,
                "Enrichr could not be reached, so GO enrichment "
                f"was not run ({exc.__class__.__name__}). "
                "Check the internet connection and try again.",
            )
    df = _GO_CACHE[key]
    df = df[df["p.adjust"] <= cutoff].sort_values("p.adjust").reset_index(drop=True)
    return GoResult(df if not df.empty else None)


def normalise_enrichr(raw: pd.DataFrame) -> pd.DataFrame:
    """Rename Enrichr columns to the clusterProfiler names used in R."""
    df = raw.rename(
        columns={
            "Term": "Description",
            "P-value": "pvalue",
            "Adjusted P-value": "p.adjust",
            "Genes": "geneID",
            "Combined Score": "CombinedScore",
        }
    ).copy()
    if "Overlap" in df.columns:
        df["Count"] = df["Overlap"].astype(str).str.split("/").str[0].astype(int)
        df["GeneRatio"] = df["Overlap"].astype(str)
    if "p.adjust" not in df.columns:
        df["p.adjust"] = pd.Series(dtype=float)
    return df


GO_TABLE_COLUMNS = ["Description", "GeneRatio", "pvalue", "p.adjust", "Count", "geneID"]


def go_table(df: pd.DataFrame) -> pd.DataFrame:
    return df[[c for c in GO_TABLE_COLUMNS if c in df.columns]]


# ── Figures (shared by the app and the report) ───────────────


def _neglog10(p: pd.Series) -> pd.Series:
    return -np.log10(p.clip(lower=1e-300))


def volcano_figure(df: pd.DataFrame, pval: float, lfc: float, n_labels: int = 5):
    sig = (df["adj.P.Val"] <= pval) & (df["logFC"].abs() >= lfc)
    direction = np.where(~sig, "NS", np.where(df["logFC"] > 0, "Up", "Down"))

    fig, ax = plt.subplots(figsize=(8, 6))
    for grp in ("NS", "Up", "Down"):
        sub = df[direction == grp]
        ax.scatter(
            sub["logFC"],
            _neglog10(sub["adj.P.Val"]),
            c=COLOURS[grp],
            s=8 if grp == "NS" else 14,
            alpha=0.5 if grp == "NS" else 0.8,
            linewidths=0,
            label=f"{grp} ({len(sub)})",
        )

    ax.axvline(lfc, ls="--", color="black", lw=0.6)
    ax.axvline(-lfc, ls="--", color="black", lw=0.6)
    ax.axhline(-np.log10(pval), ls="--", color="black", lw=0.6)

    top = pd.concat(
        [
            df[sig & (df["logFC"] > 0)].nsmallest(n_labels, "adj.P.Val"),
            df[sig & (df["logFC"] < 0)].nsmallest(n_labels, "adj.P.Val"),
        ]
    )
    texts = [
        ax.text(r["logFC"], -np.log10(max(r["adj.P.Val"], 1e-300)), r["gene"], fontsize=8)
        for _, r in top.iterrows()
    ]
    if texts:
        try:
            from adjustText import adjust_text

            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))
        except ImportError:
            pass

    ax.set_xlabel(r"$\log_2$ fold change", fontsize=11)
    ax.set_ylabel(r"$-\log_{10}$(adj. P-value)", fontsize=11)
    ax.set_title("Volcano plot", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def rank_barplot_figure(degs: pd.DataFrame, n: int = 30):
    top = degs.assign(rank=np.sign(degs["logFC"]) * _neglog10(degs["adj.P.Val"]))
    top = top.reindex(top["rank"].abs().sort_values(ascending=False).index).head(n)
    top = top.sort_values("rank")

    fig, ax = plt.subplots(figsize=(8, max(4, 0.28 * len(top) + 1.5)))
    colours = [COLOURS["Up"] if v > 0 else COLOURS["Down"] for v in top["rank"]]
    ax.barh(range(len(top)), top["rank"].values, color=colours)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["gene"].values, fontsize=8)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_xlabel(r"Rank score: sign($\log_2$FC) $\times$ $-\log_{10}$(adj. P-value)", fontsize=9)
    ax.set_title(f"Top {len(top)} DEGs by rank", fontsize=12, fontweight="bold")
    ax.legend(
        handles=[
            Patch(color=COLOURS["Up"], label="Up"),
            Patch(color=COLOURS["Down"], label="Down"),
        ],
        frameon=False,
        loc="lower right",
        fontsize=9,
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def go_dotplot_figure(df: pd.DataFrame, title: str, n: int = 15):
    top = df.nsmallest(n, "p.adjust").iloc[::-1]  # most significant at the top
    cmin, cmax = int(top["Count"].min()), int(top["Count"].max())
    span = max(cmax - cmin, 1)

    def size(c):
        return 40 + 260 * (c - cmin) / span

    fig, ax = plt.subplots(figsize=(10, max(3.5, 0.38 * len(top) + 1.5)))
    sc = ax.scatter(
        top["Count"],
        range(len(top)),
        c=top["p.adjust"],
        cmap="RdBu",
        s=size(top["Count"]),
        edgecolors="none",
    )
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([_wrap(t) for t in top["Description"]], fontsize=8)
    ax.set_xlabel("Gene count", fontsize=10)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.margins(x=0.15, y=0.06)
    cbar = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("adj. P-value", fontsize=9)
    legend_counts = sorted({cmin, (cmin + cmax) // 2, cmax})
    ax.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor="grey",
                markersize=size(v) ** 0.5,
                label=str(v),
            )
            for v in legend_counts
        ],
        title="Count",
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.18, 1),
        fontsize=8,
        title_fontsize=9,
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def message_figure(text: str):
    fig, ax = plt.subplots(figsize=(8, 2.5))
    ax.text(
        0.5,
        0.5,
        text,
        ha="center",
        va="center",
        fontsize=11,
        color="#666",
        wrap=True,
        transform=ax.transAxes,
    )
    ax.axis("off")
    return fig


def _wrap(text: str, width: int = 55) -> str:
    text = str(text)
    return text if len(text) <= width else text[: width - 1] + "…"
