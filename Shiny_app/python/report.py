"""
report.py  —  Rule 6: capture the analysis state in a downloadable folder

Builds DEG_analysis_<timestamp>.zip containing
    DEG_report_<timestamp>.html     self-contained HTML report (Jinja2)
    DEG_results_<timestamp>.csv     filtered DEG table
    GO_enrichment_<timestamp>.xlsx  GO_BP / GO_MF / GO_CC sheets

Author:      Patricia Agudelo-Romero
Institution: The Kids Research Institute Australia
"""

from __future__ import annotations

import base64
import io
import platform
import zipfile
from datetime import datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import matplotlib.pyplot as plt
import openpyxl
import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from analysis import (
    GO_LABELS,
    GO_SETS,
    GoResult,
    deg_table,
    go_dotplot_figure,
    go_table,
    rank_barplot_figure,
    volcano_figure,
)

APP_VERSION = "1.0.0"
TEMPLATE_DIR = Path(__file__).resolve().parent / "templates"


def _pkg_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not installed"


def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _go_workbook(go_results: dict[str, GoResult], enr_pval: float) -> bytes:
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    header_fill = PatternFill("solid", fgColor="2E74B5")
    header_font = Font(color="FFFFFF", bold=True, name="Arial", size=11)

    for key, res in go_results.items():
        if res.table is None:
            continue
        df = go_table(res.table)
        ws = wb.create_sheet(key)
        ws.append(list(df.columns))
        for cell in ws[1]:
            cell.fill, cell.font = header_fill, header_font
            cell.alignment = Alignment(horizontal="left")
        for row in df.itertuples(index=False):
            ws.append(list(row))
        for i, col in enumerate(df.columns, start=1):
            width = max(len(col), int(df[col].astype(str).str.len().max() or 0))
            ws.column_dimensions[get_column_letter(i)].width = min(width + 2, 60)
        ws.freeze_panes = "A2"

    if not wb.sheetnames:
        ws = wb.create_sheet("No_results")
        ws.append(["Message"])
        ws.append([f"No enriched GO terms at adj. P-value <= {enr_pval}"])
        for key, res in go_results.items():
            if res.error:
                ws.append([f"{key}: {res.error}"])

    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()


def build_report_zip(
    raw: pd.DataFrame,
    degs: pd.DataFrame,
    params: dict,
    go_results: dict[str, GoResult],
    now: datetime | None = None,
) -> bytes:
    """
    raw         validated uploaded table
    degs        genes passing the DEG filters
    params      dict with deg_pval, deg_lfc, enr_pval, input_name
    go_results  {"GO_BP": GoResult, "GO_MF": ..., "GO_CC": ...}
    """
    now = now or datetime.now()
    stamp = now.strftime("%Y%m%d_%H%M%S")

    # 1. DEG table
    csv_bytes = deg_table(degs).to_csv(index=False).encode()

    # 2. GO workbook
    xlsx_bytes = _go_workbook(go_results, params["enr_pval"])

    # 3. HTML report
    go_sections = []
    for key in GO_SETS:
        res = go_results.get(key, GoResult(None))
        title = f"{GO_LABELS[key]} (adj. P-value ≤ {params['enr_pval']})"
        go_sections.append(
            {
                "label": GO_LABELS[key],
                "n_terms": 0 if res.table is None else len(res.table),
                "img": None
                if res.table is None
                else _fig_to_b64(go_dotplot_figure(res.table, title)),
                "error": res.error,
            }
        )

    env = Environment(loader=FileSystemLoader(TEMPLATE_DIR), autoescape=select_autoescape(["html"]))
    html = env.get_template("report_template.html").render(
        date=now.strftime("%d %B %Y, %H:%M"),
        input_name=params.get("input_name", "uploaded file"),
        deg_pval=params["deg_pval"],
        deg_lfc=params["deg_lfc"],
        enr_pval=params["enr_pval"],
        n_total=len(raw),
        n_degs=len(degs),
        n_up=int((degs["logFC"] > 0).sum()),
        n_dn=int((degs["logFC"] < 0).sum()),
        volcano_b64=_fig_to_b64(volcano_figure(raw, params["deg_pval"], params["deg_lfc"])),
        barplot_b64=(_fig_to_b64(rank_barplot_figure(degs)) if len(degs) else None),
        go_sections=go_sections,
        go_libraries=", ".join(GO_SETS.values()),
        versions={
            "app": APP_VERSION,
            "python": platform.python_version(),
            "shiny": _pkg_version("shiny"),
            "gseapy": _pkg_version("gseapy"),
            "pandas": _pkg_version("pandas"),
            "matplotlib": _pkg_version("matplotlib"),
        },
    )

    # 4. Bundle
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(f"DEG_report_{stamp}.html", html.encode())
        zf.writestr(f"DEG_results_{stamp}.csv", csv_bytes)
        zf.writestr(f"GO_enrichment_{stamp}.xlsx", xlsx_bytes)
    return buf.getvalue()
