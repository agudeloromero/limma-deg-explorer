# ============================================================
# app.py  —  Limma DEG Explorer (Python Shiny)
#
# Author:      Patricia Agudelo-Romero
# Email:       patricia.agudeloromero@thekids.org.au
# Institution: The Kids Research Institute Australia
#
# Rules illustrated in this file:
#   Rule 3  — data logic (analysis.py, report.py) is separate from display
#   Rule 4  — uploads are validated at the boundary (read_limma_csv)
#   Rule 5  — every threshold is a UI control
#   Rule 6  — the analysis state is downloadable (report.py)
#   Rule 10 — the UI describes decisions, not implementation details
#
# Differences from the R implementation (Rule 2):
#   plotting        matplotlib + adjustText   (R: ggplot2 + ggrepel)
#   GO enrichment   gseapy / Enrichr API      (R: clusterProfiler, offline)
#                   run on request in a background task, because it is a
#                   network call (R recomputes automatically)
#   report          Jinja2 + base64 figures   (R: R Markdown)
#
# Run locally:   shiny run app.py --reload
# ============================================================

from __future__ import annotations

import asyncio
from datetime import datetime
from pathlib import Path

import pandas as pd
from shiny import App, reactive, render, req, ui

from analysis import (
    GO_LABELS,
    GO_SETS,
    MIN_DEGS_FOR_ENRICHMENT,
    REQUIRED_COLUMNS,
    deg_table,
    filter_degs,
    go_cache_key,
    go_dotplot_figure,
    go_table,
    is_go_cached,
    message_figure,
    rank_barplot_figure,
    read_limma_csv,
    run_go,
    volcano_figure,
)
from report import build_report_zip

APP_DIR = Path(__file__).resolve().parent
EXAMPLE_PATH = APP_DIR / "data" / "GSE60450_limma.csv"


def threshold_input(id_: str, label: str, min_: float, max_: float, value: float, step: float):
    """A slider paired with a numeric box (kept in sync by the server)."""
    return ui.layout_columns(
        ui.input_slider(id_, label, min=min_, max=max_, value=value, step=step),
        ui.div(
            {"class": "num-box"},
            ui.input_numeric(f"{id_}_num", None, min=min_, max=max_, value=value, step=step),
        ),
        col_widths=[8, 4],
    )


# ── UI ───────────────────────────────────────────────────────

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_file("file", "Upload limma results (CSV)", accept=[".csv"]),
        (
            ui.input_action_link("load_example", "or load the example dataset (GSE60450)")
            if EXAMPLE_PATH.exists()
            else None
        ),
        ui.hr(),
        ui.h5("DEG filters", class_="section dge"),
        threshold_input("deg_pval", "Adjusted P-value ≤", 0.001, 1, 0.05, 0.001),
        threshold_input("deg_lfc", "|log₂FC| ≥", 0, 4, 1.0, 0.1),
        ui.hr(),
        ui.h5("Enrichment filters", class_="section enr"),
        ui.input_select("geneset", "Gene set shown on screen", {k: GO_LABELS[k] for k in GO_SETS}),
        threshold_input("enr_pval", "Enrichment adj. P-value ≤", 0.001, 1, 0.05, 0.001),
        ui.p(
            "All three GO ontologies are included in the download. "
            "Enrichment uses the Enrichr web service (internet required).",
            class_="hint",
        ),
        ui.hr(),
        ui.h5("Report", class_="section"),
        ui.output_ui("download_ui"),
        width=320,
    ),
    ui.tags.style("""
        body { font-family: Arial, sans-serif; }
        .num-box { margin-top: 26px; }
        h5.section { font-weight: bold; color: #2C3E50; }
        h5.dge { color: #2E74B5; }  h5.enr { color: #1E8449; }
        .hint { font-size: 12px; color: #888; }
        .summary { background: #EBF3FB; padding: 10px; border-radius: 4px;
                   margin-bottom: 14px; font-size: 14px; }
        .problem { color: #B00020; font-weight: bold; padding: 8px 0; }
        .go-bar { display: flex; gap: 16px; align-items: center; margin: 12px 0; }
    """),
    ui.output_ui("validation_msg"),
    ui.output_ui("deg_summary"),
    ui.output_ui("main_panel"),
    title="Limma DEG Explorer",
    fillable=False,
)


# ── Server ───────────────────────────────────────────────────


def server(input, output, session):

    # Slider <-> numeric box synchronisation (Rule 5)
    def sync(slider_id: str):
        num_id = f"{slider_id}_num"

        @reactive.effect
        @reactive.event(input[slider_id])
        def _to_num():
            if input[num_id]() != input[slider_id]():
                ui.update_numeric(num_id, value=input[slider_id]())

        @reactive.effect
        @reactive.event(input[num_id])
        def _to_slider():
            v = input[num_id]()
            if v is not None and v != input[slider_id]():
                ui.update_slider(slider_id, value=v)

    for sid in ("deg_pval", "deg_lfc", "enr_pval"):
        sync(sid)

    # ── Data source: uploaded file or example (latest action wins)
    source = reactive.value(None)  # (path, display name)

    @reactive.effect
    @reactive.event(input.file)
    def _use_upload():
        f = input.file()
        if f:
            source.set((f[0]["datapath"], f[0]["name"]))

    if EXAMPLE_PATH.exists():

        @reactive.effect
        @reactive.event(input.load_example)
        def _use_example():
            source.set((str(EXAMPLE_PATH), EXAMPLE_PATH.name))

    # ── Rule 3: reactive pipeline  loaded -> raw_data -> degs -> enrichment
    @reactive.calc
    def loaded():
        s = source()
        return (None, []) if s is None else read_limma_csv(s[0])

    @reactive.calc
    def raw_data():
        return loaded()[0]

    @reactive.calc
    def degs():
        df = raw_data()
        req(df is not None)
        return filter_degs(df, input.deg_pval(), input.deg_lfc())

    @reactive.calc
    def go_query():
        """Genes, background and library for the on-screen enrichment."""
        return (degs()["gene"].tolist(), raw_data()["gene"].tolist(), GO_SETS[input.geneset()])

    # Enrichr is a network call, so it runs in a background task and the
    # app stays responsive for everyone while it waits.
    @reactive.extended_task
    async def go_task(genes, background, library):
        return await asyncio.to_thread(run_go, genes, background, library, 1.0)

    go_task_key = reactive.value(None)  # query behind the latest task

    @reactive.effect
    @reactive.event(input.run_go)
    def _start_go():
        go_task_key.set(go_cache_key(*go_query()))
        go_task.invoke(*go_query())

    @reactive.calc
    def enrichment():
        """(GoResult or None, message) for the current filters."""
        genes, background, library = go_query()
        if len(genes) < MIN_DEGS_FOR_ENRICHMENT:
            return None, f"At least {MIN_DEGS_FOR_ENRICHMENT} DEGs are needed for enrichment."
        if is_go_cached(genes, background, library):  # instant re-filtering
            res = run_go(genes, background, library, input.enr_pval())
            if res.table is None:
                return None, f"No enriched terms at adj. P-value ≤ {input.enr_pval()}."
            return res, None
        status = go_task.status()
        if status == "running":
            return None, "Querying Enrichr … this can take 10–30 seconds."
        if (
            status == "success"
            and go_task.result().error
            and go_task_key() == go_cache_key(genes, background, library)
        ):
            return None, go_task.result().error
        return None, "Click “Run GO enrichment” to analyse the current DEGs."

    # ── Messages and summary ─────────────────────────────────
    @render.ui
    def validation_msg():
        problems = loaded()[1]
        if problems:
            return ui.div(*[ui.p(p) for p in problems], class_="problem")
        return None

    @render.ui
    def deg_summary():
        df = raw_data()
        req(df is not None)
        d = degs()
        n_up, n_dn = int((d["logFC"] > 0).sum()), int((d["logFC"] < 0).sum())
        return ui.div(
            ui.tags.b(f"{len(d)} DEGs"),
            f" from {len(df)} genes  |  ",
            ui.span(f"▲ {n_up} up", style="color:#7B2D8B;"),
            "  ",
            ui.span(f"▼ {n_dn} down", style="color:#2E8B57;"),
            ui.span(
                f"adj.P.Val ≤ {input.deg_pval()}  |  |log₂FC| ≥ {input.deg_lfc()}",
                style="color:#888; font-size:12px; margin-left:16px;",
            ),
            class_="summary",
        )

    @render.ui
    def main_panel():
        if raw_data() is None:
            if loaded()[1]:
                return None
            return ui.div(
                {"style": "margin-top:60px; text-align:center; color:#666;"},
                ui.h4("Limma DEG Explorer", style="color:#2E74B5;"),
                ui.p(
                    "Upload a limma results table (CSV) using the panel on the left"
                    + (", or load the example dataset." if EXAMPLE_PATH.exists() else ".")
                ),
                ui.p(
                    {"style": "font-size:13px;"},
                    "Required columns: ",
                    *[x for c in REQUIRED_COLUMNS for x in (ui.code(c), " ")],
                    "— the default output of ",
                    ui.code("limma::topTable()"),
                    ".",
                ),
            )
        return ui.navset_tab(
            ui.nav_panel("DEG table", ui.output_data_frame("deg_table_out")),
            ui.nav_panel("Volcano plot", ui.output_plot("volcano", height="520px")),
            ui.nav_panel("Top ranked genes", ui.output_plot("barplot", height="620px")),
            ui.nav_panel(
                "GO enrichment",
                ui.div(
                    ui.input_task_button(
                        "run_go", "Run GO enrichment", label_busy="Querying Enrichr…"
                    ),
                    ui.span(
                        "Queries the Enrichr web service for the selected ontology "
                        "and the current DEGs. The enrichment P-value slider then "
                        "filters the result instantly.",
                        class_="hint",
                    ),
                    class_="go-bar",
                ),
                ui.output_plot("enrich_dot", height="520px"),
                ui.hr(),
                ui.output_data_frame("enrich_table"),
            ),
            id="tabs",
        )

    # ── Plots and tables ─────────────────────────────────────
    @render.data_frame
    def deg_table_out():
        return render.DataGrid(deg_table(degs()).round(6), filters=True, height="560px")

    @render.plot
    def volcano():
        df = raw_data()
        req(df is not None)
        return volcano_figure(df, input.deg_pval(), input.deg_lfc())

    @render.plot
    def barplot():
        d = degs()
        if d.empty:
            return message_figure("No DEGs at the current cutoffs.")
        return rank_barplot_figure(d)

    @render.plot
    def enrich_dot():
        res, msg = enrichment()
        if res is None:
            return message_figure(msg)
        return go_dotplot_figure(
            res.table, f"{GO_LABELS[input.geneset()]} (adj. P-value ≤ {input.enr_pval()})"
        )

    @render.data_frame
    def enrich_table():
        res, msg = enrichment()
        if res is None:
            return pd.DataFrame({"Message": [msg]})
        return render.DataGrid(go_table(res.table), filters=True)

    # ── Rule 6: download the analysis folder ─────────────────
    @render.ui
    def download_ui():
        if raw_data() is None:
            return ui.p("Load data to enable the download.", class_="hint")
        return ui.TagList(
            ui.p(
                "A .zip folder with the HTML report, the DEG table (CSV) and "
                "GO enrichment for all three ontologies (Excel).",
                class_="hint",
            ),
            ui.download_button("dl_report", "Download results folder", class_="btn-primary w-100"),
        )

    download_time = {}

    def download_name():
        download_time["now"] = datetime.now()  # same timestamp inside the zip
        return f"DEG_analysis_{download_time['now']:%Y%m%d_%H%M%S}.zip"

    @render.download(filename=download_name)
    async def dl_report():
        raw, d = raw_data(), degs()
        params = {
            "deg_pval": input.deg_pval(),
            "deg_lfc": input.deg_lfc(),
            "enr_pval": input.enr_pval(),
            "input_name": source()[1],
        }
        genes, background = d["gene"].tolist(), raw["gene"].tolist()

        with ui.Progress(min=0, max=len(GO_SETS) + 1) as p:
            go_results = {}
            for i, (key, library) in enumerate(GO_SETS.items(), start=1):
                p.set(i - 1, message=f"Running {GO_LABELS[key]} …")
                # network call in a worker thread so other users are not blocked
                go_results[key] = await asyncio.to_thread(
                    run_go, genes, background, library, params["enr_pval"]
                )
            p.set(len(GO_SETS), message="Building report …")
            # matplotlib is not thread-safe, so the figures are drawn here
            payload = build_report_zip(raw, d, params, go_results, now=download_time.get("now"))
        yield payload


app = App(app_ui, server)
