# ============================================================
# 01_data_prep.py
# Rule 7 — Validate your app against a public, citable dataset
#
# Author:      Patricia Agudelo-Romero
# Email:       patricia.agudeloromero@thekids.org.au
# Institution: The Kids Research Institute Australia
#
# Downloads GSE60450 (mouse mammary gland luminal vs basal,
# GEO) and produces a CSV in the format the Python Shiny app
# expects.
#
# Difference from R (Rule 2):
#   R uses limma-voom for DEA — the gold standard for RNA-seq.
#   Python uses a moderated t-test with BH correction via scipy
#   and statsmodels as an approximation. For publishable results,
#   use the R implementation or rpy2 to call limma from Python.
#
# Difference from R (Step 3):
#   R uses getGEOSuppFiles() which downloads the combined count
#   matrix automatically. GEOparse downloads individual sample
#   files instead, so we fetch the combined file directly from
#   the GEO FTP server.
#
# Citation:
#   GSE60450 — Mouse mammary gland RNA-seq (luminal vs basal).
#   GEO accession GSE60450.
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450
# ============================================================

import os
import urllib.request
import GEOparse
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

os.makedirs("data", exist_ok=True)


# ── 1. Download expression data from GEO ─────────────────────
print("Downloading GSE60450 metadata from GEO...")
gse = GEOparse.get_GEO("GSE60450", destdir="data/", silent=True)


# ── 2. Extract sample metadata ───────────────────────────────
pheno = pd.DataFrame({
    gsm_name: gsm.metadata
    for gsm_name, gsm in gse.gsms.items()
}).T

# Flatten list values to strings
for col in pheno.columns:
    pheno[col] = pheno[col].apply(
        lambda x: x[0] if isinstance(x, list) and len(x) == 1 else x
    )


# ── 3. Download combined count matrix directly from GEO FTP ──
# Difference from R: getGEOSuppFiles() fetches the combined
# matrix automatically. GEOparse downloads individual sample
# files, so we fetch the combined file directly from GEO FTP.
GEO_FTP = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/"
    "GSE60450/suppl/GSE60450_Lactation-GenewiseCounts.txt.gz"
)
count_file = os.path.join("data", "GSE60450_Lactation-GenewiseCounts.txt.gz")

if os.path.exists(count_file):
    print(f"Count file already present: {count_file}")
else:
    print("Downloading combined count matrix from GEO FTP...")
    print(f"  URL: {GEO_FTP}")
    urllib.request.urlretrieve(GEO_FTP, count_file)
    print(f"  Saved to: {count_file}")

print(f"Reading counts from: {count_file}")
counts = pd.read_csv(count_file, sep="\t", index_col=0)

# Drop the gene length column (first column)
counts = counts.iloc[:, 1:]
print(f"Count matrix: {counts.shape[0]} genes × {counts.shape[1]} samples")


# ── 4. Define groups from sample metadata ────────────────────
# immunophenotype is stored as one item inside the
# characteristics_ch1 list, e.g.:
#   ['strain/background: FVB/N', ..., 'immunophenotype: luminal cell population', ...]
# We search each sample's list for the immunophenotype entry.

def extract_immunophenotype(char_list):
    for item in char_list:
        if "immunophenotype" in item.lower():
            return item.split(":", 1)[-1].strip()
    return None

pheno["immunophenotype"] = pheno["characteristics_ch1"].apply(
    extract_immunophenotype
)

if pheno["immunophenotype"].isna().all():
    raise ValueError("Could not find immunophenotype in characteristics_ch1.")

# Build mapping: GSM ID → count matrix column name
# description field contains e.g.:
#   ['Sample name: MCL1-LA', 'MCL1-LA_BC2CTUACXX_GATCAG_L001_R1', ...]
# The second item is the full column name in the count matrix.
def extract_col_name(desc_list):
    if len(desc_list) >= 2:
        return desc_list[1].strip()
    return None

pheno["col_name"] = pheno["description"].apply(extract_col_name)

# Remap pheno index from GSM IDs to count matrix column names
pheno_mapped = pheno.set_index("col_name")

groups = pheno_mapped["immunophenotype"].str.contains(
    "luminal", case=False, na=False
).map({True: "luminal", False: "basal"})

print("Sample groups assigned:")
print(groups.value_counts().to_string())

luminal_idx = groups[groups == "luminal"].index.tolist()
basal_idx   = groups[groups == "basal"].index.tolist()


# ── 5. Differential expression (scipy approximation of limma) ─
# Difference from R: limma-voom applies precision weights based
# on the mean-variance relationship of counts. Here we use a
# two-sample t-test on log2-CPM values, which is an approximation.
# For true limma results, use the R implementation.

print("Computing log2-CPM and running t-tests...")

# log2-CPM normalisation
lib_sizes = counts.sum(axis=0)
cpm       = counts.divide(lib_sizes, axis=1) * 1e6
log_cpm   = np.log2(cpm + 1)

# Filter lowly expressed genes (>= 1 CPM in at least 6 samples)
keep = (cpm >= 1).sum(axis=1) >= 6
log_cpm = log_cpm[keep]
print(f"Genes after low-expression filter: {log_cpm.shape[0]}")

# Two-sample t-test: basal vs luminal
results = []
for gene in log_cpm.index:
    x_lum  = log_cpm.loc[gene, luminal_idx].values
    x_bas  = log_cpm.loc[gene, basal_idx].values
    t, p   = stats.ttest_ind(x_bas, x_lum)           # basal - luminal
    logfc  = x_bas.mean() - x_lum.mean()
    ave    = log_cpm.loc[gene].mean()
    results.append({
        "entrez_id": gene,
        "logFC":     logfc,
        "AveExpr":   ave,
        "P.Value":   p,
    })

res = pd.DataFrame(results)
res["P.Value"] = res["P.Value"].fillna(1.0)

# BH multiple testing correction
_, adj_pvals, _, _ = multipletests(res["P.Value"], method="fdr_bh")
res["adj.P.Val"] = adj_pvals
res["B"]         = 0.0          # B-statistic not available; set to 0


# ── 6. Map Entrez IDs → gene symbols ─────────────────────────
# Difference from R: R uses org.Mm.eg.db (offline).
# Here we use mygene which queries an online API.
print("Mapping Entrez IDs to gene symbols...")
try:
    import mygene
    mg     = mygene.MyGeneInfo()
    entrez = res["entrez_id"].astype(str).tolist()
    result = mg.querymany(
        entrez, scopes="entrezgene",
        fields="symbol", species="mouse",
        as_dataframe=True, verbose=False
    )
    # Keep only rows with a symbol and drop duplicates
    result = result[["symbol"]].dropna()
    result = result[~result.index.duplicated(keep="first")]
    mapping = result["symbol"].to_dict()
    res["gene"] = res["entrez_id"].astype(str).map(mapping)
    print(f"  Mapped {res['gene'].notna().sum()} of {len(res)} genes")
except Exception as e:
    print(f"  mygene mapping failed: {e}")
    print("  Falling back to Entrez IDs as gene names.")
    res["gene"] = res["entrez_id"].astype(str)

# Remove unmapped genes
res = res.dropna(subset=["gene"])
res = res[res["gene"] != ""]
print(f"Genes with symbols: {len(res)}")


# ── 7. Validate and write CSV ─────────────────────────────────
required = ["gene", "logFC", "AveExpr", "P.Value", "adj.P.Val", "B"]
assert all(c in res.columns for c in required), \
    f"Missing columns: {set(required) - set(res.columns)}"

out_path = os.path.join("data", "GSE60450_limma.csv")
res[required].sort_values("adj.P.Val").to_csv(out_path, index=False)

top = res.sort_values("adj.P.Val").iloc[0]
print(f"""
Done. Written {len(res)} rows to {out_path}
Top hit: {top['gene']}  (adj.P.Val = {top['adj.P.Val']:.2e},  logFC = {top['logFC']:.2f})

NOTE: This CSV was produced using a t-test approximation, not
limma-voom. For the paper's validation results, use the CSV
produced by the R script 01_data_prep.R — both apps accept
the same CSV format.

Validation checklist:
  [ ] Upload GSE60450_limma.csv → summary banner appears
  [ ] DEG sliders update volcano and bar plot in real time
  [ ] GO enrichment tab shows results (requires internet)
  [ ] Download Results Folder → ZIP with HTML, CSV, Excel
""")
