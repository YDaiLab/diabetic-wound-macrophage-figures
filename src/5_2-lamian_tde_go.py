"""
GO enrichment analysis on Lamian TDE results.

Filters genes by fdr.overall < 0.05, groups by branch_name + peak_time_group,
and runs ORA (MSigDB Hallmark) on each group using gseapy with a local .gmt file.
"""

import pandas as pd
import gseapy as gp
from pathlib import Path

# --- Config -----------------------------------------------------------------
INPUT_FILE = "results/lamian/rna_tde.csv"
OUTPUT_FILE = "results/lamian/rna_tde_go.csv"
GMT_FILE = "results/enrichment_results/1/GSEA/gene_sets.gmt"  # local MSigDB Hallmark

FDR_CUTOFF = 0.05
MIN_GENES = 10  # minimum gene-list size to attempt enrichment

# --- Load and filter --------------------------------------------------------
df = pd.read_csv(INPUT_FILE)
df["gene"] = df["gene"].astype(str)
sig = df[df["fdr.overall"] < FDR_CUTOFF].copy()
print(f"Total genes: {len(df):,}  |  Significant (FDR < {FDR_CUTOFF}): {len(sig)}")

# --- Run enrichment per group -----------------------------------------------
all_results = []

for (branch, peak_group), sub in sig.groupby(
    ["branch_name", "peak_time_group"], sort=False
):
    gene_list = sub["gene"].dropna().unique().tolist()
    n = len(gene_list)
    label = f"{branch} | {peak_group}"

    if n < MIN_GENES:
        print(f"  SKIP  {label}  ({n} genes < {MIN_GENES})")
        continue

    print(f"  RUN   {label}  ({n} genes) ...", end=" ", flush=True)

    try:
        enr = gp.enrich(
            gene_list=gene_list,
            gene_sets=GMT_FILE,
            outdir=None,
            cutoff=1.0,  # keep all terms; filter downstream if needed
        )
        res = enr.results
        if res is not None and not res.empty:
            res = res.copy()
            res.insert(0, "branch_name", branch)
            res.insert(1, "peak_time_group", peak_group)
            res.insert(2, "n_input_genes", n)
            all_results.append(res)
            n_sig = (res["Adjusted P-value"] < 0.05).sum()
            print(f"{len(res)} terms ({n_sig} with adj-p < 0.05)")
        else:
            print("no results")
    except Exception as e:
        print(f"ERROR: {e}")

# --- Combine and save -------------------------------------------------------
if all_results:
    out = pd.concat(all_results, ignore_index=True)
    Path(OUTPUT_FILE).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT_FILE, index=False)
    print(f"\nSaved {len(out)} rows to {OUTPUT_FILE}")
    print(f"Significant terms (adj-p < 0.05): {(out['Adjusted P-value'] < 0.05).sum()}")
else:
    print("\nNo enrichment results to save.")
