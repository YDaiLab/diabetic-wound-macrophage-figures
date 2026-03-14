import pandas as pd
import numpy as np
import gseapy as gp
from pathlib import Path

# ----------------------------------------------------------------------
input_file = "box/results/seurat/rna/findallmarkers.csv"
OUTPUT_DIR = Path("results/enrichment_results")

ORGANISM = "Mouse"
ENRICHR_LIB = "MSigDB_Hallmark_2020"  # gp.get_library_name()


# ----------------------------------------------------------------------
def run_enrichr(
    gene_list: list[str],
    clust: str,
    direction: str,
) -> pd.DataFrame | None:
    # skip tiny gene sets to avoid noise
    if len(gene_list) < 10:
        return None

    outdir = OUTPUT_DIR / f"{clust}/ORA_{direction}"
    outdir.mkdir(parents=True, exist_ok=True)
    enr = gp.enrichr(
        gene_list=gene_list,
        organism=ORGANISM,
        gene_sets=ENRICHR_LIB,
        outdir=outdir,
        cutoff=0.05,
    )

    if enr.results is not None and not enr.results.empty:
        out = enr.results.copy()
        out.insert(0, "cluster", clust)
        out.insert(1, "set", f"ORA_{direction}")
        return out
    return None


def run_gsea(rnk: pd.Series, clust: str) -> pd.DataFrame | None:
    # skip small gene sets to avoid noise
    if rnk.size < 200:
        return None

    outdir = OUTPUT_DIR / f"{clust}/GSEA"
    outdir.mkdir(parents=True, exist_ok=True)
    pre = gp.prerank(
        rnk=rnk,
        gene_sets=ENRICHR_LIB,
        min_size=10,
        max_size=5000,
        permutation_num=1000,
        seed=42,
        outdir=outdir,
        processes=4,
        format="png",
    )

    if pre.res2d is not None and not pre.res2d.empty:
        g = pre.res2d.copy()
        g.insert(0, "cluster", clust)
        g.insert(1, "set", "GSEA_preranked")
        return g
    return None


# ----------------------------------------------------------------------


# Load marker DE results
df = pd.read_csv(input_file)
df["gene"] = df["gene"].astype(str)
df["score"] = np.sign(df["avg_log2FC"]) * (-np.log10(df["p_val_adj"] + 1e-300))


all_ora = []
all_gsea = []
for clust, sub in df.groupby("cluster", sort=False):
    # ORA on up/down regulated genes
    up = sub.query("p_val_adj < 0.05 & avg_log2FC > 0")["gene"].tolist()
    down = sub.query("p_val_adj < 0.05 & avg_log2FC < 0")["gene"].tolist()

    ora_up = run_enrichr(up, clust=clust, direction="up")
    ora_dn = run_enrichr(down, clust=clust, direction="down")
    for t in (ora_up, ora_dn):
        if t is not None:
            all_ora.append(t)

    # GSEA preranked on ALL genes
    rnk = (
        sub[["gene", "score"]]
        .dropna()
        .sort_values("score", ascending=False)
        .drop_duplicates("gene", keep="first")
        .set_index("gene")["score"]
    )
    g = run_gsea(rnk, clust=clust)
    if g is not None and not g.empty:
        all_gsea.append(g)

# Write combined results
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
if all_ora:
    ora_df = pd.concat(all_ora, ignore_index=True)
    ora_df.to_csv(OUTPUT_DIR / "ORA.csv", index=False)

if all_gsea:
    gsea_df = pd.concat(all_gsea, ignore_index=True)
    gsea_df.to_csv(OUTPUT_DIR / "GSEA.csv", index=False)

print("Done. Results saved in:", OUTPUT_DIR.resolve())
