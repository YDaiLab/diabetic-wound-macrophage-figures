from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import combine_pvalues, norm
from statsmodels.stats.multitest import multipletests

# ----------------------------------------------------------------------
input_file = "results/bitfam_runs/pvals.csv"
output_file = "results/bitfam_runs/pvals_combined.csv"


# ----------------------------------------------------------------------
def stouffer_signed(pvals, effects):
    pvals = np.asarray(pvals)
    effects = np.asarray(effects)

    mask = (pvals > 0) & np.isfinite(pvals) & np.isfinite(effects)
    if not np.any(mask):
        return np.nan

    pvals, effects = pvals[mask], effects[mask]

    # convert two-sided p-values to one-sided z magnitude
    z_scores = norm.isf(pvals / 2) * np.sign(effects)
    # combine signed Zs
    z_comb = np.sum(z_scores) / np.sqrt(len(z_scores))
    # convert back to two-sided p-value
    p_comb = 2 * norm.sf(abs(z_comb))
    return p_comb


# Example tests:
# stouffer_signed1([0.01, 0.01], [1, 1])
# stouffer_signed([0.01, 0.01], [1, 1])

# stouffer_signed1([0.01, 0.01], [1, -1])
# stouffer_signed([0.01, 0.01], [1, -1])

# stouffer_signed1([0.01, 0.01], [-1, 1])
# stouffer_signed([0.01, 0.01], [-1, 1])

# stouffer_signed1([0.01, 0.01], [-1, -1])
# stouffer_signed([0.01, 0.01], [-1, -1])

df = pd.read_csv(input_file)

# Apply per gene+cluster
meta_df = (
    df.groupby(["cluster", "gene"])
    .apply(
        lambda g: pd.Series(
            {
                "runs": g["run"].nunique(),
                "p_val_list": list(g["p_val"]),
                "p_val_adj_list": list(g["p_val_adj"]),
                "avg_log2FC_list": list(g["avg_log2FC"]),
                "p_val": stouffer_signed(g["p_val"], g["avg_log2FC"]),
                "avg_log2FC": g["avg_log2FC"].mean(),
            }
        )
    )
    .reset_index()
)

meta_df["p_val_adj"] = multipletests(meta_df["p_val"], method="fdr_bh")[1]

Path(output_file).parent.mkdir(parents=True, exist_ok=True)
meta_df.to_csv(output_file, index=False)
