from pathlib import Path
import pandas as pd

# ----------------------------------------------------------------------
input_dir = "box/results/seurat/rna"
output_file = "results/de_tests.xlsx"
# ----------------------------------------------------------------------


def reorder_columns(df, first_cols=["gene", "cluster"]):
    cols = [*first_cols, *[c for c in df.columns if c not in first_cols]]
    return df[cols]


# ----------------------------------------------------------------------

# Load and reorder
df_ovr = pd.read_csv(Path(input_dir) / "findallmarkers.csv")
df_ovr = reorder_columns(df_ovr)

df_dbvnd = pd.read_csv(Path(input_dir) / "find_de.csv")
df_dbvnd = reorder_columns(df_dbvnd)

df = pd.merge(
    df_ovr,
    df_dbvnd,
    on=["gene", "cluster"],
    suffixes=("_ovr", "_db_v_nd"),
    how="outer",
)
df = reorder_columns(df)

# Save to single Excel file with multiple sheets
Path(output_file).parent.mkdir(parents=True, exist_ok=True)
with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    df_ovr.to_excel(writer, sheet_name="one_vs_rest", index=False)
    df_dbvnd.to_excel(writer, sheet_name="db_vs_nd", index=False)
    df.to_excel(writer, sheet_name="combined", index=False)
