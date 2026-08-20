# moma

Figures, supplementary data files, and the Pseudotime Graph Diffusion (PGD)
embedding for *Transcriptional regulators predicted to induce macrophage
dysregulation during impaired wound healing in diabetic mice*.

This repository turns the outputs of the analysis pipeline into the panels and
tables that appear in the manuscript. The pipeline itself — network
construction, BITFAM activity inference, Lamian trajectory analysis, and
CellOracle perturbation simulation — lives in
[macrophages](https://github.com/brandonlukas/macrophages).

## Layout

| Path | Contents |
| --- | --- |
| `make_figures/` | One directory per figure: `f1`–`f6` for main figures, `sf1`–`sf4` and `revision/` for supplementary figures. Each script writes PNG, PDF and SVG. |
| `make_tables/` | `mk_s1.R`–`mk_s8.R`, one per Supplementary Data file. |
| `src/` | Numbered post-processing steps that turn pipeline output into the tidy tables the figure scripts read. |
| `src/revision/` | Analyses added during revision: trajectory reproducibility, PAGA cross-check, cluster occupancy. |
| `pgdiffusion/` | Vendored copy of [PGDiffusion](https://github.com/brandonlukas/pgdiffusion), pinned to the exact version the published analysis ran against. See `pgdiffusion/PROVENANCE.md`. |
| `assets/` | Shared ggplot theme, embedding theme, and the cluster colour palette. |

## Running

Scripts expect to be run from the repository root, and read pipeline output
from paths configured at the top of `make_figures/cfg.R` and
`make_tables/cfg.R`. Point those at your copy of the results before running.

```bash
Rscript make_figures/run_all.R     # all main and supplementary figures
Rscript make_tables/run_all.R      # Supplementary Data 1-8
```

`run_all.R` reports per-script errors and continues, so a missing input does
not abort the whole run.

Figures are assembled into their final multi-panel layouts outside this
repository, in Inkscape.

## Note on Supplementary Data 5 and 6

`mk_s5.R` and `mk_s6.R` are numbered the opposite way round from the
Supplementary Data files as published: published Data 5 is pseudobulk
regulator expression and Data 6 is trajectory-associated regulator activity,
while the scripts produce them in the reverse order. Re-running `run_all.R`
therefore swaps the contents of those two files relative to their published
legends.

## Requirements

R with tidyverse, Seurat, ComplexHeatmap, patchwork, ggh4x, legendry, ggarrow,
ggarchery, ggrastr, shadowtext, Palo, Lamian and arrow; Python with scanpy,
anndata, torch and gseapy for the `src/` steps and PGD.

## License

MIT — see [LICENSE](LICENSE).
