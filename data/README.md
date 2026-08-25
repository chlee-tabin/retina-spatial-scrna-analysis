# Data Directory

The figure scripts read five processed objects from this directory, downloaded
from [GEO GSE322831](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE322831)
with the `GSE322831_` prefix stripped (download → rename table and copy-paste
snippet: [top-level README](../README.md#reproducing-figures)):

- `20250604_chick_RPC.h5ad`, `20250604_human_RPC.h5ad`,
  `20260528_mouse_RPC_cr9_e13e16.h5ad` (Python figure scripts)
- `20250604_01_retina.rds`, `20250604_02_fabp7.rds` (R figure scripts)

Also tracked here: `chick_W_genes.tsv` / `chick_Z_genes.tsv` (W/Z chromosome
gene lists used by the chick preprocessing provenance scripts).

Everything else that appears in this directory (analyzer `*.pkl` caches, MEX
re-export subdirectories) is generated output and is gitignored.
