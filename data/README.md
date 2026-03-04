# Data Directory

Place your spatial single-cell RNA-seq data files here.

## Expected File Format

- **Format**: AnnData (.h5ad files)
- **Required columns in `.obs`**:
  - `DV.Score`: Dorsal-ventral spatial coordinate
  - `NT.Score`: Nasal-temporal spatial coordinate
- **Gene expression**: Standard AnnData format in `.X`

## Example Files

Your data files should be named descriptively, for example:
- `20240815_fabp7.h5ad` (chick FABP7+ cells)
- `20240815_human_RPC.h5ad` (human retinal progenitor cells)
- `20240815_mouse_RPC.h5ad` (mouse retinal progenitor cells)

## Data Sources

Document your data sources here:
- **Chick**: [Source/paper citation]
- **Human**: [Source/paper citation]  
- **Mouse**: [Source/paper citation]

## Processing Notes

- Ensure spatial coordinates are properly normalized
- Verify gene symbols are standardized for your species
- Remove any unwanted cell types if focusing on specific populations

**Note**: Data files are ignored by git due to their large size. Consider using git-lfs or external storage for data sharing. 