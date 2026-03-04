# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Figure S15: Grid Size Sensitivity Analysis for 2D Map Reconstruction

# %% [markdown]
# ## Load chick RPC data

# %% tags=["cell-9"]
import sys
sys.path.append('..')
from spatial_expression_analysis import SpatialAnalysisParams, SpatialExpressionAnalyzer
import anndata as ad
# chick_adata = ad.read_h5ad("../data/20240815_fabp7.h5ad")
chick_adata = ad.read_h5ad("../data/20250604_chick_RPC.h5ad")
print(f"Loaded chick data: {chick_adata.n_obs:,} cells × {chick_adata.n_vars:,} genes")
# Derive parameters transparently based on the data characteristics
chick_params = SpatialAnalysisParams.derive_parameters_from_data(
    n_cells=chick_adata.n_obs, 
    species='chick'
)
chick_params

# %% [markdown]
# ## SF15: Grid sensitivity bin=100

# %% [markdown]
# ### Grid sensitivity - bin=100

# %% tags=["cell-12"]
chick_params.bin_size = 100

# %% tags=["cell-13"]
# Initialize analyzer
chick_analyzer = SpatialExpressionAnalyzer(chick_params)
# Run full analysis (this replaces all your manual preprocessing steps)
# Update path to point to data file location relative to notebooks/ directory
# chick_results = chick_analyzer.run_full_analysis("../data/20240815_fabp7.h5ad")
chick_results = chick_analyzer.run_full_analysis("../data/20250604_chick_RPC.h5ad")

# %% tags=["cell-14"]
import matplotlib.pyplot as plt
import os
import numpy as np
# Create output directory if it doesn't exist
FIGURES_BASE = os.path.join(os.path.dirname(__file__), "..", "..", "figures")
output_dir = os.path.join(FIGURES_BASE, "Figure_SF15", "100x100")
os.makedirs(output_dir, exist_ok=True)
key_genes = [
    "FGF8", "CYP26C1",
]
for gene in key_genes:
  if gene in chick_analyzer.gene_names:
      print(f"\nProcessing {gene}...")
      # Create figure and axis
      fig, ax = plt.subplots(figsize=(6, 6))
      # Get gene index
      gene2idx = {g: i for i, g in enumerate(chick_analyzer.gene_names)}
      idx = gene2idx[gene]
      # Get the image data
      img = chick_analyzer.images[idx].copy()
      # Apply mask if counts are available
      if chick_analyzer.counts is not None:
          mask = (chick_analyzer.counts < chick_analyzer.params.mask_count_threshold)
          img = np.ma.masked_where(mask, img)
      # Plot the image
      im = ax.imshow(
          img,
          origin="lower",
          cmap="viridis",
          aspect="equal"
      )
      # Add reference lines if available
      if chick_analyzer.dv_mid is not None and chick_analyzer.nt_mid is not None:
          ax.axhline(y=chick_analyzer.dv_mid, color="black", lw=0.5)
          ax.axvline(x=chick_analyzer.nt_mid, color="black", lw=0.5)
      # Get spatial max value
      spatial_max = float(np.max(chick_analyzer.images[idx]))
      # Set title
      ax.set_title(f"{gene} (max: {spatial_max:.2f})", fontsize=12, fontweight='bold')
      # Set axis labels
      ax.set_xlabel("Temporal ← NT → Nasal", fontsize=10)
      ax.set_ylabel("Ventral ← DV → Dorsal", fontsize=10)
      # Remove axis ticks for cleaner look
      ax.set_xticks([])
      ax.set_yticks([])
      # Add colorbar
      plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
      # Adjust layout
      plt.tight_layout()
      # Save the figure
      output_path = os.path.join(output_dir, f"{gene}_spatial_expression.png")
      plt.savefig(output_path, dpi=300, bbox_inches='tight')
      plt.close()
      print(f"  Saved to {output_path}")
      # Get and save top correlations as text
      correlations = chick_analyzer.get_gene_correlations(gene, top_n=10)
      # Save correlations to text file
      corr_output_path = os.path.join(output_dir, f"{gene}_correlations.txt")
      with open(corr_output_path, 'w') as f:
          f.write(f"Top correlations for {gene}:\n")
          for corr_gene, corr_val in correlations:
              f.write(f"  {corr_gene}: {corr_val:.3f}\n")
      print(f"  Saved correlations to {corr_output_path}")
print(f"\nAll plots saved to {output_dir}/")

# %% [markdown]
# ## SF15: Grid sensitivity bin=20

# %% [markdown]
# ### Grid sensitivity - bin=20

# %% tags=["cell-16"]
chick_params.bin_size = 20

# %% tags=["cell-17"]
# Initialize analyzer
chick_analyzer = SpatialExpressionAnalyzer(chick_params)
# Run full analysis (this replaces all your manual preprocessing steps)
# Update path to point to data file location relative to notebooks/ directory
# chick_results = chick_analyzer.run_full_analysis("../data/20240815_fabp7.h5ad")
chick_results = chick_analyzer.run_full_analysis("../data/20250604_chick_RPC.h5ad")

# %% tags=["cell-18"]
import matplotlib.pyplot as plt
import os
import numpy as np
# Create output directory if it doesn't exist
output_dir = os.path.join(FIGURES_BASE, "Figure_SF15", "20x20")
os.makedirs(output_dir, exist_ok=True)
key_genes = [
    "FGF8", "CYP26C1",
]
for gene in key_genes:
  if gene in chick_analyzer.gene_names:
      print(f"\nProcessing {gene}...")
      # Create figure and axis
      fig, ax = plt.subplots(figsize=(6, 6))
      # Get gene index
      gene2idx = {g: i for i, g in enumerate(chick_analyzer.gene_names)}
      idx = gene2idx[gene]
      # Get the image data
      img = chick_analyzer.images[idx].copy()
      # Apply mask if counts are available
      if chick_analyzer.counts is not None:
          mask = (chick_analyzer.counts < chick_analyzer.params.mask_count_threshold)
          img = np.ma.masked_where(mask, img)
      # Plot the image
      im = ax.imshow(
          img,
          origin="lower",
          cmap="viridis",
          aspect="equal"
      )
      # Add reference lines if available
      if chick_analyzer.dv_mid is not None and chick_analyzer.nt_mid is not None:
          ax.axhline(y=chick_analyzer.dv_mid, color="black", lw=0.5)
          ax.axvline(x=chick_analyzer.nt_mid, color="black", lw=0.5)
      # Get spatial max value
      spatial_max = float(np.max(chick_analyzer.images[idx]))
      # Set title
      ax.set_title(f"{gene} (max: {spatial_max:.2f})", fontsize=12, fontweight='bold')
      # Set axis labels
      ax.set_xlabel("Temporal ← NT → Nasal", fontsize=10)
      ax.set_ylabel("Ventral ← DV → Dorsal", fontsize=10)
      # Remove axis ticks for cleaner look
      ax.set_xticks([])
      ax.set_yticks([])
      # Add colorbar
      plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
      # Adjust layout
      plt.tight_layout()
      # Save the figure
      output_path = os.path.join(output_dir, f"{gene}_spatial_expression.png")
      plt.savefig(output_path, dpi=300, bbox_inches='tight')
      plt.close()
      print(f"  Saved to {output_path}")
      # Get and save top correlations as text
      correlations = chick_analyzer.get_gene_correlations(gene, top_n=10)
      # Save correlations to text file
      corr_output_path = os.path.join(output_dir, f"{gene}_correlations.txt")
      with open(corr_output_path, 'w') as f:
          f.write(f"Top correlations for {gene}:\n")
          for corr_gene, corr_val in correlations:
              f.write(f"  {corr_gene}: {corr_val:.3f}\n")
      print(f"  Saved correlations to {corr_output_path}")
print(f"\nAll plots saved to {output_dir}/")

# %% [markdown]
# ## SF15: Grid sensitivity bin=51 (final)

# %% [markdown]
# ### Grid sensitivity - bin=51

# %% tags=["cell-21"]
# You could potentially change the parameters like this:
# chick_params.bin_size = 26
chick_params.bin_size = 51
chick_params.min_gene_count = 20 # to align with human
# Do not use this normally but it is the default
# chick_params_default = SpatialAnalysisParams.get_species_defaults('chick', show_reasoning=True)

# %% [markdown]
# * This is filtering genes that have too low expression.
# * It calculates the similarity matrix.

# %% tags=["cell-23"]
# Initialize analyzer
chick_analyzer = SpatialExpressionAnalyzer(chick_params)
# Run full analysis (this replaces all your manual preprocessing steps)
# Update path to point to data file location relative to notebooks/ directory
# chick_results = chick_analyzer.run_full_analysis("../data/20240815_fabp7.h5ad")
chick_results = chick_analyzer.run_full_analysis("../data/20250604_chick_RPC.h5ad")
