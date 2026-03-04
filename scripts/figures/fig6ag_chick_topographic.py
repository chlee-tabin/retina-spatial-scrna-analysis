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
# # Figure 6A-G: 2D Topographic Maps from Chicken scRNA-seq

# %% [markdown]
# ## Load chick RPC data and configure parameters

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
# ## Create spatial analyzer (bin_size=51)

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

# %% [markdown]
# ## F6A-G: 2D topographic maps by gene group

# %% tags=["cell-42"]
import matplotlib.pyplot as plt
import os
import numpy as np
# Output directory setup
FIGURES_BASE = os.path.join(os.path.dirname(__file__), "..", "..", "figures")
output_dir = os.path.join(FIGURES_BASE, "Figure6", "F6A-G_spatial_maps")
os.makedirs(output_dir, exist_ok=True)
key_genes = [
    "TBX2", "TBX3", "TBX5", "EFNB1", "EFNB2", "EPHB2", "VAX1", "CHRDL1",
    "ALDH1A1", "ALDH1A3", "CYP26A1", "CYP26C1", "FGF8", "FOXD1", "FOXG1", "SOHO-1", "CYP1B1", "BMP2"
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

# %% tags=["cell-43"]
import matplotlib.pyplot as plt
import os
import numpy as np
# Create output directory if it doesn't exist
output_dir = os.path.join(FIGURES_BASE, "Figure_SF17", "spatial_bmp")
os.makedirs(output_dir, exist_ok=True)
key_genes = [
"BMP2",
"BMP4",
"BMP7",
"BMPR1A",
"BMPR1B",
"BMPR2",
"ACVR2A",
"ACVR2B",
"CHRD",
"NOG",
"BAMBI",
"TWSG1",
"TSKU",
"GREM1",
"GREM2",
"BMPER",
"CHRDL1",
"CHRDL2",
"FST",
"SMAD1",
"SMAD5",
"SMAD9",
"SMAD4",
"ID1",
"ID2",
"ID3",
"MSX1",
"MSX2",
"TBX5",
"VAX2",   
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

# %% tags=["cell-44"]
import matplotlib.pyplot as plt
import os
import numpy as np
# Create output directory if it doesn't exist
output_dir = os.path.join(FIGURES_BASE, "Figure_SF17", "spatial_fgf8_downstream")
os.makedirs(output_dir, exist_ok=True)
key_genes = [
#FGF8 Downstream Genes:
"DUSP6",
"SOX2",
"EGR1",
"SNAI1",
"FOS",
# Engrailed Genes:
"EN1",
"EN2",
# Sprouty Genes:
"SPRY1",
"SPRY2",
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

# %% tags=["cell-45"]
import matplotlib.pyplot as plt
import os
import numpy as np
# Create output directory if it doesn't exist
output_dir = os.path.join(FIGURES_BASE, "Figure_SF17", "spatial_fgf_family")
os.makedirs(output_dir, exist_ok=True)
key_genes = [
# FGF Genes (Fgf1-23):
"FGF1",
"FGF2",
"FGF3",
"FGF4",
"FGF5",
"FGF6",
"FGF7",
"FGF8",
"FGF9",
"FGF10",
"FGF11",
"FGF12",
"FGF13",
"FGF14",
"FGF15",
"FGF16",
"FGF17",
"FGF18",
"FGF19",
"FGF20",
"FGF21",
"FGF22",
"FGF23",
# FGF Receptor Genes (Fgfr1-R4, FgfrL1):
"FGFR1",
"FGFR2",
"FGFR3",
"FGFR4",
"FGFRL1",
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
# ## F6A-G: Additional gene groups

# %% tags=["cell-48"]
import matplotlib.pyplot as plt
import numpy as np
# key_genes = ['TBX5', 'TBX3', 'TBX2', 'HMX1', 'SOHO-1', 'FOXG1', 'FOXD1', 'FGF8'] 
key_genes = ['EPHA7', 'EPHB2', 'BMP2', 'CYP26A1', 'CYP26C1', 'CYP1B1', 'DUSP6', 'MYOF',] 
# Create figure with subplots
fig, axes = plt.subplots(2, 4, figsize=(20, 10))
axes = axes.flatten()  # Flatten for easy indexing
# Gene to index mapping
gene2idx = {gene: i for i, gene in enumerate(chick_analyzer.gene_names)}
# Calculate spatial max values for each gene
spatial_max_values = {}
for gene in key_genes:
  if gene in chick_analyzer.gene_names:
      idx = gene2idx[gene]
      spatial_max_values[gene] = float(np.max(chick_analyzer.images[idx]))
# Plot each gene
for i, gene in enumerate(key_genes):
  if gene in chick_analyzer.gene_names:
      ax = axes[i]
      idx = gene2idx[gene]
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
      # Get spatial max value for this gene
      spatial_max = spatial_max_values[gene]
      # Set title with max expression
      ax.set_title(f"{gene}\n(max: {spatial_max:.2f})", fontsize=12, fontweight='bold')
      # Remove axis ticks for cleaner look
      ax.set_xticks([])
      ax.set_yticks([])
      # Add parameter info to bottom right of each subplot
      param_text = f"p={chick_analyzer.params.percentile_clip:.2f}"
      ax.text(0.98, 0.02, param_text,
              transform=ax.transAxes,
              fontsize=8, ha='right', va='bottom',
              bbox=dict(boxstyle='round,pad=0.3',
                       facecolor='white', alpha=0.8))
      # Optionally add axis labels to first subplot only
      if i == 0:
          ax.set_xlabel("Temporal ← NT → Nasal", fontsize=10)
          ax.set_ylabel("Ventral ← DV → Dorsal", fontsize=10)
# Add overall title
# fig.suptitle(f"Spatial Expression Patterns (σ={chick_analyzer.params.smooth_sigma})",
#            fontsize=16, fontweight='bold')
# Adjust layout and save
plt.tight_layout()
plt.subplots_adjust(top=0.93)  # Make room for suptitle
plt.savefig(os.path.join(FIGURES_BASE, "Figure6", "F6A-G_composite_set2.png"), dpi=300, bbox_inches='tight')
plt.close()

# %% tags=["cell-49"]
import matplotlib.pyplot as plt
import numpy as np
# key_genes = ["FGF8", "CYP26C1", "BMP2", "TBX3", "EPHB2", "DUSP6", "FOXG1", "CYP1B1"]
key_genes = ['TBX5', 'TBX3', 'TBX2', 'HMX1', 'SOHO-1', 'FOXG1', 'FOXD1', 'CHRDL1', ] 
# key_genes = ['EPHA7', 'EPHB2', 'SOHO1', 'FOXD1', 'FGF8', 'BMP2', 'CYP26C1', 'CYP1B1'] 
# Create figure with subplots
fig, axes = plt.subplots(2, 4, figsize=(20, 10))
axes = axes.flatten()  # Flatten for easy indexing
# Gene to index mapping
gene2idx = {gene: i for i, gene in enumerate(chick_analyzer.gene_names)}
# Calculate spatial max values for each gene
spatial_max_values = {}
for gene in key_genes:
  if gene in chick_analyzer.gene_names:
      idx = gene2idx[gene]
      spatial_max_values[gene] = float(np.max(chick_analyzer.images[idx]))
# Plot each gene
for i, gene in enumerate(key_genes):
  if gene in chick_analyzer.gene_names:
      ax = axes[i]
      idx = gene2idx[gene]
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
      # Get spatial max value for this gene
      spatial_max = spatial_max_values[gene]
      # Set title with max expression
      ax.set_title(f"{gene}\n(max: {spatial_max:.2f})", fontsize=12, fontweight='bold')
      # Remove axis ticks for cleaner look
      ax.set_xticks([])
      ax.set_yticks([])
      # Add parameter info to bottom right of each subplot
      param_text = f"p={chick_analyzer.params.percentile_clip:.2f}"
      ax.text(0.98, 0.02, param_text,
              transform=ax.transAxes,
              fontsize=8, ha='right', va='bottom',
              bbox=dict(boxstyle='round,pad=0.3',
                       facecolor='white', alpha=0.8))
      # Optionally add axis labels to first subplot only
      if i == 0:
          ax.set_xlabel("Temporal ← NT → Nasal", fontsize=10)
          ax.set_ylabel("Ventral ← DV → Dorsal", fontsize=10)
# Add overall title
# fig.suptitle(f"Spatial Expression Patterns (σ={chick_analyzer.params.smooth_sigma})",
#            fontsize=16, fontweight='bold')
# Adjust layout and save
plt.tight_layout()
plt.subplots_adjust(top=0.93)  # Make room for suptitle
plt.savefig(os.path.join(FIGURES_BASE, "Figure6", "F6A-G_composite_set1.png"), dpi=300, bbox_inches='tight')
plt.close()
