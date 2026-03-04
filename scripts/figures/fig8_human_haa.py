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
# # Figure 8B-C: Human HAA Gene Correlations

# %% [markdown]
# ## F8B: Cyp26a1 top-correlated genes

# %% tags=["cell-183"]
# Self-contained cell for generating CYP26C1 correlation figure for human
# Creates a 4x3 tiled figure showing CYP26C1 and its top 11 correlated genes
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
import pickle
import warnings
warnings.filterwarnings('ignore')
# Add parent directory to path for imports
sys.path.append('..')
from spatial_expression_analysis import (
  SpatialExpressionAnalyzer,
  SpatialAnalysisParams
)
print("Loading human analyzer...")
print("=" * 60)
# Load Human analyzer with correct parameters from notebook
human_pickle = "human_analyzer_correct.pkl"
if os.path.exists(human_pickle):
  print("Loading from pickle file...")
  with open(human_pickle, 'rb') as f:
      human_analyzer = pickle.load(f)
else:
  print("Creating new analyzer...")
  # Use the exact parameters from the notebook
  human_params = SpatialAnalysisParams(
      bin_size=40,
      min_gene_count=15,  # This should keep CYP26C1
      min_cells_per_pixel=3,
      percentile_clip=0.93,
      smooth_sigma=1.0,
      mask_count_threshold=3
  )
  human_analyzer = SpatialExpressionAnalyzer(human_params)
  # Try to load from h5ad or import from export
  if os.path.exists("../data/20250604_human_RPC.h5ad"):
      human_results = human_analyzer.run_full_analysis("../data/20250604_human_RPC.h5ad")
  else:
      from spatial_expression_analysis import import_seurat_export
      human_adata = import_seurat_export("../data/20250604human.RPC", prefix="20250604human.RPC_", parallel=False)
      human_adata.write_h5ad("../data/20250604_human_RPC.h5ad")
      human_results = human_analyzer.run_full_analysis("../data/20250604_human_RPC.h5ad")
  # Save for future use
  with open(human_pickle, 'wb') as f:
      pickle.dump(human_analyzer, f)
print(f"Loaded {len(human_analyzer.gene_names)} genes")
# Check if CYP26C1 exists
target_gene = "CYP26A1"
if target_gene not in human_analyzer.gene_names:
  # Try different case variations
  for variant in ["Cyp26a1", "cyp26a1", "CYP26a1"]:
      if variant in human_analyzer.gene_names:
          target_gene = variant
          break
  else:
      print(f"ERROR: {target_gene} not found in human data!")
      print("Available genes starting with CYP:")
      for gene in sorted(human_analyzer.gene_names):
          if gene.upper().startswith("CYP"):
              print(f"  {gene}")
      raise ValueError(f"{target_gene} not found")
print(f"\nAnalyzing {target_gene}...")
# Get top correlated genes
correlations = human_analyzer.get_gene_correlations(target_gene, top_n=12)  # Get 12 total (main + 11 others)
print(f"\nTop correlations for {target_gene}:")
for i, (gene, corr) in enumerate(correlations):
  print(f"  {i+1}. {gene}: {corr:.3f}")
# Create output directory
FIGURES_BASE = os.path.join(os.path.dirname(__file__), "..", "..", "figures")
output_dir = os.path.join(FIGURES_BASE, "Figure8")
os.makedirs(output_dir, exist_ok=True)
print("\n" + "=" * 60)
print("Creating correlation figure...")
print("=" * 60)
# Create 4x3 figure (12 panels total)
fig = plt.figure(figsize=(12, 9))
# Create grid with minimal spacing
gs = gridspec.GridSpec(3, 4, figure=fig,
                     hspace=0.3,  # Some vertical spacing for titles
                     wspace=0.1,  # Minimal horizontal spacing
                     left=0.05, right=0.95,
                     top=0.92, bottom=0.05)
# Gene to index mapping
gene2idx = {gene: i for i, gene in enumerate(human_analyzer.gene_names)}
# Plot each gene
for plot_idx, (gene, corr) in enumerate(correlations[:12]):
  # Calculate position in grid
  row = plot_idx // 4
  col = plot_idx % 4
  ax = fig.add_subplot(gs[row, col])
  # Get gene image
  idx = gene2idx[gene]
  img = human_analyzer.images[idx].copy()
  # Apply mask if available
  if human_analyzer.counts is not None:
      mask = (human_analyzer.counts < human_analyzer.params.mask_count_threshold)
      img[mask] = np.nan
  # Calculate max intensity
  spatial_max = float(np.nanmax(img))
  # Plot the image with correct origin
  im = ax.imshow(img,
                 origin='lower',  # Important: match notebook orientation
                 cmap='viridis',
                 aspect='equal',
                 interpolation='nearest')
  # Add midlines if available
  if hasattr(human_analyzer, 'dv_mid') and human_analyzer.dv_mid is not None and \
     hasattr(human_analyzer, 'nt_mid') and human_analyzer.nt_mid is not None:
      ax.axhline(y=human_analyzer.dv_mid, color="white", lw=0.5, alpha=0.5)
      ax.axvline(x=human_analyzer.nt_mid, color="white", lw=0.5, alpha=0.5)
  # Format title
  if plot_idx == 0:
      # Main gene - different formatting
      ax.set_title(f"{gene}\n(max: {spatial_max:.2f})",
                  fontsize=10, fontweight='bold', color='darkred')
  else:
      # Correlated genes - show correlation value
      ax.set_title(f"{gene}\nr={corr:.3f}, max: {spatial_max:.2f}",
                  fontsize=9)
  # Remove axes
  ax.set_xticks([])
  ax.set_yticks([])
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.spines['bottom'].set_visible(False)
  ax.spines['left'].set_visible(False)
# Add main title
fig.suptitle(f'Human {target_gene} and Top 11 Correlated Genes',
           fontsize=14, fontweight='bold', y=0.98)
# Add parameter caption at bottom
param_text = (f"Parameters: percentile={human_analyzer.params.percentile_clip:.2f}, "
            f"smooth={human_analyzer.params.smooth_sigma:.1f}, "
            f"bin_size={human_analyzer.params.bin_size}")
fig.text(0.5, 0.01, param_text, ha='center', fontsize=8, style='italic', color='gray')
# Save the figure
output_path = f"{output_dir}/F8B_CYP26A1.png"
fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nFigure saved to: {output_path}")
plt.show()
# Also save correlation data to text file
corr_output_path = f"{output_dir}/F8B_CYP26A1.txt"
with open(corr_output_path, 'w') as f:
  f.write(f"Top correlations for human {target_gene}:\n")
  f.write("=" * 50 + "\n")
  for i, (gene, corr) in enumerate(correlations[:12]):
      if i == 0:
          f.write(f"{i+1:2d}. {gene:15s}: r = {corr:.4f} (MAIN GENE)\n")
      else:
          f.write(f"{i+1:2d}. {gene:15s}: r = {corr:.4f}\n")
  f.write("\n" + "=" * 50 + "\n")
  f.write(f"Analysis parameters:\n")
  f.write(f"  Percentile clip: {human_analyzer.params.percentile_clip:.2f}\n")
  f.write(f"  Smoothing sigma: {human_analyzer.params.smooth_sigma:.1f}\n")
  f.write(f"  Bin size: {human_analyzer.params.bin_size}\n")
  f.write(f"  Min gene count: {human_analyzer.params.min_gene_count}\n")
print(f"Correlation data saved to: {corr_output_path}")
print("\n" + "=" * 60)
print("Figure generation complete!")
print("=" * 60)

# %% [markdown]
# ## F8C: Cyp26c1 top-correlated genes

# %% tags=["cell-185"]
# Self-contained cell for generating CYP26C1 correlation figure for human
# Creates a 4x3 tiled figure showing CYP26C1 and its top 11 correlated genes
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
import pickle
import warnings
warnings.filterwarnings('ignore')
# Add parent directory to path for imports
sys.path.append('..')
from spatial_expression_analysis import (
  SpatialExpressionAnalyzer,
  SpatialAnalysisParams
)
print("Loading human analyzer...")
print("=" * 60)
# Load Human analyzer with correct parameters from notebook
human_pickle = "human_analyzer_correct.pkl"
if os.path.exists(human_pickle):
  print("Loading from pickle file...")
  with open(human_pickle, 'rb') as f:
      human_analyzer = pickle.load(f)
else:
  print("Creating new analyzer...")
  # Use the exact parameters from the notebook
  human_params = SpatialAnalysisParams(
      bin_size=40,
      min_gene_count=15,  # This should keep CYP26C1
      min_cells_per_pixel=3,
      percentile_clip=0.93,
      smooth_sigma=1.0,
      mask_count_threshold=3
  )
  human_analyzer = SpatialExpressionAnalyzer(human_params)
  # Try to load from h5ad or import from export
  if os.path.exists("../data/20250604_human_RPC.h5ad"):
      human_results = human_analyzer.run_full_analysis("../data/20250604_human_RPC.h5ad")
  else:
      from spatial_expression_analysis import import_seurat_export
      human_adata = import_seurat_export("../data/20250604human.RPC", prefix="20250604human.RPC_", parallel=False)
      human_adata.write_h5ad("../data/20250604_human_RPC.h5ad")
      human_results = human_analyzer.run_full_analysis("../data/20250604_human_RPC.h5ad")
  # Save for future use
  with open(human_pickle, 'wb') as f:
      pickle.dump(human_analyzer, f)
print(f"Loaded {len(human_analyzer.gene_names)} genes")
# Check if CYP26A1 exists
target_gene = "CYP26C1"
if target_gene not in human_analyzer.gene_names:
  # Try different case variations
  for variant in ["Cyp26c1", "cyp26c1", "CYP26c1"]:
      if variant in human_analyzer.gene_names:
          target_gene = variant
          break
  else:
      print(f"ERROR: {target_gene} not found in human data!")
      print("Available genes starting with CYP:")
      for gene in sorted(human_analyzer.gene_names):
          if gene.upper().startswith("CYP"):
              print(f"  {gene}")
      raise ValueError(f"{target_gene} not found")
print(f"\nAnalyzing {target_gene}...")
# Get top correlated genes
correlations = human_analyzer.get_gene_correlations(target_gene, top_n=12)  # Get 12 total (main + 11 others)
print(f"\nTop correlations for {target_gene}:")
for i, (gene, corr) in enumerate(correlations):
  print(f"  {i+1}. {gene}: {corr:.3f}")
# Create output directory
output_dir = os.path.join(FIGURES_BASE, "Figure8")
os.makedirs(output_dir, exist_ok=True)
print("\n" + "=" * 60)
print("Creating correlation figure...")
print("=" * 60)
# Create 4x3 figure (12 panels total)
fig = plt.figure(figsize=(12, 9))
# Create grid with minimal spacing
gs = gridspec.GridSpec(3, 4, figure=fig,
                     hspace=0.3,  # Some vertical spacing for titles
                     wspace=0.1,  # Minimal horizontal spacing
                     left=0.05, right=0.95,
                     top=0.92, bottom=0.05)
# Gene to index mapping
gene2idx = {gene: i for i, gene in enumerate(human_analyzer.gene_names)}
# Plot each gene
for plot_idx, (gene, corr) in enumerate(correlations[:12]):
  # Calculate position in grid
  row = plot_idx // 4
  col = plot_idx % 4
  ax = fig.add_subplot(gs[row, col])
  # Get gene image
  idx = gene2idx[gene]
  img = human_analyzer.images[idx].copy()
  # Apply mask if available
  if human_analyzer.counts is not None:
      mask = (human_analyzer.counts < human_analyzer.params.mask_count_threshold)
      img[mask] = np.nan
  # Calculate max intensity
  spatial_max = float(np.nanmax(img))
  # Plot the image with correct origin
  im = ax.imshow(img,
                 origin='lower',  # Important: match notebook orientation
                 cmap='viridis',
                 aspect='equal',
                 interpolation='nearest')
  # Add midlines if available
  if hasattr(human_analyzer, 'dv_mid') and human_analyzer.dv_mid is not None and \
     hasattr(human_analyzer, 'nt_mid') and human_analyzer.nt_mid is not None:
      ax.axhline(y=human_analyzer.dv_mid, color="white", lw=0.5, alpha=0.5)
      ax.axvline(x=human_analyzer.nt_mid, color="white", lw=0.5, alpha=0.5)
  # Format title
  if plot_idx == 0:
      # Main gene - different formatting
      ax.set_title(f"{gene}\n(max: {spatial_max:.2f})",
                  fontsize=10, fontweight='bold', color='darkred')
  else:
      # Correlated genes - show correlation value
      ax.set_title(f"{gene}\nr={corr:.3f}, max: {spatial_max:.2f}",
                  fontsize=9)
  # Remove axes
  ax.set_xticks([])
  ax.set_yticks([])
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.spines['bottom'].set_visible(False)
  ax.spines['left'].set_visible(False)
# Add main title
fig.suptitle(f'Human {target_gene} and Top 11 Correlated Genes',
           fontsize=14, fontweight='bold', y=0.98)
# Add parameter caption at bottom
param_text = (f"Parameters: percentile={human_analyzer.params.percentile_clip:.2f}, "
            f"smooth={human_analyzer.params.smooth_sigma:.1f}, "
            f"bin_size={human_analyzer.params.bin_size}")
fig.text(0.5, 0.01, param_text, ha='center', fontsize=8, style='italic', color='gray')
# Save the figure
output_path = f"{output_dir}/F8C_CYP26C1.png"
fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nFigure saved to: {output_path}")
plt.show()
# Also save correlation data to text file
corr_output_path = f"{output_dir}/F8C_CYP26C1.txt"
with open(corr_output_path, 'w') as f:
  f.write(f"Top correlations for human {target_gene}:\n")
  f.write("=" * 50 + "\n")
  for i, (gene, corr) in enumerate(correlations[:12]):
      if i == 0:
          f.write(f"{i+1:2d}. {gene:15s}: r = {corr:.4f} (MAIN GENE)\n")
      else:
          f.write(f"{i+1:2d}. {gene:15s}: r = {corr:.4f}\n")
  f.write("\n" + "=" * 50 + "\n")
  f.write(f"Analysis parameters:\n")
  f.write(f"  Percentile clip: {human_analyzer.params.percentile_clip:.2f}\n")
  f.write(f"  Smoothing sigma: {human_analyzer.params.smooth_sigma:.1f}\n")
  f.write(f"  Bin size: {human_analyzer.params.bin_size}\n")
  f.write(f"  Min gene count: {human_analyzer.params.min_gene_count}\n")
print(f"Correlation data saved to: {corr_output_path}")
print("\n" + "=" * 60)
print("Figure generation complete!")
print("=" * 60)
