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
# # Figures S20-S22: Fgf/BMP Pathway Maps in Mouse and Human Retinas

# %% [markdown]
# ## SF22B: Human BMP signaling

# %% [markdown]
# ## Figure 94

# %% tags=["cell-187"]
# Self-contained cell for generating BMP signaling genes figure for human
# Creates a multi-panel figure organized by functional categories
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
# Define BMP signaling genes by functional category
BMP_GENES = {
  'Ligands': ['BMP2', 'BMP4', 'BMP7'],
  'Receptors': ['BMPR1A', 'BMPR1B', 'BMPR2', 'ACVR2A', 'ACVR2B'],
  'Antagonists/Modulators': ['CHRD', 'NOG', 'BAMBI', 'TWSG1', 'TSKU', 'GREM1', 'GREM2',
                             'BMPER', 'CHRDL1', 'CHRDL2', 'FST'],
  'Intracellular Effectors': ['SMAD1', 'SMAD5', 'SMAD9', 'SMAD4'],
  'Downstream Targets': ['ID1', 'ID2', 'ID3', 'MSX1', 'MSX2', 'TBX5'],
  'Counter-Marker': ['VAX2']
}
# Flatten the gene list for easy access
ALL_BMP_GENES = []
for category, genes in BMP_GENES.items():
  ALL_BMP_GENES.extend(genes)
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
  human_params = SpatialAnalysisParams(
      bin_size=40,
      min_gene_count=15,
      min_cells_per_pixel=3,
      percentile_clip=0.93,
      smooth_sigma=1.0,
      mask_count_threshold=3
  )
  human_analyzer = SpatialExpressionAnalyzer(human_params)
  if os.path.exists("../data/20250604_human_RPC.h5ad"):
      human_results = human_analyzer.run_full_analysis("../data/20250604_human_RPC.h5ad")
  else:
      from spatial_expression_analysis import import_seurat_export
      human_adata = import_seurat_export("../data/20250604human.RPC", prefix="20250604human.RPC_", parallel=False)
      human_adata.write_h5ad("../data/20250604_human_RPC.h5ad")
      human_results = human_analyzer.run_full_analysis("../data/20250604_human_RPC.h5ad")
  with open(human_pickle, 'wb') as f:
      pickle.dump(human_analyzer, f)
print(f"Loaded {len(human_analyzer.gene_names)} genes")
# Create output directories
FIGURES_BASE = os.path.join(os.path.dirname(__file__), "..", "..", "figures")
output_dir = os.path.join(FIGURES_BASE, "Figure_SF22")
individual_dir = f"{output_dir}/individual_genes"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(individual_dir, exist_ok=True)
# Gene to index mapping
gene2idx = {gene: i for i, gene in enumerate(human_analyzer.gene_names)}
def find_gene_variant(gene_name, gene_list):
  """Find gene in list with case variations."""
  variations = [
      gene_name,
      gene_name.upper(),
      gene_name.lower(),
      gene_name.capitalize(),
      gene_name.replace('L', 'l')  # For CHRDL1/2
  ]
  for variant in variations:
      if variant in gene_list:
          return variant
  return None
def plot_gene(ax, analyzer, gene_name, show_title=True, title_size=8):
  """Plot a single gene's spatial expression."""
  # Find the gene
  found_gene = find_gene_variant(gene_name, analyzer.gene_names)
  if found_gene:
      idx = gene2idx[found_gene]
      img = analyzer.images[idx].copy()
      # Apply mask
      if analyzer.counts is not None:
          mask = (analyzer.counts < analyzer.params.mask_count_threshold)
          img[mask] = np.nan
      # Calculate max intensity
      spatial_max = float(np.nanmax(img))
      # Plot the image
      im = ax.imshow(img,
                    origin='lower',
                    cmap='viridis',
                    aspect='equal',
                    interpolation='nearest')
      # Add midlines
      if hasattr(analyzer, 'dv_mid') and analyzer.dv_mid is not None and \
         hasattr(analyzer, 'nt_mid') and analyzer.nt_mid is not None:
          ax.axhline(y=analyzer.dv_mid, color="white", lw=0.5, alpha=0.5)
          ax.axvline(x=analyzer.nt_mid, color="white", lw=0.5, alpha=0.5)
      if show_title:
          ax.set_title(f"{gene_name}\n(max: {spatial_max:.2f})", fontsize=title_size)
      return True, found_gene, spatial_max
  else:
      # Gene not found
      ax.text(0.5, 0.5, f'{gene_name}\nNot found', ha='center', va='center',
             transform=ax.transAxes, fontsize=7, color='gray')
      if show_title:
          ax.set_title(gene_name, fontsize=title_size)
      return False, None, 0.0
print("\n" + "=" * 60)
print("Creating BMP signaling figure with categories...")
print("=" * 60)
# Create main composite figure - arrange genes in a grid
# Calculate total needed panels
total_genes = len(ALL_BMP_GENES)
ncols = 7  # 7 columns for a good layout
nrows = int(np.ceil(total_genes / ncols))
fig = plt.figure(figsize=(21, nrows * 3))
# Create main grid
gs = gridspec.GridSpec(nrows, ncols, figure=fig,
                     hspace=0.4,  # More space for category labels
                     wspace=0.15,
                     left=0.03, right=0.97,
                     top=0.94, bottom=0.02)
# Track genes found and not found
genes_found = []
genes_not_found = []
# Plot genes by category with visual separation
plot_idx = 0
for category, genes in BMP_GENES.items():
  # Add category label
  start_idx = plot_idx
  for gene in genes:
      row = plot_idx // ncols
      col = plot_idx % ncols
      ax = fig.add_subplot(gs[row, col])
      # Use different color for category title
      if gene == genes[0]:  # First gene in category
          title_color = 'darkred'
          title_weight = 'bold'
      else:
          title_color = 'black'
          title_weight = 'normal'
      found, variant_name, spatial_max = plot_gene(ax, human_analyzer, gene, show_title=False)
      # Custom title with category indication
      if gene == genes[0]:
          ax.set_title(f"[{category}]\n{gene}\n(max: {spatial_max:.2f})" if found else f"[{category}]\n{gene}\nNot found",
                      fontsize=7, color=title_color, weight=title_weight)
      else:
          ax.set_title(f"{gene}\n(max: {spatial_max:.2f})" if found else f"{gene}\nNot found",
                      fontsize=7)
      # Remove axes
      ax.set_xticks([])
      ax.set_yticks([])
      for spine in ax.spines.values():
          spine.set_visible(False)
      if found:
          genes_found.append((gene, variant_name))
          # Save individual image
          individual_fig, individual_ax = plt.subplots(1, 1, figsize=(5, 5))
          plot_gene(individual_ax, human_analyzer, gene, show_title=False)
          individual_ax.set_title(f"{gene} - {category}", fontsize=12, fontweight='bold')
          individual_ax.axis('off')
          individual_path = f"{individual_dir}/{gene}_{category.replace('/', '_')}.png"
          individual_fig.savefig(individual_path, dpi=150, bbox_inches='tight', facecolor='white')
          plt.close(individual_fig)
      else:
          genes_not_found.append(gene)
      plot_idx += 1
# Hide any remaining empty subplots
for idx in range(plot_idx, nrows * ncols):
  row = idx // ncols
  col = idx % ncols
  ax = fig.add_subplot(gs[row, col])
  ax.axis('off')
# Add main title
fig.suptitle('BMP Signaling Pathway Components in Human Retinal RPCs',
           fontsize=16, fontweight='bold', y=0.98)
# Save the composite figure
composite_path = f"{output_dir}/SF22_bmp_signaling.png"
fig.savefig(composite_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nComposite figure saved to: {composite_path}")
plt.show()
# Create a cleaner category-separated version
print("\n" + "=" * 60)
print("Creating category-separated version...")
print("=" * 60)
# Calculate layout for category-separated figure
fig2 = plt.figure(figsize=(20, 24))
# Create custom layout with sections for each category
gs_main = gridspec.GridSpec(6, 1, figure=fig2, hspace=0.3, top=0.96, bottom=0.02)
category_axes = {}
y_offset = 0
for cat_idx, (category, genes) in enumerate(BMP_GENES.items()):
  n_genes = len(genes)
  if n_genes == 0:
      continue
  # Create sub-gridspec for this category
  if n_genes <= 4:
      cat_ncols = n_genes
      cat_nrows = 1
  elif n_genes <= 8:
      cat_ncols = 4
      cat_nrows = 2
  else:
      cat_ncols = 6
      cat_nrows = int(np.ceil(n_genes / cat_ncols))
  # Create sub-grid for this category
  gs_cat = gridspec.GridSpecFromSubplotSpec(cat_nrows, cat_ncols,
                                            subplot_spec=gs_main[cat_idx],
                                            wspace=0.1, hspace=0.2)
  # Add category title
  fig2.text(0.02, 1 - (cat_idx + 0.5) / 6, category,
           fontsize=12, fontweight='bold', color='darkred',
           rotation=90, va='center')
  # Plot genes in this category
  for gene_idx, gene in enumerate(genes):
      row = gene_idx // cat_ncols
      col = gene_idx % cat_ncols
      ax = fig2.add_subplot(gs_cat[row, col])
      found, variant_name, spatial_max = plot_gene(ax, human_analyzer, gene, show_title=True, title_size=8)
      # Remove axes
      ax.set_xticks([])
      ax.set_yticks([])
      for spine in ax.spines.values():
          spine.set_visible(False)
# Add main title
fig2.suptitle('BMP Signaling Pathway Components in Human Retinal RPCs (By Category)',
            fontsize=16, fontweight='bold')
# Save the category-separated figure
category_path = f"{output_dir}/SF22_bmp_by_category.png"
fig2.savefig(category_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Category-separated figure saved to: {category_path}")
plt.show()
# Save summary to text file
summary_path = f"{output_dir}/SF22_bmp_summary.txt"
with open(summary_path, 'w') as f:
  f.write("BMP Signaling Pathway Components Analysis\n")
  f.write("=" * 60 + "\n\n")
  f.write(f"Total genes analyzed: {len(ALL_BMP_GENES)}\n")
  f.write(f"Genes found: {len(genes_found)}\n")
  f.write(f"Genes not found: {len(genes_not_found)}\n\n")
  f.write("Genes by Category:\n")
  f.write("-" * 40 + "\n")
  for category, genes in BMP_GENES.items():
      f.write(f"\n{category} ({len(genes)} genes):\n")
      for gene in genes:
          found = any(g[0] == gene for g in genes_found)
          status = "✓" if found else "✗"
          f.write(f"  {status} {gene}\n")
  f.write("\n" + "=" * 60 + "\n")
  f.write("Analysis parameters:\n")
  f.write(f"  Percentile clip: {human_analyzer.params.percentile_clip:.2f}\n")
  f.write(f"  Smoothing sigma: {human_analyzer.params.smooth_sigma:.1f}\n")
  f.write(f"  Bin size: {human_analyzer.params.bin_size}\n")
  f.write(f"  Min gene count: {human_analyzer.params.min_gene_count}\n")
print(f"\nSummary saved to: {summary_path}")
print("\n" + "=" * 60)
print("Figure generation complete!")
print(f"Main figure: {composite_path}")
print(f"Category figure: {category_path}")
print(f"Individual images: {individual_dir}/")
print("=" * 60)
# Print summary
print("\nSummary:")
print(f"  Genes found: {len(genes_found)}/{len(ALL_BMP_GENES)}")
if genes_not_found:
  print(f"  Genes not found: {', '.join(genes_not_found)}")

# %% [markdown]
# ## SF21B: Human FGF8 downstream

# %% [markdown]
# ## Figure 95

# %% tags=["cell-189"]
# Self-contained cell for generating FGF8 downstream and related genes figure for human
# Creates a multi-panel figure organized by functional categories
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
# Define FGF8-related genes by functional category
FGF8_GENES = {
  'FGF8 Ligand': ['FGF8'],
  'FGF8 Downstream Genes': ['DUSP6', 'SOX2', 'EGR1', 'SNAI1', 'FOS'],  # SNAIL -> SNAI1
  'Engrailed Genes': ['EN1', 'EN2'],
  'Sprouty Genes': ['SPRY1', 'SPRY2']  # Standard gene nomenclature
}
# Flatten the gene list for easy access
ALL_FGF8_GENES = []
for category, genes in FGF8_GENES.items():
  ALL_FGF8_GENES.extend(genes)
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
  human_params = SpatialAnalysisParams(
      bin_size=40,
      min_gene_count=15,
      min_cells_per_pixel=3,
      percentile_clip=0.93,
      smooth_sigma=1.0,
      mask_count_threshold=3
  )
  human_analyzer = SpatialExpressionAnalyzer(human_params)
  if os.path.exists("../data/20250604_human_RPC.h5ad"):
      human_results = human_analyzer.run_full_analysis("../data/20250604_human_RPC.h5ad")
  else:
      from spatial_expression_analysis import import_seurat_export
      human_adata = import_seurat_export("../data/20250604human.RPC", prefix="20250604human.RPC_", parallel=False)
      human_adata.write_h5ad("../data/20250604_human_RPC.h5ad")
      human_results = human_analyzer.run_full_analysis("../data/20250604_human_RPC.h5ad")
  with open(human_pickle, 'wb') as f:
      pickle.dump(human_analyzer, f)
print(f"Loaded {len(human_analyzer.gene_names)} genes")
# Create output directories
output_dir = os.path.join(FIGURES_BASE, "Figure_SF21")
individual_dir = f"{output_dir}/individual_genes"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(individual_dir, exist_ok=True)
# Gene to index mapping
gene2idx = {gene: i for i, gene in enumerate(human_analyzer.gene_names)}
def find_gene_variant(gene_name, gene_list):
  """Find gene in list with case variations and common aliases."""
  # Handle common aliases
  aliases = {
      'SNAIL': 'SNAI1',
      'Sprouty 1': 'SPRY1',
      'Sprouty 2': 'SPRY2',
      'SPROUTY1': 'SPRY1',
      'SPROUTY2': 'SPRY2'
  }
  if gene_name in aliases:
      gene_name = aliases[gene_name]
  variations = [
      gene_name,
      gene_name.upper(),
      gene_name.lower(),
      gene_name.capitalize(),
      gene_name.replace('SPRY', 'Spry'),
      gene_name.replace('SPRY', 'SPROUTY'),
      gene_name.replace('SPRY', 'Sprouty'),
      f"SPROUTY{gene_name[-1]}" if gene_name.startswith('SPRY') else gene_name
  ]
  for variant in variations:
      if variant in gene_list:
          return variant
  return None
def plot_gene(ax, analyzer, gene_name, show_title=True, title_size=8):
  """Plot a single gene's spatial expression."""
  # Find the gene
  found_gene = find_gene_variant(gene_name, analyzer.gene_names)
  if found_gene:
      idx = gene2idx[found_gene]
      img = analyzer.images[idx].copy()
      # Apply mask
      if analyzer.counts is not None:
          mask = (analyzer.counts < analyzer.params.mask_count_threshold)
          img[mask] = np.nan
      # Calculate max intensity
      spatial_max = float(np.nanmax(img))
      # Plot the image
      im = ax.imshow(img,
                    origin='lower',
                    cmap='viridis',
                    aspect='equal',
                    interpolation='nearest')
      # Add midlines
      if hasattr(analyzer, 'dv_mid') and analyzer.dv_mid is not None and \
         hasattr(analyzer, 'nt_mid') and analyzer.nt_mid is not None:
          ax.axhline(y=analyzer.dv_mid, color="white", lw=0.5, alpha=0.5)
          ax.axvline(x=analyzer.nt_mid, color="white", lw=0.5, alpha=0.5)
      if show_title:
          ax.set_title(f"{gene_name}\n(max: {spatial_max:.2f})", fontsize=title_size)
      return True, found_gene, spatial_max
  else:
      # Gene not found
      ax.text(0.5, 0.5, f'{gene_name}\nNot found', ha='center', va='center',
             transform=ax.transAxes, fontsize=7, color='gray')
      if show_title:
          ax.set_title(gene_name, fontsize=title_size)
      return False, None, 0.0
print("\n" + "=" * 60)
print("Creating FGF8 pathway figure with categories...")
print("=" * 60)
# Create main composite figure - arrange genes in a grid
# Calculate total needed panels
total_genes = len(ALL_FGF8_GENES)
ncols = 5  # 5 columns for this smaller set
nrows = int(np.ceil(total_genes / ncols))
fig = plt.figure(figsize=(15, nrows * 3))
# Create main grid
gs = gridspec.GridSpec(nrows, ncols, figure=fig,
                     hspace=0.4,
                     wspace=0.15,
                     left=0.03, right=0.97,
                     top=0.92, bottom=0.02)
# Track genes found and not found
genes_found = []
genes_not_found = []
# Plot genes by category with visual separation
plot_idx = 0
for category, genes in FGF8_GENES.items():
  for gene in genes:
      row = plot_idx // ncols
      col = plot_idx % ncols
      ax = fig.add_subplot(gs[row, col])
      # Use different color for category title
      if gene == genes[0]:  # First gene in category
          title_color = 'darkblue' if 'FGF8' in category else 'darkred'
          title_weight = 'bold'
      else:
          title_color = 'black'
          title_weight = 'normal'
      found, variant_name, spatial_max = plot_gene(ax, human_analyzer, gene, show_title=False)
      # Custom title with category indication
      if gene == genes[0]:
          ax.set_title(f"[{category}]\n{gene}\n(max: {spatial_max:.2f})" if found else f"[{category}]\n{gene}\nNot found",
                      fontsize=8, color=title_color, weight=title_weight)
      else:
          ax.set_title(f"{gene}\n(max: {spatial_max:.2f})" if found else f"{gene}\nNot found",
                      fontsize=8)
      # Remove axes
      ax.set_xticks([])
      ax.set_yticks([])
      for spine in ax.spines.values():
          spine.set_visible(False)
      if found:
          genes_found.append((gene, variant_name))
          # Save individual image
          individual_fig, individual_ax = plt.subplots(1, 1, figsize=(5, 5))
          plot_gene(individual_ax, human_analyzer, gene, show_title=False)
          individual_ax.set_title(f"{gene} - {category}", fontsize=12, fontweight='bold')
          individual_ax.axis('off')
          individual_path = f"{individual_dir}/{gene}_{category.replace(' ', '_')}.png"
          individual_fig.savefig(individual_path, dpi=150, bbox_inches='tight', facecolor='white')
          plt.close(individual_fig)
      else:
          genes_not_found.append(gene)
      plot_idx += 1
# Hide any remaining empty subplots
for idx in range(plot_idx, nrows * ncols):
  row = idx // ncols
  col = idx % ncols
  ax = fig.add_subplot(gs[row, col])
  ax.axis('off')
# Add main title
fig.suptitle('FGF8 Signaling and Related Genes in Human Retinal RPCs',
           fontsize=16, fontweight='bold', y=0.96)
# Save the composite figure
composite_path = f"{output_dir}/SF21_fgf8_pathway.png"
fig.savefig(composite_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nComposite figure saved to: {composite_path}")
plt.show()
# Create a cleaner category-separated version
print("\n" + "=" * 60)
print("Creating category-separated version...")
print("=" * 60)
# Calculate layout for category-separated figure
fig2 = plt.figure(figsize=(16, 12))
# Create custom layout with sections for each category
gs_main = gridspec.GridSpec(4, 1, figure=fig2, hspace=0.3, top=0.94, bottom=0.02)
for cat_idx, (category, genes) in enumerate(FGF8_GENES.items()):
  n_genes = len(genes)
  if n_genes == 0:
      continue
  # Create sub-gridspec for this category
  if n_genes <= 3:
      cat_ncols = n_genes
      cat_nrows = 1
  elif n_genes <= 6:
      cat_ncols = 3
      cat_nrows = 2
  else:
      cat_ncols = 4
      cat_nrows = int(np.ceil(n_genes / cat_ncols))
  # Create sub-grid for this category
  gs_cat = gridspec.GridSpecFromSubplotSpec(cat_nrows, cat_ncols,
                                            subplot_spec=gs_main[cat_idx],
                                            wspace=0.1, hspace=0.2)
  # Add category title
  fig2.text(0.02, 1 - (cat_idx + 0.5) / 4, category,
           fontsize=11, fontweight='bold',
           color='darkblue' if 'FGF8' in category else 'darkred',
           rotation=90, va='center')
  # Plot genes in this category
  for gene_idx, gene in enumerate(genes):
      row = gene_idx // cat_ncols
      col = gene_idx % cat_ncols
      ax = fig2.add_subplot(gs_cat[row, col])
      found, variant_name, spatial_max = plot_gene(ax, human_analyzer, gene, show_title=True, title_size=9)
      # Remove axes
      ax.set_xticks([])
      ax.set_yticks([])
      for spine in ax.spines.values():
          spine.set_visible(False)
# Add main title
fig2.suptitle('FGF8 Signaling and Related Genes in Human Retinal RPCs (By Category)',
            fontsize=16, fontweight='bold')
# Save the category-separated figure
category_path = f"{output_dir}/SF21_fgf8_by_category.png"
fig2.savefig(category_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Category-separated figure saved to: {category_path}")
plt.show()
# Also create a correlation analysis with FGF8
print("\n" + "=" * 60)
print("Creating FGF8 correlation analysis...")
print("=" * 60)
# Check if FGF8 exists and get correlations
fgf8_variant = find_gene_variant('FGF8', human_analyzer.gene_names)
if fgf8_variant:
  # Get correlations with all other genes in the set
  correlations = human_analyzer.get_gene_correlations(fgf8_variant, top_n=len(human_analyzer.gene_names))
  # Filter for our genes of interest
  gene_correlations = {}
  for gene in ALL_FGF8_GENES:
      if gene != 'FGF8':
          variant = find_gene_variant(gene, human_analyzer.gene_names)
          if variant:
              for corr_gene, corr_val in correlations:
                  if corr_gene == variant:
                      gene_correlations[gene] = corr_val
                      break
  # Sort by correlation
  sorted_correlations = sorted(gene_correlations.items(), key=lambda x: x[1], reverse=True)
  # Create correlation plot
  fig3, ax = plt.subplots(figsize=(10, 6))
  genes_list = [g for g, _ in sorted_correlations]
  corr_values = [c for _, c in sorted_correlations]
  bars = ax.barh(range(len(genes_list)), corr_values)
  # Color bars by category
  for i, (gene, _) in enumerate(sorted_correlations):
      for category, cat_genes in FGF8_GENES.items():
          if gene in cat_genes:
              if 'Downstream' in category:
                  bars[i].set_color('steelblue')
              elif 'Engrailed' in category:
                  bars[i].set_color('darkgreen')
              elif 'Sprouty' in category:
                  bars[i].set_color('darkorange')
              break
  ax.set_yticks(range(len(genes_list)))
  ax.set_yticklabels(genes_list)
  ax.set_xlabel('Correlation with FGF8', fontsize=12)
  ax.set_title('Correlation of FGF8-Related Genes with FGF8', fontsize=14, fontweight='bold')
  ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
  ax.grid(axis='x', alpha=0.3)
  # Add legend
  from matplotlib.patches import Patch
  legend_elements = [
      Patch(facecolor='steelblue', label='FGF8 Downstream'),
      Patch(facecolor='darkgreen', label='Engrailed'),
      Patch(facecolor='darkorange', label='Sprouty')
  ]
  ax.legend(handles=legend_elements, loc='lower right')
  plt.tight_layout()
  correlation_path = f"{output_dir}/SF21_fgf8_correlations.png"
  fig3.savefig(correlation_path, dpi=300, bbox_inches='tight', facecolor='white')
  print(f"Correlation plot saved to: {correlation_path}")
  plt.show()
# Save summary to text file
summary_path = f"{output_dir}/SF21_fgf8_summary.txt"
with open(summary_path, 'w') as f:
  f.write("FGF8 Signaling and Related Genes Analysis\n")
  f.write("=" * 60 + "\n\n")
  f.write(f"Total genes analyzed: {len(ALL_FGF8_GENES)}\n")
  f.write(f"Genes found: {len(genes_found)}\n")
  f.write(f"Genes not found: {len(genes_not_found)}\n\n")
  f.write("Genes by Category:\n")
  f.write("-" * 40 + "\n")
  for category, genes in FGF8_GENES.items():
      f.write(f"\n{category} ({len(genes)} genes):\n")
      for gene in genes:
          found = any(g[0] == gene for g in genes_found)
          status = "✓" if found else "✗"
          f.write(f"  {status} {gene}\n")
  if fgf8_variant and gene_correlations:
      f.write("\n" + "-" * 40 + "\n")
      f.write("Correlations with FGF8:\n")
      for gene, corr in sorted_correlations:
          f.write(f"  {gene:10s}: r = {corr:.4f}\n")
  f.write("\n" + "=" * 60 + "\n")
  f.write("Analysis parameters:\n")
  f.write(f"  Percentile clip: {human_analyzer.params.percentile_clip:.2f}\n")
  f.write(f"  Smoothing sigma: {human_analyzer.params.smooth_sigma:.1f}\n")
  f.write(f"  Bin size: {human_analyzer.params.bin_size}\n")
  f.write(f"  Min gene count: {human_analyzer.params.min_gene_count}\n")
print(f"\nSummary saved to: {summary_path}")
print("\n" + "=" * 60)
print("Figure generation complete!")
print(f"Main figure: {composite_path}")
print(f"Category figure: {category_path}")
if fgf8_variant:
  print(f"Correlation plot: {correlation_path}")
print(f"Individual images: {individual_dir}/")
print("=" * 60)
# Print summary
print("\nSummary:")
print(f"  Genes found: {len(genes_found)}/{len(ALL_FGF8_GENES)}")
if genes_not_found:
  print(f"  Genes not found: {', '.join(genes_not_found)}")
if fgf8_variant and gene_correlations:
  print("\nTop correlations with FGF8:")
  for gene, corr in sorted_correlations[:5]:
      print(f"  {gene}: r = {corr:.3f}")

# %% [markdown]
# ## SF20B: Human FGF family

# %% [markdown]
# ## Figure 96

# %% tags=["cell-191"]
# Self-contained cell for generating comprehensive FGF family genes figure for human
# Creates a multi-panel figure for all FGF ligands and receptors
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
# Define FGF family genes by category
# Note: FGF11-14 are actually a subfamily, FGF15 in mouse = FGF19 in humans
FGF_FAMILY = {
  'FGF Ligands (1-10)': ['FGF1', 'FGF2', 'FGF3', 'FGF4', 'FGF5',
                         'FGF6', 'FGF7', 'FGF8', 'FGF9', 'FGF10'],
  'FGF Ligands (11-14)': ['FGF11', 'FGF12', 'FGF13', 'FGF14'],  # FGF11-14 subfamily
  'FGF Ligands (16-23)': ['FGF16', 'FGF17', 'FGF18', 'FGF19',  # FGF19 is human equivalent of mouse FGF15
                          'FGF20', 'FGF21', 'FGF22', 'FGF23'],
  'FGF Receptors': ['FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FGFRL1']
}
# Flatten the gene list for easy access
ALL_FGF_GENES = []
for category, genes in FGF_FAMILY.items():
  ALL_FGF_GENES.extend(genes)
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
  human_params = SpatialAnalysisParams(
      bin_size=40,
      min_gene_count=15,
      min_cells_per_pixel=3,
      percentile_clip=0.93,
      smooth_sigma=1.0,
      mask_count_threshold=3
  )
  human_analyzer = SpatialExpressionAnalyzer(human_params)
  if os.path.exists("../data/20250604_human_RPC.h5ad"):
      human_results = human_analyzer.run_full_analysis("../data/20250604_human_RPC.h5ad")
  else:
      from spatial_expression_analysis import import_seurat_export
      human_adata = import_seurat_export("../data/20250604human.RPC", prefix="20250604human.RPC_", parallel=False)
      human_adata.write_h5ad("../data/20250604_human_RPC.h5ad")
      human_results = human_analyzer.run_full_analysis("../data/20250604_human_RPC.h5ad")
  with open(human_pickle, 'wb') as f:
      pickle.dump(human_analyzer, f)
print(f"Loaded {len(human_analyzer.gene_names)} genes")
# Create output directories
output_dir = os.path.join(FIGURES_BASE, "Figure_SF20")
individual_dir = f"{output_dir}/individual_genes"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(individual_dir, exist_ok=True)
# Gene to index mapping
gene2idx = {gene: i for i, gene in enumerate(human_analyzer.gene_names)}
def find_gene_variant(gene_name, gene_list):
  """Find gene in list with case variations."""
  variations = [
      gene_name,
      gene_name.upper(),
      gene_name.lower(),
      gene_name.capitalize(),
      gene_name.replace('FGF', 'Fgf'),
      gene_name.replace('FGFR', 'Fgfr'),
      gene_name.replace('L1', 'l1')  # For FGFRL1
  ]
  for variant in variations:
      if variant in gene_list:
          return variant
  return None
def plot_gene(ax, analyzer, gene_name, show_title=True, title_size=7):
  """Plot a single gene's spatial expression."""
  # Find the gene
  found_gene = find_gene_variant(gene_name, analyzer.gene_names)
  if found_gene:
      idx = gene2idx[found_gene]
      img = analyzer.images[idx].copy()
      # Apply mask
      if analyzer.counts is not None:
          mask = (analyzer.counts < analyzer.params.mask_count_threshold)
          img[mask] = np.nan
      # Calculate max intensity
      spatial_max = float(np.nanmax(img))
      # Plot the image
      im = ax.imshow(img,
                    origin='lower',
                    cmap='viridis',
                    aspect='equal',
                    interpolation='nearest')
      # Add midlines
      if hasattr(analyzer, 'dv_mid') and analyzer.dv_mid is not None and \
         hasattr(analyzer, 'nt_mid') and analyzer.nt_mid is not None:
          ax.axhline(y=analyzer.dv_mid, color="white", lw=0.5, alpha=0.5)
          ax.axvline(x=analyzer.nt_mid, color="white", lw=0.5, alpha=0.5)
      if show_title:
          ax.set_title(f"{gene_name}\n(max: {spatial_max:.2f})", fontsize=title_size)
      return True, found_gene, spatial_max
  else:
      # Gene not found
      ax.text(0.5, 0.5, f'{gene_name}\nNot found', ha='center', va='center',
             transform=ax.transAxes, fontsize=6, color='gray')
      if show_title:
          ax.set_title(gene_name, fontsize=title_size)
      return False, None, 0.0
print("\n" + "=" * 60)
print("Creating comprehensive FGF family figure...")
print("=" * 60)
# Create main composite figure - arrange genes in a grid
# Calculate total needed panels
total_genes = len(ALL_FGF_GENES)
ncols = 7  # 7 columns for good layout with many genes
nrows = int(np.ceil(total_genes / ncols))
fig = plt.figure(figsize=(21, nrows * 3))
# Create main grid
gs = gridspec.GridSpec(nrows, ncols, figure=fig,
                     hspace=0.4,
                     wspace=0.15,
                     left=0.03, right=0.97,
                     top=0.93, bottom=0.02)
# Track genes found and not found
genes_found = []
genes_not_found = []
genes_with_expression = []  # Track which genes have notable expression
# Plot genes by category
plot_idx = 0
for category, genes in FGF_FAMILY.items():
  for gene in genes:
      row = plot_idx // ncols
      col = plot_idx % ncols
      ax = fig.add_subplot(gs[row, col])
      # Use different colors for different categories
      if gene == genes[0]:  # First gene in category
          if 'Receptor' in category:
              title_color = 'darkgreen'
          elif '1-10' in category:
              title_color = 'darkblue'
          elif '11-14' in category:
              title_color = 'darkorange'
          else:  # 16-23
              title_color = 'darkred'
          title_weight = 'bold'
      else:
          title_color = 'black'
          title_weight = 'normal'
      found, variant_name, spatial_max = plot_gene(ax, human_analyzer, gene, show_title=False)
      # Track genes with notable expression (arbitrary threshold)
      if found and spatial_max > 0.1:
          genes_with_expression.append((gene, spatial_max))
      # Custom title with category indication
      if gene == genes[0]:
          title_text = f"[{category.split('(')[0].strip()}]\n{gene}"
          if found:
              title_text += f"\n(max: {spatial_max:.2f})"
          else:
              title_text += "\nNot found"
          ax.set_title(title_text, fontsize=7, color=title_color, weight=title_weight)
      else:
          ax.set_title(f"{gene}\n(max: {spatial_max:.2f})" if found else f"{gene}\nNot found",
                      fontsize=7)
      # Remove axes
      ax.set_xticks([])
      ax.set_yticks([])
      for spine in ax.spines.values():
          spine.set_visible(False)
      if found:
          genes_found.append((gene, variant_name))
          # Save individual image
          individual_fig, individual_ax = plt.subplots(1, 1, figsize=(5, 5))
          plot_gene(individual_ax, human_analyzer, gene, show_title=False)
          individual_ax.set_title(f"{gene} - {category}", fontsize=12, fontweight='bold')
          individual_ax.axis('off')
          individual_path = f"{individual_dir}/{gene}.png"
          individual_fig.savefig(individual_path, dpi=150, bbox_inches='tight', facecolor='white')
          plt.close(individual_fig)
      else:
          genes_not_found.append(gene)
      plot_idx += 1
# Hide any remaining empty subplots
for idx in range(plot_idx, nrows * ncols):
  row = idx // ncols
  col = idx % ncols
  ax = fig.add_subplot(gs[row, col])
  ax.axis('off')
# Add main title
fig.suptitle('Comprehensive FGF Family Expression in Human Retinal RPCs',
           fontsize=16, fontweight='bold', y=0.96)
# Save the composite figure
composite_path = f"{output_dir}/SF20_fgf_family.png"
fig.savefig(composite_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nComposite figure saved to: {composite_path}")
plt.show()
# Create a focused figure showing only expressed genes
print("\n" + "=" * 60)
print("Creating expressed genes figure...")
print("=" * 60)
# Sort genes by expression level
genes_with_expression.sort(key=lambda x: x[1], reverse=True)
if genes_with_expression:
  n_expressed = len(genes_with_expression)
  expressed_ncols = min(5, n_expressed)
  expressed_nrows = int(np.ceil(n_expressed / expressed_ncols))
  fig2 = plt.figure(figsize=(expressed_ncols * 3, expressed_nrows * 3))
  gs2 = gridspec.GridSpec(expressed_nrows, expressed_ncols, figure=fig2,
                         hspace=0.3, wspace=0.15,
                         left=0.05, right=0.95,
                         top=0.92, bottom=0.05)
  for idx, (gene, max_val) in enumerate(genes_with_expression):
      row = idx // expressed_ncols
      col = idx % expressed_ncols
      ax = fig2.add_subplot(gs2[row, col])
      found, variant_name, spatial_max = plot_gene(ax, human_analyzer, gene, show_title=True, title_size=9)
      # Remove axes
      ax.set_xticks([])
      ax.set_yticks([])
      for spine in ax.spines.values():
          spine.set_visible(False)
  # Hide remaining subplots
  for idx in range(n_expressed, expressed_nrows * expressed_ncols):
      row = idx // expressed_ncols
      col = idx % expressed_ncols
      ax = fig2.add_subplot(gs2[row, col])
      ax.axis('off')
  fig2.suptitle('FGF Family Members with Notable Expression in Human Retinal RPCs',
                fontsize=14, fontweight='bold')
  expressed_path = f"{output_dir}/SF20_fgf_expressed_only.png"
  fig2.savefig(expressed_path, dpi=300, bbox_inches='tight', facecolor='white')
  print(f"Expressed genes figure saved to: {expressed_path}")
  plt.show()
# Create expression heatmap summary
print("\n" + "=" * 60)
print("Creating expression summary heatmap...")
print("=" * 60)
# Collect expression data for all found genes
expression_data = []
gene_names = []
categories = []
for category, genes in FGF_FAMILY.items():
  for gene in genes:
      variant = find_gene_variant(gene, human_analyzer.gene_names)
      if variant:
          idx = gene2idx[variant]
          img = human_analyzer.images[idx]
          max_expr = float(np.nanmax(img))
          mean_expr = float(np.nanmean(img[~np.isnan(img)]))
          expression_data.append([max_expr, mean_expr])
          gene_names.append(gene)
          categories.append(category.split('(')[0].strip())
if expression_data:
  # Create heatmap
  fig3, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, len(gene_names) * 0.3))
  # Max expression heatmap
  max_data = np.array([d[0] for d in expression_data]).reshape(-1, 1)
  im1 = ax1.imshow(max_data, cmap='YlOrRd', aspect='auto')
  ax1.set_yticks(range(len(gene_names)))
  ax1.set_yticklabels(gene_names, fontsize=8)
  ax1.set_xticks([0])
  ax1.set_xticklabels(['Max Expression'])
  ax1.set_title('Maximum Expression', fontsize=10, fontweight='bold')
  # Add values
  for i, val in enumerate(max_data):
      ax1.text(0, i, f'{val[0]:.2f}', ha='center', va='center', fontsize=7)
  # Mean expression heatmap
  mean_data = np.array([d[1] for d in expression_data]).reshape(-1, 1)
  im2 = ax2.imshow(mean_data, cmap='YlOrRd', aspect='auto')
  ax2.set_yticks(range(len(gene_names)))
  ax2.set_yticklabels(gene_names, fontsize=8)
  ax2.set_xticks([0])
  ax2.set_xticklabels(['Mean Expression'])
  ax2.set_title('Mean Expression', fontsize=10, fontweight='bold')
  # Add values
  for i, val in enumerate(mean_data):
      ax2.text(0, i, f'{val[0]:.3f}', ha='center', va='center', fontsize=7)
  fig3.suptitle('FGF Family Expression Summary', fontsize=12, fontweight='bold')
  plt.tight_layout()
  heatmap_path = f"{output_dir}/SF20_fgf_expression_summary.png"
  fig3.savefig(heatmap_path, dpi=300, bbox_inches='tight', facecolor='white')
  print(f"Expression heatmap saved to: {heatmap_path}")
  plt.show()
# Print final summary
print("\n" + "=" * 60)
print("FGF Family Analysis Complete!")
print(f"  Genes found: {len(genes_found)}/{len(ALL_FGF_GENES)}")
print(f"  Genes with notable expression: {len(genes_with_expression)}")
if genes_not_found:
  print(f"  Genes not found: {', '.join(genes_not_found)}")
print(f"Output directory: {output_dir}")
print("=" * 60)
