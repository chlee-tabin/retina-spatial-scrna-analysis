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
# # Figures S20A-S22A: Fgf/BMP Pathway Maps in the Mouse Retina
# Uses data/20260528_mouse_RPC_cr9_e13e16.h5ad from GEO GSE322831 and the
# published Figure 7 mouse parameters via the shared analyzer helper.

import matplotlib; matplotlib.use("Agg")

# %% [markdown]
# ## SF20A: FGF family

# %%
# Self-contained cell for generating comprehensive FGF family genes figure for MOUSE
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

# Put the repo root on sys.path for imports.
try:
    REPO = Path(__file__).resolve().parents[2]
except NameError:  # jupytext notebook: start Jupyter from the repo root
    REPO = Path.cwd()
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from scripts.realign_mouse.build_mouse_pathway_pickle import build_analyzer
FIGURES_BASE = os.path.join(str(REPO), "figures")

from spatial_expression_analysis import (
  SpatialExpressionAnalyzer,
  SpatialAnalysisParams
)

# Define FGF family genes by category (mouse nomenclature)
# Note: FGF15 in mouse = FGF19 in humans
FGF_FAMILY = {
  'FGF Ligands (1-10)': ['Fgf1', 'Fgf2', 'Fgf3', 'Fgf4', 'Fgf5',
                         'Fgf6', 'Fgf7', 'Fgf8', 'Fgf9', 'Fgf10'],
  'FGF Ligands (11-14)': ['Fgf11', 'Fgf12', 'Fgf13', 'Fgf14'],  # FGF11-14 subfamily
  'FGF Ligands (15-23)': ['Fgf15', 'Fgf16', 'Fgf17', 'Fgf18',  # FGF15 is mouse-specific
                          'Fgf20', 'Fgf21', 'Fgf22', 'Fgf23'],
  'FGF Receptors': ['Fgfr1', 'Fgfr2', 'Fgfr3', 'Fgfr4', 'Fgfrl1']
}

# Flatten the gene list
ALL_FGF_GENES = []
for category, genes in FGF_FAMILY.items():
  ALL_FGF_GENES.extend(genes)

print("Loading mouse analyzer...")
print("=" * 60)

# Load or build from the deposited CR9 E13-E16 object with published parameters.
mouse_analyzer = build_analyzer()

print(f"Loaded {len(mouse_analyzer.gene_names)} genes")

# Create output directories
output_dir = os.path.join(FIGURES_BASE, "Figure_SF20")
individual_dir = f"{output_dir}/SF20A_individual_genes"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(individual_dir, exist_ok=True)

# Gene to index mapping
gene2idx = {gene: i for i, gene in enumerate(mouse_analyzer.gene_names)}

def find_gene_variant(gene_name, gene_list):
  """Find gene in list with case variations for mouse genes."""
  variations = [
      gene_name,
      gene_name.upper(),
      gene_name.lower(),
      gene_name.capitalize(),
      gene_name[0].upper() + gene_name[1:].lower() if len(gene_name) > 1 else gene_name
  ]

  # Special handling for Fgfrl1
  if 'fgfrl' in gene_name.lower():
      variations.extend(['Fgfrl1', 'FGFRL1', 'FgfrL1'])

  for variant in variations:
      if variant in gene_list:
          return variant
  return None

def plot_gene(ax, analyzer, gene_name, show_title=True, title_size=7):
  """Plot a single gene's spatial expression."""
  found_gene = find_gene_variant(gene_name, analyzer.gene_names)

  if found_gene:
      idx = gene2idx[found_gene]
      img = analyzer.images[idx].copy()

      if analyzer.counts is not None:
          mask = (analyzer.counts < analyzer.params.mask_count_threshold)
          img[mask] = np.nan

      spatial_max = float(np.nanmax(img))

      im = ax.imshow(img,
                    origin='lower',
                    cmap='viridis',
                    aspect='equal',
                    interpolation='nearest')

      if hasattr(analyzer, 'dv_mid') and analyzer.dv_mid is not None and \
         hasattr(analyzer, 'nt_mid') and analyzer.nt_mid is not None:
          ax.axhline(y=analyzer.dv_mid, color="white", lw=0.5, alpha=0.5)
          ax.axvline(x=analyzer.nt_mid, color="white", lw=0.5, alpha=0.5)

      if show_title:
          ax.set_title(f"{gene_name}\n(max: {spatial_max:.2f})", fontsize=title_size)

      return True, found_gene, spatial_max
  else:
      ax.text(0.5, 0.5, f'{gene_name}\nNot found', ha='center', va='center',
             transform=ax.transAxes, fontsize=6, color='gray')
      if show_title:
          ax.set_title(gene_name, fontsize=title_size)
      return False, None, 0.0

print("\n" + "=" * 60)
print("Creating comprehensive FGF family figure for mouse...")
print("=" * 60)

# Create main composite figure
total_genes = len(ALL_FGF_GENES)
ncols = 7  # 7 columns for good layout
nrows = int(np.ceil(total_genes / ncols))

fig = plt.figure(figsize=(21, nrows * 3))

gs = gridspec.GridSpec(nrows, ncols, figure=fig,
                     hspace=0.4,
                     wspace=0.15,
                     left=0.03, right=0.97,
                     top=0.93, bottom=0.02)

# Track genes
genes_found = []
genes_not_found = []
genes_with_expression = []

# Plot genes by category
plot_idx = 0
for category, genes in FGF_FAMILY.items():
  for gene in genes:
      row = plot_idx // ncols
      col = plot_idx % ncols

      ax = fig.add_subplot(gs[row, col])

      # Color coding for categories
      if gene == genes[0]:  # First gene in category
          if 'Receptor' in category:
              title_color = 'darkgreen'
          elif '1-10' in category:
              title_color = 'darkblue'
          elif '11-14' in category:
              title_color = 'darkorange'
          else:  # 15-23
              title_color = 'darkred'
          title_weight = 'bold'
      else:
          title_color = 'black'
          title_weight = 'normal'

      found, variant_name, spatial_max = plot_gene(ax, mouse_analyzer, gene, show_title=False)

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

      ax.set_xticks([])
      ax.set_yticks([])
      for spine in ax.spines.values():
          spine.set_visible(False)

      if found:
          genes_found.append((gene, variant_name))

          # Save individual image
          individual_fig, individual_ax = plt.subplots(1, 1, figsize=(5, 5))
          plot_gene(individual_ax, mouse_analyzer, gene, show_title=False)
          individual_ax.set_title(f"{gene} - {category} (Mouse)", fontsize=12, fontweight='bold')
          individual_ax.axis('off')

          individual_path = f"{individual_dir}/{gene}.png"
          individual_fig.savefig(individual_path, dpi=150, bbox_inches='tight', facecolor='white')
          plt.close(individual_fig)
      else:
          genes_not_found.append(gene)

      plot_idx += 1

# Hide remaining empty subplots
for idx in range(plot_idx, nrows * ncols):
  row = idx // ncols
  col = idx % ncols
  ax = fig.add_subplot(gs[row, col])
  ax.axis('off')

# figure-level suptitle removed (titleless-figure standard; caption -> LEGEND.md;
# the baked-in header overlapped the top-row panel titles on the short grids)

composite_path = f"{output_dir}/SF20A_fgf_family.png"
fig.savefig(composite_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nMain figure saved to: {composite_path}")

plt.show()

# Create a focused figure showing only expressed genes
if genes_with_expression:
  print("\n" + "=" * 60)
  print("Creating expressed genes figure...")
  print("=" * 60)

  genes_with_expression.sort(key=lambda x: x[1], reverse=True)

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

      found, variant_name, spatial_max = plot_gene(ax, mouse_analyzer, gene, show_title=True, title_size=9)

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

  fig2.suptitle('FGF Family Members with Notable Expression in Mouse Retinal RPCs',
                fontsize=14, fontweight='bold')

  expressed_path = f"{output_dir}/SF20A_fgf_family_expressed_only.png"
  fig2.savefig(expressed_path, dpi=300, bbox_inches='tight', facecolor='white')
  print(f"Expressed genes figure saved to: {expressed_path}")

  plt.show()

# Save summary
summary_path = f"{output_dir}/SF20A_fgf_family_summary.txt"
with open(summary_path, 'w') as f:
  f.write("Comprehensive FGF Family Analysis - MOUSE\n")
  f.write("=" * 60 + "\n\n")

  f.write(f"Total genes analyzed: {len(ALL_FGF_GENES)}\n")
  f.write(f"Genes found: {len(genes_found)}\n")
  f.write(f"Genes not found: {len(genes_not_found)}\n")
  f.write(f"Genes with expression >0.1: {len(genes_with_expression)}\n\n")

  for category, genes in FGF_FAMILY.items():
      f.write(f"\n{category} ({len(genes)} genes):\n")
      for gene in genes:
          found = any(g[0] == gene for g in genes_found)
          if found:
              variant = find_gene_variant(gene, mouse_analyzer.gene_names)
              if variant:
                  idx = gene2idx[variant]
                  img = mouse_analyzer.images[idx]
                  max_expr = float(np.nanmax(img))
                  f.write(f"  ✓ {gene:8s} (max: {max_expr:.3f})\n")
          else:
              f.write(f"  ✗ {gene:8s}\n")

  if genes_with_expression:
      f.write("\n" + "-" * 40 + "\n")
      f.write("Top Expressed FGF Family Members:\n")
      for i, (gene, max_val) in enumerate(genes_with_expression[:10]):
          f.write(f"  {i+1:2d}. {gene:8s}: max = {max_val:.3f}\n")

  f.write("\n" + "=" * 60 + "\n")
  f.write("Analysis parameters:\n")
  f.write(f"  Percentile clip: {mouse_analyzer.params.percentile_clip:.2f}\n")
  f.write(f"  Smoothing sigma: {mouse_analyzer.params.smooth_sigma:.1f}\n")
  f.write(f"  Bin size: {mouse_analyzer.params.bin_size}\n")
  f.write(f"  Min gene count: {mouse_analyzer.params.min_gene_count}\n")

  f.write("\n" + "-" * 40 + "\n")
  f.write("Notes:\n")
  f.write("- FGF15 in mouse corresponds to FGF19 in humans\n")
  f.write("- FGF11-14 form a distinct subfamily with intracellular functions\n")
  f.write("- Fgfrl1 is FGF receptor-like 1, lacks tyrosine kinase domain\n")

print(f"Summary saved to: {summary_path}")

print("\n" + "=" * 60)
print("Figure generation complete!")
print(f"  Main figure: {composite_path}")
if genes_with_expression:
  print(f"  Expressed genes figure: {expressed_path}")
print(f"  Individual images: {individual_dir}/")
print(f"  Summary: {summary_path}")
print("=" * 60)

# Print summary
print(f"\nTotal FGF genes: {len(ALL_FGF_GENES)}")
print(f"Genes found: {len(genes_found)}/{len(ALL_FGF_GENES)}")
print(f"Genes with notable expression: {len(genes_with_expression)}")

if genes_with_expression:
  print("\nTop 5 expressed FGF family members:")
  for gene, max_val in genes_with_expression[:5]:
      print(f"  {gene}: max = {max_val:.3f}")

if genes_not_found:
  print(f"\nGenes not found ({len(genes_not_found)} total):")
  for gene in genes_not_found[:10]:
      print(f"  {gene}")
  if len(genes_not_found) > 10:
      print(f"  ... and {len(genes_not_found) - 10} more")



# %% [markdown]
# ## SF21A: FGF8 downstream

# %%
# Self-contained cell for generating FGF8 downstream and related genes figure for MOUSE
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

# Put the repo root on sys.path for imports.
try:
    REPO = Path(__file__).resolve().parents[2]
except NameError:  # jupytext notebook: start Jupyter from the repo root
    REPO = Path.cwd()
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from scripts.realign_mouse.build_mouse_pathway_pickle import build_analyzer
FIGURES_BASE = os.path.join(str(REPO), "figures")

from spatial_expression_analysis import (
  SpatialExpressionAnalyzer,
  SpatialAnalysisParams
)

# Define FGF8-related genes by functional category (mouse nomenclature)
FGF8_GENES = {
  'FGF8 Ligand': ['Fgf8'],
  'FGF8 Downstream Genes': ['Dusp6', 'Sox2', 'Egr1', 'Snai1', 'Fos'],  # SNAIL -> Snai1
  'Engrailed Genes': ['En1', 'En2'],
  'Sprouty Genes': ['Spry1', 'Spry2']  # Mouse nomenclature
}

# Flatten the gene list
ALL_FGF8_GENES = []
for category, genes in FGF8_GENES.items():
  ALL_FGF8_GENES.extend(genes)

print("Loading mouse analyzer...")
print("=" * 60)

# Load or build from the deposited CR9 E13-E16 object with published parameters.
mouse_analyzer = build_analyzer()

print(f"Loaded {len(mouse_analyzer.gene_names)} genes")

# Create output directories
output_dir = os.path.join(FIGURES_BASE, "Figure_SF21")
individual_dir = f"{output_dir}/SF21A_individual_genes"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(individual_dir, exist_ok=True)

# Gene to index mapping
gene2idx = {gene: i for i, gene in enumerate(mouse_analyzer.gene_names)}

def find_gene_variant(gene_name, gene_list):
  """Find gene in list with case variations for mouse genes."""
  # Handle common aliases
  aliases = {
      'SNAIL': 'Snai1',
      'Snail': 'Snai1',
      'Sprouty 1': 'Spry1',
      'Sprouty 2': 'Spry2',
      'SPROUTY1': 'Spry1',
      'SPROUTY2': 'Spry2'
  }

  if gene_name in aliases:
      gene_name = aliases[gene_name]

  variations = [
      gene_name,
      gene_name.upper(),
      gene_name.lower(),
      gene_name.capitalize(),
      gene_name[0].upper() + gene_name[1:].lower() if len(gene_name) > 1 else gene_name
  ]

  # Special variations for Sprouty genes
  if 'Spry' in gene_name:
      variations.extend(['Sprouty' + gene_name[-1], 'SPROUTY' + gene_name[-1]])

  for variant in variations:
      if variant in gene_list:
          return variant
  return None

def plot_gene(ax, analyzer, gene_name, show_title=True, title_size=8):
  """Plot a single gene's spatial expression."""
  found_gene = find_gene_variant(gene_name, analyzer.gene_names)

  if found_gene:
      idx = gene2idx[found_gene]
      img = analyzer.images[idx].copy()

      if analyzer.counts is not None:
          mask = (analyzer.counts < analyzer.params.mask_count_threshold)
          img[mask] = np.nan

      spatial_max = float(np.nanmax(img))

      im = ax.imshow(img,
                    origin='lower',
                    cmap='viridis',
                    aspect='equal',
                    interpolation='nearest')

      if hasattr(analyzer, 'dv_mid') and analyzer.dv_mid is not None and \
         hasattr(analyzer, 'nt_mid') and analyzer.nt_mid is not None:
          ax.axhline(y=analyzer.dv_mid, color="white", lw=0.5, alpha=0.5)
          ax.axvline(x=analyzer.nt_mid, color="white", lw=0.5, alpha=0.5)

      if show_title:
          ax.set_title(f"{gene_name}\n(max: {spatial_max:.2f})", fontsize=title_size)

      return True, found_gene, spatial_max
  else:
      ax.text(0.5, 0.5, f'{gene_name}\nNot found', ha='center', va='center',
             transform=ax.transAxes, fontsize=7, color='gray')
      if show_title:
          ax.set_title(gene_name, fontsize=title_size)
      return False, None, 0.0

print("\n" + "=" * 60)
print("Creating Mouse FGF8 pathway figure...")
print("=" * 60)

# Create main composite figure
total_genes = len(ALL_FGF8_GENES)
ncols = 5  # 5 columns for this smaller set
nrows = int(np.ceil(total_genes / ncols))

fig = plt.figure(figsize=(15, nrows * 3))

gs = gridspec.GridSpec(nrows, ncols, figure=fig,
                     hspace=0.4,
                     wspace=0.15,
                     left=0.03, right=0.97,
                     top=0.92, bottom=0.02)

# Track genes
genes_found = []
genes_not_found = []
genes_with_expression = []

# Plot genes by category
plot_idx = 0
for category, genes in FGF8_GENES.items():
  for gene in genes:
      row = plot_idx // ncols
      col = plot_idx % ncols

      ax = fig.add_subplot(gs[row, col])

      if gene == genes[0]:  # First gene in category
          title_color = 'darkblue' if 'FGF8' in category else 'darkred'
          title_weight = 'bold'
      else:
          title_color = 'black'
          title_weight = 'normal'

      found, variant_name, spatial_max = plot_gene(ax, mouse_analyzer, gene, show_title=False)

      if found and spatial_max > 0.1:
          genes_with_expression.append((gene, spatial_max))

      if gene == genes[0]:
          ax.set_title(f"[{category}]\n{gene}\n(max: {spatial_max:.2f})" if found else f"[{category}]\n{gene}\nNot found",
                      fontsize=8, color=title_color, weight=title_weight)
      else:
          ax.set_title(f"{gene}\n(max: {spatial_max:.2f})" if found else f"{gene}\nNot found",
                      fontsize=8)

      ax.set_xticks([])
      ax.set_yticks([])
      for spine in ax.spines.values():
          spine.set_visible(False)

      if found:
          genes_found.append((gene, variant_name))

          # Save individual image
          individual_fig, individual_ax = plt.subplots(1, 1, figsize=(5, 5))
          plot_gene(individual_ax, mouse_analyzer, gene, show_title=False)
          individual_ax.set_title(f"{gene} - {category} (Mouse)", fontsize=12, fontweight='bold')
          individual_ax.axis('off')

          individual_path = f"{individual_dir}/{gene}_{category.replace(' ', '_')}.png"
          individual_fig.savefig(individual_path, dpi=150, bbox_inches='tight', facecolor='white')
          plt.close(individual_fig)
      else:
          genes_not_found.append(gene)

      plot_idx += 1

# Hide remaining empty subplots
for idx in range(plot_idx, nrows * ncols):
  row = idx // ncols
  col = idx % ncols
  ax = fig.add_subplot(gs[row, col])
  ax.axis('off')

# figure-level suptitle removed (titleless-figure standard; caption -> LEGEND.md;
# the baked-in header overlapped the top-row panel titles on the short grids)

composite_path = f"{output_dir}/SF21A_fgf8_downstream.png"
fig.savefig(composite_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nMain figure saved to: {composite_path}")

plt.show()

# Create FGF8 correlation analysis
print("\n" + "=" * 60)
print("Creating FGF8 correlation analysis...")
print("=" * 60)

# Check if Fgf8 exists and get correlations
fgf8_variant = find_gene_variant('Fgf8', mouse_analyzer.gene_names)
if fgf8_variant:
  # Get correlations with all genes
  correlations = mouse_analyzer.get_gene_correlations(fgf8_variant, top_n=len(mouse_analyzer.gene_names))

  # Filter for our genes of interest
  gene_correlations = {}
  for gene in ALL_FGF8_GENES:
      if gene != 'Fgf8':
          variant = find_gene_variant(gene, mouse_analyzer.gene_names)
          if variant:
              for corr_gene, corr_val in correlations:
                  if corr_gene == variant:
                      gene_correlations[gene] = corr_val
                      break

  if gene_correlations:
      # Sort by correlation
      sorted_correlations = sorted(gene_correlations.items(), key=lambda x: x[1], reverse=True)

      # Create correlation plot
      fig2, ax = plt.subplots(figsize=(10, 6))

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
      ax.set_xlabel('Correlation with Fgf8', fontsize=12)
      ax.set_title('Correlation of FGF8-Related Genes with Fgf8 (Mouse)', fontsize=14, fontweight='bold')
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

      correlation_path = f"{output_dir}/SF21A_fgf8_downstream_correlations.png"
      fig2.savefig(correlation_path, dpi=300, bbox_inches='tight', facecolor='white')
      print(f"Correlation plot saved to: {correlation_path}")

      plt.show()

# Save summary
summary_path = f"{output_dir}/SF21A_fgf8_downstream_summary.txt"
with open(summary_path, 'w') as f:
  f.write("FGF8 Signaling and Related Genes - MOUSE\n")
  f.write("=" * 60 + "\n\n")

  f.write(f"Total genes analyzed: {len(ALL_FGF8_GENES)}\n")
  f.write(f"Genes found: {len(genes_found)}\n")
  f.write(f"Genes not found: {len(genes_not_found)}\n")
  f.write(f"Genes with expression >0.1: {len(genes_with_expression)}\n\n")

  for category, genes in FGF8_GENES.items():
      f.write(f"\n{category} ({len(genes)} genes):\n")
      for gene in genes:
          found = any(g[0] == gene for g in genes_found)
          if found:
              variant = find_gene_variant(gene, mouse_analyzer.gene_names)
              if variant:
                  idx = gene2idx[variant]
                  img = mouse_analyzer.images[idx]
                  max_expr = float(np.nanmax(img))
                  f.write(f"  ✓ {gene:10s} (max: {max_expr:.3f})\n")
          else:
              f.write(f"  ✗ {gene:10s}\n")

  if fgf8_variant and gene_correlations:
      f.write("\n" + "-" * 40 + "\n")
      f.write("Correlations with Fgf8:\n")
      for gene, corr in sorted_correlations:
          f.write(f"  {gene:10s}: r = {corr:.4f}\n")

  if genes_with_expression:
      f.write("\n" + "-" * 40 + "\n")
      f.write("Top Expressed Genes:\n")
      genes_with_expression.sort(key=lambda x: x[1], reverse=True)
      for i, (gene, max_val) in enumerate(genes_with_expression[:5]):
          f.write(f"  {i+1}. {gene:10s}: max = {max_val:.3f}\n")

print(f"Summary saved to: {summary_path}")

print("\n" + "=" * 60)
print("Figure generation complete!")
print(f"  Main figure: {composite_path}")
if fgf8_variant and gene_correlations:
  print(f"  Correlation plot: {correlation_path}")
print(f"  Individual images: {individual_dir}/")
print(f"  Summary: {summary_path}")
print("=" * 60)

# Print summary
print(f"\nGenes found: {len(genes_found)}/{len(ALL_FGF8_GENES)}")
if genes_with_expression:
  print(f"Genes with notable expression: {len(genes_with_expression)}")
  print("\nTop expressed genes:")
  genes_with_expression.sort(key=lambda x: x[1], reverse=True)
  for gene, max_val in genes_with_expression[:5]:
      print(f"  {gene}: max = {max_val:.3f}")

if fgf8_variant and gene_correlations:
  print("\nTop correlations with Fgf8:")
  for gene, corr in sorted_correlations[:3]:
      print(f"  {gene}: r = {corr:.3f}")

if genes_not_found:
  print(f"\nGenes not found: {', '.join(genes_not_found)}")


# %% [markdown]
# ## SF22A: BMP signaling

# %%
# Self-contained cell for generating BMP signaling genes figure for MOUSE
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

# Put the repo root on sys.path for imports.
try:
    REPO = Path(__file__).resolve().parents[2]
except NameError:  # jupytext notebook: start Jupyter from the repo root
    REPO = Path.cwd()
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from scripts.realign_mouse.build_mouse_pathway_pickle import build_analyzer
FIGURES_BASE = os.path.join(str(REPO), "figures")

from spatial_expression_analysis import (
  SpatialExpressionAnalyzer,
  SpatialAnalysisParams
)

# Define BMP signaling genes by functional category
BMP_GENES = {
  'Ligands': ['Bmp2', 'Bmp4', 'Bmp7'],
  'Receptors': ['Bmpr1a', 'Bmpr1b', 'Bmpr2', 'Acvr2a', 'Acvr2b'],
  'Antagonists/Modulators': ['Chrd', 'Nog', 'Bambi', 'Twsg1', 'Tsku', 'Grem1', 'Grem2',
                             'Bmper', 'Chrdl1', 'Chrdl2', 'Fst'],
  'Intracellular Effectors': ['Smad1', 'Smad5', 'Smad9', 'Smad4'],
  'Downstream Targets': ['Id1', 'Id2', 'Id3', 'Msx1', 'Msx2', 'Tbx5'],
  'Counter-Marker': ['Vax2']
}

# Flatten the gene list
ALL_BMP_GENES = []
for category, genes in BMP_GENES.items():
  ALL_BMP_GENES.extend(genes)

print("Loading mouse analyzer...")
print("=" * 60)

# Load or build from the deposited CR9 E13-E16 object with published parameters.
mouse_analyzer = build_analyzer()

print(f"Loaded {len(mouse_analyzer.gene_names)} genes")

# Create output directories
output_dir = os.path.join(FIGURES_BASE, "Figure_SF22")
individual_dir = f"{output_dir}/SF22A_individual_genes"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(individual_dir, exist_ok=True)

# Gene to index mapping
gene2idx = {gene: i for i, gene in enumerate(mouse_analyzer.gene_names)}

def find_gene_variant(gene_name, gene_list):
  """Find gene in list with case variations for mouse genes."""
  variations = [
      gene_name,
      gene_name.upper(),
      gene_name.lower(),
      gene_name.capitalize(),
      gene_name[0].upper() + gene_name[1:].lower() if len(gene_name) > 1 else gene_name
  ]

  if gene_name.lower().startswith('bmpr'):
      variations.append('Bmpr' + gene_name[-2:].lower())
  if gene_name.lower().startswith('acvr'):
      variations.append('Acvr' + gene_name[-2:].lower())

  for variant in variations:
      if variant in gene_list:
          return variant
  return None

def plot_gene(ax, analyzer, gene_name, show_title=True, title_size=8):
  """Plot a single gene's spatial expression."""
  found_gene = find_gene_variant(gene_name, analyzer.gene_names)

  if found_gene:
      idx = gene2idx[found_gene]
      img = analyzer.images[idx].copy()

      if analyzer.counts is not None:
          mask = (analyzer.counts < analyzer.params.mask_count_threshold)
          img[mask] = np.nan

      spatial_max = float(np.nanmax(img))

      im = ax.imshow(img,
                    origin='lower',
                    cmap='viridis',
                    aspect='equal',
                    interpolation='nearest')

      if hasattr(analyzer, 'dv_mid') and analyzer.dv_mid is not None and \
         hasattr(analyzer, 'nt_mid') and analyzer.nt_mid is not None:
          ax.axhline(y=analyzer.dv_mid, color="white", lw=0.5, alpha=0.5)
          ax.axvline(x=analyzer.nt_mid, color="white", lw=0.5, alpha=0.5)

      if show_title:
          ax.set_title(f"{gene_name}\n(max: {spatial_max:.2f})", fontsize=title_size)

      return True, found_gene, spatial_max
  else:
      ax.text(0.5, 0.5, f'{gene_name}\nNot found', ha='center', va='center',
             transform=ax.transAxes, fontsize=7, color='gray')
      if show_title:
          ax.set_title(gene_name, fontsize=title_size)
      return False, None, 0.0

print("\n" + "=" * 60)
print("Creating Mouse BMP signaling figure...")
print("=" * 60)

# Create main composite figure
total_genes = len(ALL_BMP_GENES)
ncols = 7
nrows = int(np.ceil(total_genes / ncols))

fig = plt.figure(figsize=(21, nrows * 3))

gs = gridspec.GridSpec(nrows, ncols, figure=fig,
                     hspace=0.4,
                     wspace=0.15,
                     left=0.03, right=0.97,
                     top=0.94, bottom=0.02)

# Track genes
genes_found = []
genes_not_found = []
genes_with_expression = []

# Plot genes by category
plot_idx = 0
for category, genes in BMP_GENES.items():
  for gene in genes:
      row = plot_idx // ncols
      col = plot_idx % ncols

      ax = fig.add_subplot(gs[row, col])

      if gene == genes[0]:
          title_color = 'darkred'
          title_weight = 'bold'
      else:
          title_color = 'black'
          title_weight = 'normal'

      found, variant_name, spatial_max = plot_gene(ax, mouse_analyzer, gene, show_title=False)

      if found and spatial_max > 0.1:
          genes_with_expression.append((gene, spatial_max))

      if gene == genes[0]:
          ax.set_title(f"[{category}]\n{gene}\n(max: {spatial_max:.2f})" if found else f"[{category}]\n{gene}\nNot found",
                      fontsize=7, color=title_color, weight=title_weight)
      else:
          ax.set_title(f"{gene}\n(max: {spatial_max:.2f})" if found else f"{gene}\nNot found",
                      fontsize=7)

      ax.set_xticks([])
      ax.set_yticks([])
      for spine in ax.spines.values():
          spine.set_visible(False)

      if found:
          genes_found.append((gene, variant_name))

          # Save individual image
          individual_fig, individual_ax = plt.subplots(1, 1, figsize=(5, 5))
          plot_gene(individual_ax, mouse_analyzer, gene, show_title=False)
          individual_ax.set_title(f"{gene} - {category} (Mouse)", fontsize=12, fontweight='bold')
          individual_ax.axis('off')

          individual_path = f"{individual_dir}/{gene}_{category.replace('/', '_')}.png"
          individual_fig.savefig(individual_path, dpi=150, bbox_inches='tight', facecolor='white')
          plt.close(individual_fig)
      else:
          genes_not_found.append(gene)

      plot_idx += 1

# Hide remaining empty subplots
for idx in range(plot_idx, nrows * ncols):
  row = idx // ncols
  col = idx % ncols
  ax = fig.add_subplot(gs[row, col])
  ax.axis('off')

# figure-level suptitle removed (titleless-figure standard; caption -> LEGEND.md;
# kept consistent with SF20A/SF21A which had the header-overlap issue)

composite_path = f"{output_dir}/SF22A_bmp_signaling.png"
fig.savefig(composite_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nMain figure saved to: {composite_path}")

plt.show()

# Save summary
summary_path = f"{output_dir}/SF22A_bmp_signaling_summary.txt"
with open(summary_path, 'w') as f:
  f.write("BMP Signaling Pathway Components - MOUSE\n")
  f.write("=" * 60 + "\n\n")

  f.write(f"Total genes analyzed: {len(ALL_BMP_GENES)}\n")
  f.write(f"Genes found: {len(genes_found)}\n")
  f.write(f"Genes not found: {len(genes_not_found)}\n")
  f.write(f"Genes with expression >0.1: {len(genes_with_expression)}\n\n")

  for category, genes in BMP_GENES.items():
      f.write(f"\n{category} ({len(genes)} genes):\n")
      for gene in genes:
          found = any(g[0] == gene for g in genes_found)
          if found:
              variant = find_gene_variant(gene, mouse_analyzer.gene_names)
              if variant:
                  idx = gene2idx[variant]
                  img = mouse_analyzer.images[idx]
                  max_expr = float(np.nanmax(img))
                  f.write(f"  ✓ {gene:10s} (max: {max_expr:.3f})\n")
          else:
              f.write(f"  ✗ {gene:10s}\n")

  if genes_with_expression:
      f.write("\n" + "-" * 40 + "\n")
      f.write("Top Expressed BMP Genes:\n")
      genes_with_expression.sort(key=lambda x: x[1], reverse=True)
      for i, (gene, max_val) in enumerate(genes_with_expression[:10]):
          f.write(f"  {i+1:2d}. {gene:10s}: max = {max_val:.3f}\n")

print(f"Summary saved to: {summary_path}")

print("\n" + "=" * 60)
print("Figure generation complete!")
print(f"  Main figure: {composite_path}")
print(f"  Individual images: {individual_dir}/")
print(f"  Summary: {summary_path}")
print("=" * 60)

# Print summary
print(f"\nGenes found: {len(genes_found)}/{len(ALL_BMP_GENES)}")
if genes_with_expression:
  print(f"Genes with notable expression: {len(genes_with_expression)}")
  print("\nTop 5 expressed BMP genes:")
  genes_with_expression.sort(key=lambda x: x[1], reverse=True)
  for gene, max_val in genes_with_expression[:5]:
      print(f"  {gene}: max = {max_val:.3f}")
if genes_not_found:
  print(f"\nGenes not found: {', '.join(genes_not_found)}")

