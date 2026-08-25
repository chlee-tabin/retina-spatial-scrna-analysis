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
# # Figure 7A-G: Cross-Species 2D Topographic Maps

# %% [markdown]
# ## F7A-G: Cross-species orthologous gene expression (3 species x 7 groups)

# %% tags=["cell-181"]
# Self-contained cell for generating orthologous gene expression figure
# This cell loads all necessary data and creates a tiled figure without requiring pre-run cells
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
# Put the repo root on sys.path for the module import
try:
    REPO = Path(__file__).resolve().parents[2]  # repo root; keeps the script runnable from any CWD
except NameError:  # running as a notebook kernel (jupytext): no __file__ — start Jupyter from the repo root
    REPO = Path.cwd()
    print(f"NOTE: no __file__ (notebook mode); assuming repo root = {REPO}")
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from spatial_expression_analysis import (
  SpatialExpressionAnalyzer,
  SpatialAnalysisParams,
  load_or_build_analyzer
)
# Define the genes list - 24 genes total, will be arranged in 4 rows of 6
GENE_LIST = [
  'Aldh1a1', 'Aldh1a3', 'Cyp26a1', 'Cyp26c1', 'Fgf8', 'Tbx2',
  'Tbx3', 'Tbx5', 'Efnb1', 'Efnb2', 'Ephb2', 'Vax1',  # Note: EphrinB1/B2 -> Efnb1/Efnb2, Vax -> Vax1
  'Chrdl1', 'Foxd1', 'Foxg1', 'Soho1', 'Hmx1', 'Cyp1b1',  # Ventroptin -> Chrdl1, Soho-1 -> Soho1
  'Bmp2', 'Bambi', 'Bmpr1b', 'Hmx1', 'Nr2f1', 'Nr2f2'  # HMX1 appears twice (rows 17 and 22)
]
# Species order for display
SPECIES_ORDER = ['chick', 'human', 'mouse']
SPECIES_TITLES = {'chick': 'Chick', 'human': 'Human', 'mouse': 'Mouse'}
print("Loading analyzers for each species...")
print("=" * 60)
# Load or create analyzers for each species
analyzers = {}
# Load Chick analyzer with CORRECT parameters and ERROR HANDLING
print("Loading chick analyzer...")
chick_pickle = f"{REPO}/data/chick_analyzer_correct.pkl"  # one cache per repo, not per CWD
chick_params = SpatialAnalysisParams(
    bin_size=51,
    min_gene_count=20,  # Aligned with other species
    min_cells_per_pixel=3,
    percentile_clip=0.93,
    smooth_sigma=1.0,
    mask_count_threshold=3
)
# The deposited h5ad from GEO GSE322831 is the supported input (see README);
# the cache self-invalidates when the source h5ad or parameters change.
chick_analyzer = load_or_build_analyzer(chick_pickle, chick_params,
                                       f"{REPO}/data/20250604_chick_RPC.h5ad", label="chick")
analyzers['chick'] = chick_analyzer
print(f"  Loaded {len(chick_analyzer.gene_names)} genes")
# Load Human analyzer with CORRECT parameters and ERROR HANDLING
print("Loading human analyzer...")
human_pickle = f"{REPO}/data/human_analyzer_correct.pkl"  # one cache per repo, not per CWD
human_params = SpatialAnalysisParams(
    bin_size=40,  # Changed from 20 to 40 as per notebook
    min_gene_count=15,  # This should keep CYP26C1
    min_cells_per_pixel=3,
    percentile_clip=0.93,
    smooth_sigma=1.0,
    mask_count_threshold=3
)
# The deposited h5ad from GEO GSE322831 is the supported input (see README);
# the cache self-invalidates when the source h5ad or parameters change.
human_analyzer = load_or_build_analyzer(human_pickle, human_params,
                                       f"{REPO}/data/20250604_human_RPC.h5ad", label="human")
analyzers['human'] = human_analyzer
print(f"  Loaded {len(human_analyzer.gene_names)} genes")
# Load Mouse analyzer with CORRECT parameters and ERROR HANDLING
print("Loading mouse analyzer...")
mouse_pickle = f"{REPO}/data/mouse_analyzer_correct.pkl"  # one cache per repo, not per CWD
mouse_params = SpatialAnalysisParams(
    bin_size=51,
    min_gene_count=30,  # matches human/chick floors scaled per species
    # Harmonized with the chick/human rows of this same figure (and with SF19's
    # mouse analyzer): 3 / 0.93 / 3. The previous 5 / 0.95 / 5 came from "high
    # density defaults" — a >50k-cell bucket no mouse object ever occupied
    # (old arm: 25,202 cells; CR9 e13e16: 26,505 — both medium density,
    # which derives exactly these values).
    min_cells_per_pixel=3,
    percentile_clip=0.93,
    smooth_sigma=1.0,
    mask_count_threshold=3
)
# The deposited h5ad from GEO GSE322831 is the supported input (see README);
# the cache self-invalidates when the source h5ad or parameters change.
mouse_analyzer = load_or_build_analyzer(mouse_pickle, mouse_params,
                                       f"{REPO}/data/20260528_mouse_RPC_cr9_e13e16.h5ad", label="mouse")
analyzers['mouse'] = mouse_analyzer
print(f"  Loaded {len(mouse_analyzer.gene_names)} genes")
def get_gene_image(analyzer, gene_name, species):
  """Get spatial expression image for a gene, handling case variations."""
  # Create gene to index mapping
  gene2idx = {gene: i for i, gene in enumerate(analyzer.gene_names)}
  # Try different case variations
  gene_variations = [
      gene_name,
      gene_name.upper(),
      gene_name.lower(),
      gene_name.capitalize()
  ]
  # Special handling for some gene name variations
  if gene_name.lower() == 'efnb1':
      gene_variations.extend(['Ephrinb1', 'EphrinB1', 'EFNB1'])
  elif gene_name.lower() == 'efnb2':
      gene_variations.extend(['Ephrinb2', 'EphrinB2', 'EFNB2'])
  elif gene_name.lower() == 'vax1':
      # Case/orthography variants of VAX1 only. 'Vax2' was previously in this
      # list and is a DIFFERENT gene: a VAX1 miss would silently render VAX2
      # under a panel titled 'Vax1'.
      gene_variations.extend(['Vax', 'VAX', 'VAX1'])
  elif gene_name.lower() == 'chrdl1':
      gene_variations.extend(['Ventroptin', 'VENTROPTIN'])
  elif gene_name.lower() == 'soho1':
      gene_variations.extend(['Soho-1', 'SOHO-1', 'Soho'])
  for variant in gene_variations:
      if variant in analyzer.gene_names:
          idx = gene2idx[variant]
          img = analyzer.images[idx].copy()
          # Apply mask if available
          if analyzer.counts is not None:
              mask = (analyzer.counts < analyzer.params.mask_count_threshold)
              img[mask] = np.nan
          return img, variant
  return None, None
# Create output directory
FIGURES_BASE = os.path.join(str(REPO), "figures")  # REPO is notebook-safe; a bare __file__ here is not
output_dir = os.path.join(FIGURES_BASE, "Figure7")
os.makedirs(output_dir, exist_ok=True)
os.makedirs(os.path.join(output_dir, "variants"), exist_ok=True)
print("\n" + "=" * 60)
print("Creating tiled figure...")
print("=" * 60)
# Setup figure - 6 columns x 4 rows of gene groups
# Adjust aspect ratio to be less vertically compressed
fig = plt.figure(figsize=(15, 20))  # Changed aspect ratio - taller and narrower
# Create grid with MINIMAL spacing for tight layout
# Adjusted top margin since we're removing the title
gs = gridspec.GridSpec(12, 6, figure=fig,
                     hspace=0.05,  # Very minimal vertical spacing
                     wspace=0.02,  # Very minimal horizontal spacing
                     left=0.02, right=0.98,  # Minimal margins
                     top=0.99, bottom=0.01)  # More space at top since no title
# Track which genes were found
genes_found = {species: [] for species in SPECIES_ORDER}
genes_not_found = {species: [] for species in SPECIES_ORDER}
# Process each gene
for gene_idx, gene in enumerate(GENE_LIST):
  # Calculate position in grid
  col = gene_idx % 6
  row_group = gene_idx // 6
  print(f"Processing gene {gene_idx+1}/{len(GENE_LIST)}: {gene}")
  # Process each species
  for species_idx, species in enumerate(SPECIES_ORDER):
      # Calculate subplot position
      # Each gene group spans 3 rows (one for each species)
      row = row_group * 3 + species_idx
      ax = fig.add_subplot(gs[row, col])
      # Get gene image
      analyzer = analyzers[species]
      img, found_gene = get_gene_image(analyzer, gene, species)
      if img is not None:
          # Plot the image with CORRECT ORIGIN
          im = ax.imshow(img,
                        origin='lower',  # IMPORTANT: Match the original notebook
                        cmap='viridis',
                        aspect='equal',
                        interpolation='nearest')
          # Add max value annotation with adaptive formatting
          max_val = np.nanmax(img)
          # Use adaptive formatting: avoids "0.0" for small values
          if max_val >= 1:
              max_str = f'{max_val:.2f}'
          elif max_val >= 0.01:
              max_str = f'{max_val:.3f}'
          else:
              max_str = f'{max_val:.2e}'  # Scientific notation for very small values
          ax.text(0.95, 0.05, max_str,
                 transform=ax.transAxes,
                 fontsize=6, color='white',
                 ha='right', va='bottom',
                 bbox=dict(boxstyle='round,pad=0.3',
                          facecolor='black', alpha=0.5))
          # Add midlines if available
          if hasattr(analyzer, 'dv_mid') and analyzer.dv_mid is not None and hasattr(analyzer, 'nt_mid') and analyzer.nt_mid is not None:
              ax.axhline(y=analyzer.dv_mid, color="white", lw=0.5, alpha=0.5)
              ax.axvline(x=analyzer.nt_mid, color="white", lw=0.5, alpha=0.5)
          # Format title - CONSISTENT formatting (not bold for any species)
          if species_idx == 0:  # Only add gene name for top species
              ax.set_title(f"{gene}\n{SPECIES_TITLES[species]}", fontsize=8, pad=1)
          else:
              ax.set_title(f"{SPECIES_TITLES[species]}", fontsize=7, pad=1)
          genes_found[species].append(found_gene)
          # Save individual image with CORRECT ORIGIN
          individual_fig, individual_ax = plt.subplots(1, 1, figsize=(5, 5))
          individual_ax.imshow(img,
                             origin='lower',  # IMPORTANT: Match the original notebook
                             cmap='viridis',
                             aspect='equal',
                             interpolation='nearest')
          # Add max value annotation to individual images too
          individual_ax.text(0.95, 0.05, f'Max: {max_str}',
                           transform=individual_ax.transAxes,
                           fontsize=10, color='white',
                           ha='right', va='bottom',
                           bbox=dict(boxstyle='round,pad=0.5',
                                    facecolor='black', alpha=0.6))
          if hasattr(analyzer, 'dv_mid') and analyzer.dv_mid is not None and hasattr(analyzer, 'nt_mid') and analyzer.nt_mid is not None:
              individual_ax.axhline(y=analyzer.dv_mid, color="white", lw=0.5, alpha=0.5)
              individual_ax.axvline(x=analyzer.nt_mid, color="white", lw=0.5, alpha=0.5)
          individual_ax.set_title(f"{gene} - {SPECIES_TITLES[species]}", fontsize=12, fontweight='bold')
          individual_ax.axis('off')
          # Save individual images as both PDF and PNG
          individual_path_pdf = f"{output_dir}/variants/{gene}_{species}.pdf"
          individual_fig.savefig(individual_path_pdf, bbox_inches='tight')
          individual_path_png = f"{output_dir}/variants/{gene}_{species}.png"
          individual_fig.savefig(individual_path_png, dpi=150, bbox_inches='tight')
          plt.close(individual_fig)
      else:
          # Gene not found - show empty plot with message
          ax.text(0.5, 0.5, 'Not found', ha='center', va='center',
                 transform=ax.transAxes, fontsize=6, color='gray')
          if species_idx == 0:
              ax.set_title(f"{gene}\n{SPECIES_TITLES[species]}", fontsize=8, pad=1)
          else:
              ax.set_title(f"{SPECIES_TITLES[species]}", fontsize=7, pad=1)
          genes_not_found[species].append(gene)
      # Remove axes
      ax.set_xticks([])
      ax.set_yticks([])
      ax.spines['top'].set_visible(False)
      ax.spines['right'].set_visible(False)
      ax.spines['bottom'].set_visible(False)
      ax.spines['left'].set_visible(False)
# Save the composite figure as PDF (vector format)
composite_path_pdf = f"{output_dir}/F7A-G_cross_species.pdf"
fig.savefig(composite_path_pdf, bbox_inches='tight', facecolor='white')
print(f"\nComposite figure saved to: {composite_path_pdf}")
# Also save high-res PNG as backup
composite_path_png = f"{output_dir}/F7A-G_cross_species.png"
fig.savefig(composite_path_png, dpi=300, bbox_inches='tight', facecolor='white')
print(f"PNG version saved to: {composite_path_png}")
plt.show()
# Print summary
print("\n" + "=" * 60)
print("Summary:")
print("=" * 60)
for species in SPECIES_ORDER:
  print(f"\n{SPECIES_TITLES[species]}:")
  print(f"  Genes found: {len(genes_found[species])}/{len(GENE_LIST)}")
  if genes_not_found[species]:
      print(f"  Genes not found: {', '.join(genes_not_found[species])}")
print("\n" + "=" * 60)
print("Figure generation complete!")
print(f"Main figure (PDF): {composite_path_pdf}")
print(f"Main figure (PNG): {composite_path_png}")
print(f"Individual images: {output_dir}/variants/")
print("=" * 60)
