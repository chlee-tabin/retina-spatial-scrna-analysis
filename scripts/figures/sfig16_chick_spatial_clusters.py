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
# # Figure S16: Identification of Spatial Expression Clusters in Chicken Retina

# %% [markdown]
# ## Select top 20 spatial anchors

# %% tags=["cell-54"]
import sys
sys.path.append('..')
from spatial_expression_analysis import (
    SpatialAnalysisParams, SpatialExpressionAnalyzer,
    reload_control_genes, get_fixed_anchors, MarkerSelectorPy,
    cluster_by_anchor, create_cluster_summary, interrogate_cluster,
    CONTROL_GENES
)
# To customize anchors, edit `control_genes.yaml` and call reload_control_genes()
reload_control_genes()

# %% tags=["cell-55"]
chick_fixed_anchors_dict = get_fixed_anchors(chick_analyzer.gene_info, 'chick')
print("Fixed anchors from control_genes.yaml:")
for category, gene in chick_fixed_anchors_dict.items():
    expression = chick_analyzer.gene_info[chick_analyzer.gene_info['name'] == gene]['tot'].iloc[0]
    print(f"  {category}: {gene} (expression: {expression:.1f})")
# Extract gene names for marker selection
chick_fixed_anchors = list(chick_fixed_anchors_dict.values())
print(f"\nUsing {len(chick_fixed_anchors)} anchor genes for marker selection")
# Set up marker selector (replaces your MarkerSelectorPy setup)
valid_mask = chick_analyzer.gene_info['tot'] >= 2000
chick_selector = MarkerSelectorPy(
    gene_names=chick_analyzer.gene_names,
    Q=chick_analyzer.similarity_matrix,
    valid_mask=valid_mask,
    verbose=1,
    rownorm=False
)
# Select markers (replaces your select_markers call)
chick_markers = chick_selector.select_markers(K=20, fixed_anchors=chick_fixed_anchors)
print("Selected markers:")
for i, marker in enumerate(chick_markers, 1):
    expr = chick_analyzer.gene_info[chick_analyzer.gene_info['name'] == marker]['tot'].iloc[0]
    print(f"{i:2d}. {marker} ({expr:.1f})")

# %% [markdown]
# ## Build gene-to-label mapping

# %% tags=["cell-58"]
# First, create the gene-to-label mapping from control genes
chick_gene_to_label = {}
chick_control_genes = CONTROL_GENES.get('chick', {})
for category, genes in chick_control_genes.items():
    for gene in genes:
        if gene in chick_analyzer.gene_names:
            chick_gene_to_label[gene] = category

# %% [markdown]
# ## Find neighbors to anchors (network construction)

# %% tags=["cell-60"]
chick_neighbors = chick_selector.find_neighbors_to_anchors(
    m=10, # maximum candidate neighbors ("shortlist" size)
    max_rank_fraction=50/len(chick_analyzer.gene_names), # acceptable level of rank quality. Only top 50 as rank
    gene_to_label=chick_gene_to_label
)
# Keep only markers with neighbors
chick_kept_markers = [nbr.name for nbr in chick_neighbors if len(nbr.neighbors) > 0]
print(f"Kept {len(chick_kept_markers)} markers with clear neighbors")
# Cluster genes (it generates a background expression pattern so incorporates the gene expression level)
chick_Qavg = (chick_analyzer.similarity_matrix * chick_analyzer.gene_info['tot'].values.reshape(-1, 1)).sum(axis=0)
chick_Qavg = chick_Qavg / chick_Qavg.sum()
chick_cluster_idx, chick_clusters, chick_d2, chick_anchor_mapping = cluster_by_anchor(
    chick_selector, 
    fixed_anchors=chick_kept_markers, 
    null_q=chick_Qavg, 
    gene_to_label=chick_gene_to_label,
    verbose=1
)
# Add this new cell for enhanced interrogation
from spatial_expression_analysis import create_cluster_summary, interrogate_cluster
# Create summary table
chick_summary = create_cluster_summary(chick_anchor_mapping, chick_clusters, chick_markers)
print("\n=== CHICK CLUSTER SUMMARY ===")
print(chick_summary)

# %% [markdown]
# ## SF16: 20-anchor spatial pattern tiled figure

# %% tags=["cell-64"]
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

FIGURES_BASE = os.path.join(os.path.dirname(__file__), "..", "..", "figures")
# The 20 anchors in order
key_genes = chick_markers[:20]
# Get the anchors that actually formed clusters from the clustering results
cluster_forming_anchors = set(chick_kept_markers)
# Create a dictionary mapping anchor names to cluster sizes from chick_summary
anchor_to_cluster_size = {}
if 'chick_summary' in globals():
  for _, row in chick_summary.iterrows():
      if row['anchor_name'] != 'NULL':
          anchor_to_cluster_size[row['anchor_name']] = row['cluster_size']
print(f"Cluster-forming anchors from data: {len(cluster_forming_anchors)} out of {len(key_genes)}")
# Create figure with subplots
fig, axes = plt.subplots(4, 5, figsize=(12.5, 10))
axes = axes.flatten()
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
      # Check if this anchor formed a cluster and get cluster size
      forms_cluster = gene in cluster_forming_anchors
      # Create title with cluster size information
      if forms_cluster and gene in anchor_to_cluster_size:
          cluster_size = anchor_to_cluster_size[gene]
          # Format with comma for numbers >= 1000
          if cluster_size >= 1000:
              size_str = f"{cluster_size:,}"
          else:
              size_str = str(cluster_size)
          title_text = f"{gene} ({size_str})\nmax: {spatial_max:.2f}"
      else:
          # No cluster formed - just show gene name
          title_text = f"{gene}\nmax: {spatial_max:.2f}"
      ax.set_title(title_text, fontsize=11, fontweight='bold')
      # Remove axis ticks for cleaner look
      ax.set_xticks([])
      ax.set_yticks([])
      # Add red border for cluster-forming anchors
      if forms_cluster:
          # Create a rectangle patch for the border
          rect = patches.Rectangle((0, 0), 1, 1, linewidth=3,
                                  edgecolor='red', facecolor='none',
                                  transform=ax.transAxes)
          ax.add_patch(rect)
      # Add parameter info to bottom right of each subplot
      param_text = f"p={chick_analyzer.params.percentile_clip:.2f}"
      ax.text(0.98, 0.02, param_text,
              transform=ax.transAxes,
              fontsize=8, ha='right', va='bottom',
              bbox=dict(boxstyle='round,pad=0.3',
                       facecolor='white', alpha=0.8))
      # Add spatial category label if available
      if gene in chick_gene_to_label:
          category = chick_gene_to_label[gene]
          ax.text(0.02, 0.98, category,
                 transform=ax.transAxes,
                 fontsize=8, ha='left', va='top',
                 bbox=dict(boxstyle='round,pad=0.3',
                          facecolor='yellow', alpha=0.5))
      # Optionally add axis labels to first subplot only
      if i == 0:
          ax.set_xlabel("Temporal ← NT → Nasal", fontsize=10)
          ax.set_ylabel("Ventral ← DV → Dorsal", fontsize=10)
# Add overall title with legend
# fig.suptitle("Spatial Expression Patterns of Top 20 Anchor Genes\n" +
#           "Red border = forms cluster, (n) = number of genes in cluster",
#           fontsize=14, fontweight='bold')
# Adjust layout and save
plt.tight_layout()
plt.subplots_adjust(top=0.90)  # Make room for suptitle with legend
os.makedirs(os.path.join(FIGURES_BASE, "Figure_SF16"), exist_ok=True)
plt.savefig(os.path.join(FIGURES_BASE, "Figure_SF16", "SF16_chick_spatial_clusters.png"), dpi=300, bbox_inches='tight')
plt.show()
# Print detailed summary
print("\n" + "="*60)
print("CLUSTER FORMATION SUMMARY")
print("="*60)
print(f"Total anchors analyzed: {len(key_genes)}")
print(f"Anchors that formed clusters: {sum(1 for g in key_genes if g in cluster_forming_anchors)}")
print(f"Anchors that didn't form clusters: {sum(1 for g in key_genes if g not in cluster_forming_anchors)}")
print("\nCluster sizes for the 20 anchors:")
for gene in key_genes:
  if gene in anchor_to_cluster_size:
      size = anchor_to_cluster_size[gene]
      size_str = f"{size:,}" if size >= 1000 else str(size)
      category = chick_gene_to_label.get(gene, "Unknown")
      print(f"  {gene:15} → {size_str:>6} genes ({category})")
  else:
      print(f"  {gene:15} → No cluster formed")
# Summary statistics
if anchor_to_cluster_size:
  sizes = [anchor_to_cluster_size[g] for g in key_genes if g in anchor_to_cluster_size]
  print(f"\nCluster size statistics for these 20 anchors:")
  print(f"  Mean: {np.mean(sizes):.1f} genes")
  print(f"  Median: {np.median(sizes):.1f} genes")
  print(f"  Range: {min(sizes)} - {max(sizes):,} genes")

# %% [markdown]
# ## Table S1: Anchor gene pairs and cluster members

# %% [markdown]
# ### Supplemental Table 1

# %% tags=["cell-66"]
import pandas as pd
# Create detailed list - ONLY genes actually in clusters
detailed_data = []
for _, row in chick_summary.iterrows():
  if row['anchor_name'] != 'NULL':
      anchor = row['anchor_name']
      category = row['category']
      cluster_size = row['cluster_size']
      if anchor in chick_analyzer.gene_names:
          # Get ALL cluster members
          cluster_genes = chick_clusters.get(anchor, [])
          # Get enough correlations to cover all cluster members
          correlations = chick_analyzer.get_gene_correlations(anchor,
                                                             top_n=len(chick_analyzer.gene_names))
          corr_dict = dict(correlations)
          # Add each cluster member with its correlation
          for gene in cluster_genes:
              if gene in corr_dict:
                  # Find rank in correlation list
                  rank = next((i for i, (g, _) in enumerate(correlations, 1) if g == gene), None)
                  detailed_data.append({
                      'Anchor': anchor,
                      'Category': category,
                      'Cluster_Size': cluster_size,
                      'Correlation_Rank': rank,
                      'Gene': gene,
                      'Correlation': corr_dict[gene]
                  })
# Create DataFrame and sort
detailed_table = pd.DataFrame(detailed_data)
detailed_table = detailed_table.sort_values(['Anchor', 'Correlation_Rank'])
# Save to TSV
os.makedirs(os.path.join(FIGURES_BASE, "Tables"), exist_ok=True)
detailed_table.to_csv(os.path.join(FIGURES_BASE, "Tables", "TableS1_chick_cluster_members.tsv"), sep='\t', index=False)
print(f"Saved {len(detailed_table)} cluster member entries")
print(f"This should match total of all cluster sizes: {chick_summary[chick_summary['anchor_name'] != 'NULL']['cluster_size'].sum()}")
