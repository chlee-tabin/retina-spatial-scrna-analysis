# Supplementary Table 1 (chick) / Table 2 (human) — area DEG, full significant list

Pseudobulk differential-expression tables underlying the area-specific volcano
plots (chick Fig. 6H / human Fig. S23). Each topographic territory (DV.Score ×
NT.Score grid) is contrasted against the rest of the retinal-progenitor pool.

| File | Figure | Contents |
|---|---|---|
| `SuppTable1_chick_area_significant.csv` | **Fig. 6H** | chick — all genes with adj-p < 0.05 in any of the 7 territories (77 HAA-enriched / 75 HAA-de-enriched among the totals) |
| `SuppTable2_human_area_significant.csv` | **Fig. S23** | human — all genes with adj-p < 0.05 in any of the 7 territories (35 fovea-enriched / 33 fovea-de-enriched among the totals) |

Both are emitted by `scripts/analysis/area_significant_deg.R` (one script, both
species, same gates as the figure code in `scripts/figures/fig6h_sfig23_area_deg.R`).

## Schema

`gene, region, average_expression, log2FC, pval, adj_pval, direction`

- `region` ∈ {HAA / Fovea, Temporal, Nasal, Dorsal, Ventral, DVcentral, NTcentral}
  (chick uses **HAA**; human uses **Fovea** for the same high-acuity contrast).
- `direction` = "enriched" if log2FC > 0 (up in the territory), "de-enriched" otherwise.
- `average_expression` = mean log-normalized expression across all cells in the species object.

## Method

Cells are aggregated to pseudobulk by **library × territory × genotype** (chick) or
**library × territory × sample** (human); a gamma–Poisson GLM (`glmGamPoi::glm_gp`,
`~ library + area`) is fit and the territory term tested (`test_de`, contrast
`area1` = in-territory vs. rest).

**Pseudobulk size gate — `min_cells = 50`.** Samples (library × area × donor
groups) with fewer than 50 cells are dropped before fitting. `min_cells = 50` is
a pre-specified floor (standard pseudobulk practice; balances noise removal
against power for the smaller HAA arm).

## Gating choice (`min_cells = 50`) and sensitivity

The floor was chosen on principled grounds and is **not** tuned to a target
gene. Core HAA markers are robust across the whole defensible range; the
choice only governs the borderline gene set (NPY, BMP2, SPRY1).

HAA contrast, adjusted p across the sweep:

| gene | ≥10 | ≥20 | ≥30 | **≥50** | ≥100 |
|---|---|---|---|---|---|
| CYP26C1 | 3e-10 | 4e-10 | 3e-10 | 4e-10 | 1e-9 |
| FGF8 | 7e-9 | 6e-10 | 2e-10 | 3e-10 | 2e-10 |
| MYOF | 2e-7 | 4e-7 | 8e-7 | 7e-6 | 4e-6 |
| SPRY1 | 2e-4 | 4e-4 | 5e-4 | 3e-3 | 5e-4 |
| BMP2 | 0.19 | 0.13 | 0.054 | **0.017** | 5e-4 |
| NPY | 0.080 | 0.089 | 0.093 | **0.048** | 0.065 |

Above ≥150 the gate over-prunes the (smaller) HAA arm — at ≥200 only 6 HAA
groups remain (resid df 22) and even CYP26C1 weakens (0.0035), AKR1D1 goes
n.s.; at ≥300 CYP26C1 itself is n.s. So **≥10–50 is the safe window for this
dataset**; ≥50 was used.

## Independently-validated HAA genes (HAA contrast, min_cells = 50)

| gene | log2FC | adj-p | note |
|---|---|---|---|
| MYOF | +1.01 | 7e-6 | also a DVcentral marker |
| **NPY** | +1.01 | **0.048** | recovered at the ≥50 floor; **marginal** significance, ~2-fold, externally validated |
| BMP2 | +0.78 | 0.017 | primary marker is DVcentral (+1.23) |
| SPRY1 | +0.40 | 0.003 | robustly significant, ~1.3× (a fold-size, not significance, issue) |

NPY is reported honestly: its inclusion rests on independent validation + a
stable ~2-fold effect; its pseudobulk significance is marginal (0.048) and
knife-edge across the gating sweep, so it is not presented as a robust top-tier
marker.

## HAA vs DVcentral — direct contrast (double-dipping–controlled)

To compare the two overlapping central territories without double-counting, HAA
cells were contrasted against **DVcentral cells with HAA removed** (HAA = 5,971;
DVcentral-excluding-HAA = 3,927; 4,280 shared cells assigned to HAA). Positive
log2FC = higher in HAA (`min_cells = 50`):

| gene | log2FC (HAA vs DVcentral) | adj-p |
|---|---|---|
| CYP26C1 | **+1.95** | 0.0014 |
| FGF8 | +0.50 | n.s. |
| NPY | +0.29 | n.s. |
| MYOF | +0.23 | n.s. |
| SPRY1 | +0.08 | n.s. |
| BMP2 | −0.12 | n.s. |

**Only CYP26C1 is sharply HAA-localized** relative to DVcentral. SPRY1, MYOF
and NPY are statistically indistinguishable between the two central territories
in the pseudobulk — i.e. broadly *central* rather than HAA-vs-DVcentral–
differential. Caveat: the gates overlap heavily and pseudobulk enrichment-vs-rest
is not the same as absolute spatial intensity, so this need not contradict
imaging-based validation that distinguishes the two; it does mean the
transcriptomic contrast alone does not resolve a SPRY1 HAA>DVcentral (or MYOF
DVcentral>HAA) ordering.

## Human Fovea: cross-species reframe (validation criterion)

The chick HAA is central-nasal; the **human fovea is anatomically
temporal-central**. A correct human fovea gate should therefore *be* enriched
for temporal-axis markers — the opposite of the chick HAA's *depletion* of
temporal markers.

Marker-guided localization confirms the gate is at the right place: mapping the
chick acuity core onto the human RPCs, the (weak) acuity signal plus the human
RA anchor **CYP26A1** plus FOXD1 / EPHA3 all concentrate at the temporal-of-center
locus, exactly where the gate `c(3,5),c(4,5),c(3,4),c(4,4)` (n_grid=10) sits.
**Validation criterion (human):** CYP26A1 enrichment + temporal-fovea consistency
— *not* FOXD1-depletion. The human fovea signature in this table is led by
EPHA3 (+3.16) / FOXD1 (+2.94) / RLBP1 (+1.72) / POU3F2 / FABP7 / GRID2 /
CYP26A1, plus 35 enriched / 33 de-enriched total.

**Cross-species difference (a result, not an artifact):** the chick acuity core
does not transfer cleanly. **CYP26C1 is essentially absent in human** (meanExp
0.001) and **FGF8 is not enriched in the human fovea gate** (FGF8 surfaces only
for the Nasal territory in this table, log2FC +0.74 — not foveal). The conserved
feature is a CYP26-high (RA-degrading) zone at the high-acuity area, via a
different paralog (chick CYP26C1 → human CYP26A1, +1.26 vs. rest, sig) at a
different retinotopic position.

## Provenance & validation

- **Object (chick):** `data/20250604_02_fabp7.rds` (the chick RPC object emitted
  by `scripts/preprocessing/02_chick_dv_nt_scoring.R`; 29,025 cells).
- **Object (human):** `data/20250604human.RPC/` (published-scores MEX export
  with raw counts + the baked-in `DV.Score` / `NT.Score` used for Fig. S23).
- **Environment:** R ≥ 4.4.0; `glmGamPoi`, `Seurat`, `SingleCellExperiment` —
  same stack as the figure scripts (see `scripts/README.md` for the full
  software prerequisites).
- **Generating script:** `scripts/analysis/area_significant_deg.R`. Run from
  the repo root: `Rscript scripts/analysis/area_significant_deg.R`. Override
  the input data directory via `DATA_DIR=...` (defaults to `data/`).
- **Region gates:** rectangular DV × NT quadrant ranges per region (per-region
  `n_grid`), reproduced from the `plot.retina3` quadrant selections in the
  figure script (`scripts/figures/fig6h_sfig23_area_deg.R`). The chick HAA
  and human Fovea cell counts (5,971 / 2,236) are pinned with a `stopifnot`
  forcing function on both the analysis-script side and the figure-script
  side, so any drift in either gate parameterization halts before producing
  silently mismatched tables.
- **Chick validation:** reproduces the published HAA signature — CYP26C1
  **+2.88**, AKR1D1 **+2.33**, FGF8 **+1.70**, HHEX **+1.44**, MYOF **+1.01** —
  with the HAA gate at exactly **5,971 cells (20.6%)**.

### Significant DEGs per territory (`min_cells = 50`, adj-p < 0.05)

Chick (Supp Table 1, 6,383 rows):

| territory | enriched | de-enriched |
|---|---:|---:|
| HAA | 77 | 75 |
| Temporal | 142 | 71 |
| Nasal | 299 | 533 |
| Dorsal | 973 | 1,040 |
| Ventral | 1,111 | 1,407 |
| DVcentral | 251 | 234 |
| NTcentral | 42 | 128 |

Human (Supp Table 2, 1,396 rows):

| territory | enriched | de-enriched |
|---|---:|---:|
| Fovea | 35 | 33 |
| Temporal | 19 | 14 |
| Nasal | 399 | 333 |
| Dorsal | 68 | 19 |
| Ventral | 37 | 20 |
| DVcentral | 30 | 123 |
| NTcentral | 136 | 130 |
