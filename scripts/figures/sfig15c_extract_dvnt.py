#!/usr/bin/env python3
"""
Extract DV.Score / NT.Score + selected gene expression from the chick RPC h5ad
for the 'flower / orange-peel' reprojection exploration.

Inputs: data/20250604_chick_RPC.h5ad (GEO GSE322831).
Standalone outputs: figures/Figure_SF15/dvnt_genes.tsv and extract_summary.txt.
"""
import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

REPO = Path(__file__).resolve().parents[2]
H5 = REPO / "data/20250604_chick_RPC.h5ad"
OUTDIR = REPO / "figures/Figure_SF15"

TARGETS   = ["FGF8", "CYP26C1", "BMP2"]
LANDMARKS = ["TBX5", "TBX3", "TBX2", "ALDH1A1", "EFNB1", "EFNB2",   # dorsal
             "VAX1", "CHRDL1", "ALDH1A3", "CYP26A1",               # ventral
             "FOXG1", "HMX1", "SOHO1", "EFNA5", "EFNA2",           # nasal
             "FOXD1", "EPHA3",                                      # temporal
             "CYP1B1", "BAMBI", "GDF6", "NR2F2", "CRABP1", "ID2",   # extras
             "ID3", "NTN1", "GDF10", "SIX3", "SIX6", "NTN4"]

_LINES = []
def log(s):
    s = str(s)
    print(s, flush=True)
    _LINES.append(s)

def get_matrix(A, genes):
    """Return (dense ndarray [n_obs x len(genes)], source_label). Explicit choice."""
    sub = A[:, genes]
    src = None
    M = None
    for cand in ["data", "lognorm", "logcounts", "lognormalized", "normalized"]:
        if cand in A.layers:
            M, src = sub.layers[cand], f"layer:{cand}"
            break
    if M is None:
        M, src = sub.X, "X"
    Md = M.toarray() if sp.issparse(M) else np.asarray(M)
    Md = np.asarray(Md, dtype=float)
    if src == "X":
        mx = float(np.nanmax(Md)) if Md.size else 0.0
        is_int = bool(np.allclose(Md, np.round(Md)))
        log(f"matrix: X[:,genes] max={mx:.3f} integer_like={is_int}")
        if is_int:
            ncol = None
            for c in ["nCount_RNA", "nCount_originalexp", "total_counts", "nCount_SCT"]:
                if c in A.obs.columns:
                    ncol = c
                    break
            if ncol is not None:
                tot = pd.to_numeric(A.obs[ncol], errors="coerce").to_numpy().astype(float)
                tot[tot <= 0] = np.nan
                Md = np.log1p(Md / tot[:, None] * 1e4)
                src = f"X_counts->log1p_CP10k[{ncol}]"
                log(f"matrix: X looked like counts; normalized to log1p(CP10k) via obs['{ncol}']")
            else:
                Md = np.log1p(Md)
                src = "X_counts->log1p(raw)"
                log("matrix: X looked like counts and no depth column; used log1p(raw) [NOT depth-normalized]")
        else:
            log("matrix: X assumed already normalized; used as-is")
    else:
        log(f"matrix: using {src} (assumed normalized log)")
    return Md, src

def extract_dataframe(h5ad=H5):
    """Read the published scores and expression, preserving the source extraction."""
    _LINES.clear()
    log(f"h5ad: {h5ad}")
    A = ad.read_h5ad(h5ad)
    log(f"shape: n_obs={A.n_obs}  n_vars={A.n_vars}")
    log(f"layers: {list(A.layers.keys())}   raw_present={A.raw is not None}")

    score_cols = [c for c in A.obs.columns if "core" in c]  # *.Score / *.score
    log(f"obs Score-like columns: {score_cols}")
    for req in ["DV.Score", "NT.Score"]:
        if req not in A.obs.columns:
            log(f"FATAL: required obs column '{req}' not found.")
            log(f"all obs columns: {list(A.obs.columns)}")
            sys.exit(2)

    annot_cols = [c for c in A.obs.columns
                  if c.lower() in ("annotation", "celltype", "cell_type",
                                   "seurat_clusters", "clusters", "library", "technology")]
    log(f"annotation/grouping columns kept: {annot_cols}")

    vset = set(map(str, A.var_names))
    present_t = [g for g in TARGETS if g in vset]
    missing_t = [g for g in TARGETS if g not in vset]
    present_l = [g for g in LANDMARKS if g in vset]
    missing_l = [g for g in LANDMARKS if g not in vset]
    log(f"TARGETS present={present_t}  MISSING={missing_t}")
    log(f"LANDMARKS present={present_l}  MISSING={missing_l}")
    if missing_t:
        # try case-insensitive recovery for reporting
        upper = {v.upper(): v for v in map(str, A.var_names)}
        for g in missing_t:
            if g.upper() in upper:
                log(f"  note: '{g}' absent but case-variant present as '{upper[g.upper()]}'")
    if not present_t:
        log("FATAL: none of the target genes present; aborting.")
        sys.exit(3)

    genes = present_t + present_l
    Md, msrc = get_matrix(A, genes)
    expr = pd.DataFrame(Md, columns=genes, index=A.obs_names.astype(str))

    keep_obs = ["DV.Score", "NT.Score"]
    for c in ["Dorsal.Score1", "Ventral.Score1", "Nasal.Score1", "Temporal.Score1"]:
        if c in A.obs.columns:
            keep_obs.append(c)
    keep_obs += annot_cols
    obs = A.obs[keep_obs].copy()
    obs.index = obs.index.astype(str)

    n_before = A.n_obs
    df = obs.join(expr, how="inner")
    log(f"row-count check: obs n={n_before}  joined n={df.shape[0]}  preserved={df.shape[0] == n_before}")
    if df.shape[0] != n_before:
        log("WARNING: join changed row count — index mismatch between obs and expr.")

    for c in ["DV.Score", "NT.Score"]:
        v = pd.to_numeric(df[c], errors="coerce").to_numpy()
        qs = np.nanpercentile(v, [1, 5, 25, 50, 75, 95, 99])
        log(f"{c}: n_nan={int(np.isnan(v).sum())} min={np.nanmin(v):.3f} max={np.nanmax(v):.3f} "
            f"mean={np.nanmean(v):.3f} sd={np.nanstd(v):.3f} "
            f"q[1,5,25,50,75,95,99]={np.round(qs,3).tolist()}")

    dv = pd.to_numeric(df["DV.Score"], errors="coerce").to_numpy()
    nt = pd.to_numeric(df["NT.Score"], errors="coerce").to_numpy()
    dvn = (dv - np.nanmean(dv)) / np.nanstd(dv)
    ntn = (nt - np.nanmean(nt)) / np.nanstd(nt)
    both_hi = float(np.nanmean((np.abs(dvn) > 1.0) & (np.abs(ntn) > 1.0)))
    one_hi  = float(np.nanmean((np.abs(dvn) > 1.0) & (np.abs(ntn) <= 0.3)))
    pear    = float(np.corrcoef(dv[~np.isnan(dv) & ~np.isnan(nt)], nt[~np.isnan(dv) & ~np.isnan(nt)])[0, 1])
    log(f"joint support: corner_frac(both|z|>1)={both_hi:.3f}  axis_edge_frac(one|z|>1,other<=0.3)={one_hi:.3f}  pearson(DV,NT)={pear:.3f}")
    log(f"interpretation hint: corner_frac << axis_edge_frac => diamond/cross support (curvature-deficit corners) ; expect periphery compression.")

    log(f"matrix_source={msrc}")
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", type=Path, default=H5)
    ap.add_argument("--outdir", type=Path, default=OUTDIR)
    args = ap.parse_args()
    df = extract_dataframe(args.h5ad)
    args.outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = args.outdir / "dvnt_genes.tsv"
    df.to_csv(out_tsv, sep="\t")
    log(f"WROTE {out_tsv}  ({df.shape[0]} rows x {df.shape[1]} cols)")
    (args.outdir / "extract_summary.txt").write_text("\n".join(_LINES) + "\n")


if __name__ == "__main__":
    main()
