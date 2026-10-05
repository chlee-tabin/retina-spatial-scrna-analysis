#!/usr/bin/env python3
# =============================================================================
# mouse_cr9_03_assemble_h5ad.py  --  Stage 3 of the CR9 mouse atlas
#
# Assemble an export_mex() output directory (raw_counts + normalized_data MEX,
# metadata.tsv, *_embeddings.tsv, features.tsv, barcodes.tsv) into an .h5ad
# that matches the existing 20250604_mouse_RPC.h5ad convention so downstream
# figure scripts work unchanged:
#     .X            = log-normalized data            (float)
#     .raw.X        = raw counts                     (.layers left empty)
#     .obsm["X_<red>"] for each <red>_embeddings.tsv (X_pca, X_harmony, X_umap.harmony)
#     .obs          = metadata.tsv  (+ `age` alias = `stage` if age absent)
#     .var_names    = gene symbols (features.tsv)
#
# Usage:
#   python mouse_cr9_03_assemble_h5ad.py --dir <export_dir> --prefix <p> --out <file.h5ad>
# =============================================================================
import argparse, os, glob
import numpy as np, pandas as pd
import scipy.io as sio, scipy.sparse as sp
import anndata as ad


def load_mex(d, prefix):
    feats = pd.read_csv(os.path.join(d, f"{prefix}features.tsv"), header=None)[0].astype(str).tolist()
    bcs   = pd.read_csv(os.path.join(d, f"{prefix}barcodes.tsv"), header=None)[0].astype(str).tolist()
    # MEX matrices are genes x cells -> transpose to cells x genes
    counts = sio.mmread(os.path.join(d, f"{prefix}raw_counts.mtx.gz")).T.tocsr()
    data   = sio.mmread(os.path.join(d, f"{prefix}normalized_data.mtx.gz")).T.tocsr()
    assert counts.shape == (len(bcs), len(feats)), (counts.shape, len(bcs), len(feats))
    assert data.shape   == counts.shape, (data.shape, counts.shape)

    meta = pd.read_csv(os.path.join(d, f"{prefix}metadata.tsv"), sep="\t", index_col=0)
    meta.index = meta.index.astype(str)
    assert len(meta) == len(bcs), (len(meta), len(bcs))
    meta = meta.loc[bcs]                       # align to matrix cell order

    var = pd.DataFrame(index=pd.Index(feats, name=None))
    A = ad.AnnData(X=data, obs=meta, var=var)
    A.raw = ad.AnnData(X=counts, obs=meta[[]].copy(), var=var.copy())

    for f in sorted(glob.glob(os.path.join(d, f"{prefix}*_embeddings.tsv"))):
        red = os.path.basename(f)[len(prefix):-len("_embeddings.tsv")]
        emb = pd.read_csv(f, sep="\t", index_col=0)
        emb.index = emb.index.astype(str)
        A.obsm[f"X_{red}"] = emb.loc[bcs].values.astype(np.float32)

    if "age" not in A.obs.columns and "stage" in A.obs.columns:
        A.obs["age"] = A.obs["stage"].values
    return A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    A = load_mex(a.dir, a.prefix)
    print("assembled:", A.shape, "| obsm:", list(A.obsm.keys()))
    print("X  log-norm? min %.3f max %.3f integer=%s" %
          (float(A.X.min()), float(A.X.max()),
           bool(np.allclose(A.X[:200].toarray(), np.round(A.X[:200].toarray())))))
    print("raw counts integer=%s max %.1f" %
          (bool(np.allclose(A.raw.X[:200].toarray(), np.round(A.raw.X[:200].toarray()))),
           float(A.raw.X.max())))
    for c in ["DV.Score", "NT.Score", "library", "stage", "age", "annotation", "annotation_class"]:
        print("  obs has", c, "?", c in A.obs.columns)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    A.write_h5ad(a.out)
    print("wrote", a.out, "(%.0f MB)" % (os.path.getsize(a.out) / 1e6))


if __name__ == "__main__":
    main()
