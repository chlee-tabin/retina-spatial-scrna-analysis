#!/usr/bin/env python3
"""Build/cache the published CR9 mouse analyzer for Figures S20A-S22A.

Defaults read data/20260528_mouse_RPC_cr9_e13e16.h5ad (GEO GSE322831).
The parameters (5 cells per pixel / clip 0.93 / mask 5) are the ones the published
SF20A-SF22A panels were rendered with (June 2026). They differ from Figure 7's
mouse row, which was harmonized to 3 / 0.93 / 3 afterwards, so this cache is
kept separate from Figure 7's mouse_analyzer_correct.pkl. The command also writes a
companion per-gene image-maximum TSV, separate from Supplementary Table 3.

Usage: python scripts/realign_mouse/build_mouse_pathway_pickle.py [h5ad] [tag]
"""
import argparse
import sys
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from spatial_expression_analysis import SpatialAnalysisParams, load_or_build_analyzer

H5AD = REPO / "data/20260528_mouse_RPC_cr9_e13e16.h5ad"
CACHE = REPO / "data/mouse_pathway_analyzer_mc5.pkl"  # not Fig 7's cache: different params

def build_analyzer(h5ad=H5AD, cache_path=CACHE):
    params = SpatialAnalysisParams(
        bin_size=51,
        min_gene_count=30,
        min_cells_per_pixel=5,
        percentile_clip=0.93,
        smooth_sigma=1.0,
        mask_count_threshold=5,
    )
    Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
    return load_or_build_analyzer(cache_path, params, h5ad, label="mouse")

# Per-gene 2D-image maxima for the union of the three pathway gene sets.
FGF_FAMILY = ['Fgf1','Fgf2','Fgf3','Fgf4','Fgf5','Fgf6','Fgf7','Fgf8','Fgf9','Fgf10',
              'Fgf11','Fgf12','Fgf13','Fgf14','Fgf15','Fgf16','Fgf17','Fgf18',
              'Fgf20','Fgf21','Fgf22','Fgf23','Fgfr1','Fgfr2','Fgfr3','Fgfr4','Fgfrl1']
FGF8 = ['Fgf8','Dusp6','Sox2','Egr1','Snai1','Fos','En1','En2','Spry1','Spry2']
BMP  = ['Bmp2','Bmp4','Bmp7','Bmpr1a','Bmpr1b','Bmpr2','Acvr2a','Acvr2b','Chrd','Nog',
        'Bambi','Twsg1','Tsku','Grem1','Grem2','Bmper','Chrdl1','Chrdl2','Fst',
        'Smad1','Smad5','Smad9','Smad4','Id1','Id2','Id3','Msx1','Msx2','Tbx5','Vax2']
GENES = list(dict.fromkeys(FGF_FAMILY + FGF8 + BMP))

def write_gene_maxima(az, out_path):
    g2i = {g: i for i, g in enumerate(az.gene_names)}
    def variant(g):
        for v in (g, g.upper(), g.lower(), g.capitalize()):
            if v in g2i:
                return v
        return None
    with open(out_path, "w") as f:
        f.write("gene\tfound_as\tmax2d\n")
        for g in GENES:
            v = variant(g)
            if v is not None:
                mx = float(np.nanmax(az.images[g2i[v]]))
                f.write(f"{g}\t{v}\t{mx:.6f}\n")
            else:
                f.write(f"{g}\tNA\tNA\n")
    print(f"[build] wrote {out_path} ({len(GENES)} genes)", flush=True)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("h5ad", nargs="?", type=Path, default=H5AD)
    ap.add_argument("tag", nargs="?", default="cr9_e13e16")
    ap.add_argument("--cache", type=Path, default=CACHE)
    ap.add_argument("--outdir", type=Path, default=REPO / "figures/Tables")
    args = ap.parse_args()
    az = build_analyzer(args.h5ad, args.cache)
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_gene_maxima(az, args.outdir / f"mouse_pathway_maxes_{args.tag}.tsv")

if __name__ == "__main__":
    main()
