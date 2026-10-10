#!/usr/bin/env python3
"""
Figure S15C: curvature-aware reprojection of the chick DV/NT topographic maps.

Run without flags using data/20250604_chick_RPC.h5ad from GEO GSE322831.
Defaults reproduce the published, titleless 12-gene grid with locked geometry.

Exploratory, mostly-aesthetic transform that:
  (1) keeps the pole-proximal HAA center (origin) ~undistorted,
  (2) expands the nasal/temporal periphery (de-saturates module-score compression),
  (3) opens orange-peel relief "rips" (gores) that widen toward the rim,
to make the flat (DV, NT) plane resemble a flat-mounted whole-mount.

Model (single, internally consistent option):
  - center at the HAA pole; scale each axis robustly; go to polar (r, theta).
  - treat score-radius as orthographic foreshortening of a spherical cap & invert:
        rho = arcsin( r * sin(rho_max) )          [periphery expands as r->1]
  - lay out azimuthal-equidistant: display radius R = rho (colatitude).
  - cut K radial slits; within each gore open a wedge gap carrying the curvature
        deficit  Delta(rho) = 2*pi*(1 - sin(rho)/rho)  (0 at pole, grows to rim).
  - anisotropy: rho_max larger toward NT so nasal/temporal petals run longer.

Renders BINNED + SMOOTHED maps (matching the manuscript's spatial pipeline):
  square panel = pcolormesh heatmap on the (NT,DV) grid;
  flower panel = the SAME grid warped, drawn as filled quads, with quads that
                 cross a rip (beyond a small joined center) dropped -> real slits.

Orientation matches the FISH flat-mount convention: nasal=right, dorsal=top,
temporal=left, ventral=bottom; HAA at center. Every knob is explicit.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
PUBLISHED_GENES = "ALDH1A1,TBX5,TBX3,TBX2,BAMBI,FGF8,CYP26C1,ALDH1A3,CYP1B1,FOXG1,HMX1,FOXD1"

DEFAULT_PARAMS = dict(
    pole=(0.0, 0.0),               # HAA pole = score origin (balanced D/V & N/T module scores)
    symmetric=True,                # equalize N/T and D/V extents via per-direction scaling
    scale="p99",                   # used only when symmetric=False
    rho_max_dv_deg=64.0,           # cap half-angle toward dorsal/ventral
    rho_max_nt_deg=82.0,           # cap half-angle toward nasal/temporal (stretched wider)
    dewarp="arcsin",               # "arcsin" | "pow" | "none"
    pow_p=1.6,
    cut_angles_deg=(0, 90, 180, 270),  # slits at cardinals -> 4 lobes at diagonals (DN,DT,VT,VN)
    gap_gain=1.0,                 # deficit-mode rip width
    gap_mode="deficit",             # "linear" (visible V-notch from rip-start) | "deficit" (curvature-true)
    gap_frac=0.5,                  # linear mode: max gap as fraction of the half-gore width (at rim)
    rho_join_deg=14.0,             # keep lobes joined within this colatitude (solid HAA center) = rip-start
    missing_nasal_frac=0.13,        # per-layer outward extrusion depth of the empty "missing" nasal cap
    missing_layers=2,              # radial layers of missing nasal tiles (ventro-nasal gets +1)
    stretch_bumps=((180.0, 0.60, 50.0), (225.0, 0.85, 50.0), (45.0, 0.40, 50.0)),              # selective directional radial stretch: [(angle_deg, amplitude, width_deg), ...]
    rot_deg=0.0,
)

def _resolve_pole(dv, nt, pole):
    if pole == "median":
        return float(np.nanmedian(dv)), float(np.nanmedian(nt))
    return float(pole[0]), float(pole[1])

def _resolve_scale(x, y, scale):
    if scale == "p99":
        return (float(np.nanpercentile(np.abs(x), 99)) or 1.0,
                float(np.nanpercentile(np.abs(y), 99)) or 1.0)
    return float(scale[0]), float(scale[1])

def _sym_scale(v, q=97.0):
    """Scale + and - sides independently to their q-th percentile -> ~[-1,1] symmetric.
    Equalizes opposing petals (e.g. nasal vs temporal) despite asymmetric score ranges."""
    v = np.asarray(v, float)
    pos = v[v > 0]; neg = -v[v < 0]
    sp = (float(np.nanpercentile(pos, q)) if pos.size else 1.0) or 1.0
    sn = (float(np.nanpercentile(neg, q)) if neg.size else 1.0) or 1.0
    return np.where(v >= 0, v / sp, v / sn)

def flower_transform(dv, nt, **p):
    """Map (DV.Score, NT.Score) -> (X, Y). Returns (X, Y, diag) with per-point gore/rho."""
    P = {**DEFAULT_PARAMS, **p}
    dv = np.asarray(dv, float); nt = np.asarray(nt, float)
    dv0, nt0 = _resolve_pole(dv, nt, P["pole"])
    x = nt - nt0          # NT -> horizontal (nasal +x, right)
    y = dv - dv0          # DV -> vertical   (dorsal +y, up)
    if P.get("symmetric", True):
        xs, ys = _sym_scale(x), _sym_scale(y)
        sx = sy = None
    else:
        sx, sy = _resolve_scale(x, y, P["scale"])
        xs, ys = x / sx, y / sy
    th = np.arctan2(ys, xs) + np.deg2rad(P["rot_deg"])
    r = np.clip(np.hypot(xs, ys), 0.0, 1.0)

    w_nt = np.cos(th) ** 2
    rho_max = np.deg2rad(P["rho_max_nt_deg"]) * w_nt + np.deg2rad(P["rho_max_dv_deg"]) * (1.0 - w_nt)
    if P["dewarp"] == "arcsin":
        rho = np.arcsin(np.clip(r * np.sin(rho_max), -1.0, 1.0))
    elif P["dewarp"] == "pow":
        rho = rho_max * (r ** P["pow_p"])
    else:
        rho = rho_max * r
    R = rho.copy()

    cuts = np.deg2rad(np.sort(np.asarray(P["cut_angles_deg"], float) % 360.0))
    K = len(cuts)
    phi = th % (2 * np.pi)
    gore = np.zeros_like(r, dtype=int)
    if K > 0:
        rj = np.deg2rad(P["rho_join_deg"])
        if P.get("gap_mode", "linear") == "deficit":
            with np.errstate(invalid="ignore", divide="ignore"):
                deficit = 2 * np.pi * (1.0 - np.where(rho > 1e-9, np.sin(rho) / rho, 1.0))
            deficit = np.clip(deficit * P["gap_gain"], 0.0, 2 * np.pi * 0.9)
            delta = deficit / (2 * K)
        else:  # "linear": V-notch that opens from the rip-start (rho_join) out to the rim
            half_gore = np.pi / K
            frac = np.clip((rho - rj) / np.maximum(rho_max - rj, 1e-6), 0.0, 1.0)
            delta = P.get("gap_frac", 0.5) * half_gore * frac
        thm = th % (2 * np.pi)
        thm = np.where(thm < cuts[0], thm + 2 * np.pi, thm)
        edges = np.concatenate([cuts, [cuts[0] + 2 * np.pi]])
        gore = np.clip(np.searchsorted(edges, thm, side="right") - 1, 0, K - 1)
        c0 = edges[gore]; c1 = edges[gore + 1]
        width = c1 - c0
        t = (thm - c0) / np.where(width > 0, width, 1.0)
        span = np.clip(width - 2 * delta, 1e-3, None)
        phi = c0 + delta + t * span
    # selective directional stretch (anatomical direction = th): elongate chosen sectors
    if P.get("stretch_bumps"):
        st = np.ones_like(R)
        for ang, amp, wid in P["stretch_bumps"]:
            d = np.angle(np.exp(1j * (th - np.deg2rad(ang))))   # signed angular diff in (-pi, pi]
            st = st + amp * np.exp(-0.5 * (d / np.deg2rad(wid)) ** 2)
        R = R * st
    X = R * np.cos(phi); Y = R * np.sin(phi)
    return X, Y, dict(r=r, rho=rho, phi=phi, gore=gore, pole=(dv0, nt0), scale=(sx, sy))

def _bin_sum_count(nt, dv, vals, nt_edges, dv_edges):
    ni = len(nt_edges) - 1; nj = len(dv_edges) - 1
    ix = np.clip(np.digitize(nt, nt_edges) - 1, 0, ni - 1)
    iy = np.clip(np.digitize(dv, dv_edges) - 1, 0, nj - 1)
    flat = iy * ni + ix
    cnt = np.bincount(flat, minlength=ni * nj).astype(float)
    ssum = np.bincount(flat, weights=np.nan_to_num(vals), minlength=ni * nj)
    return ssum.reshape(nj, ni), cnt.reshape(nj, ni)   # [dv, nt]

def _smooth_masked(img, mask, sigma):
    from scipy.ndimage import gaussian_filter
    filled = np.where(mask, 0.0, np.nan_to_num(img))
    wsm = gaussian_filter((~mask).astype(float), sigma)
    vsm = gaussian_filter(filled, sigma)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(wsm > 1e-6, vsm / wsm, np.nan)
    out[mask] = np.nan
    return out

def _pstr(P):
    Q = {**DEFAULT_PARAMS, **P}
    gap = f"{Q['gap_frac']}" if Q.get('gap_mode', 'linear') == 'linear' else f"{Q['gap_gain']}"
    s = (f"rho_max(DV/NT)={Q['rho_max_dv_deg']:.0f}/{Q['rho_max_nt_deg']:.0f} "
         f"cuts={list(Q['cut_angles_deg'])} join={Q['rho_join_deg']:.0f} gap={Q.get('gap_mode','linear')}:{gap}")
    if Q.get('stretch_bumps'):
        s += " stretch=" + str([(int(a), amp) for a, amp, _ in Q['stretch_bumps']])
    if Q.get('missing_nasal_frac', 0) > 0:
        s += f" missNasal={Q['missing_nasal_frac']}"
    return s

def render_binned(df, genes, params, out_png, nbins=64, min_cells=6, smooth_sigma=1.3, title="", show_title=True):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    dv = pd.to_numeric(df["DV.Score"], errors="coerce").to_numpy()
    nt = pd.to_numeric(df["NT.Score"], errors="coerce").to_numpy()
    nt_edges = np.linspace(np.nanmin(nt), np.nanmax(nt), nbins + 1)
    dv_edges = np.linspace(np.nanmin(dv), np.nanmax(dv), nbins + 1)

    NTv, DVv = np.meshgrid(nt_edges, dv_edges)               # (nbins+1, nbins+1) vertices
    Xv, Yv, diagv = flower_transform(DVv.ravel(), NTv.ravel(), **params)
    Xv = Xv.reshape(NTv.shape); Yv = Yv.reshape(NTv.shape)
    gore_v = diagv["gore"].reshape(NTv.shape)
    rho_v = diagv["rho"].reshape(NTv.shape)
    rho_join = np.deg2rad({**DEFAULT_PARAMS, **params}["rho_join_deg"])

    _, cnt = _bin_sum_count(nt, dv, np.ones_like(dv), nt_edges, dv_edges)
    mask = cnt < min_cells

    ng = len(genes)
    fig, axes = plt.subplots(2, ng, figsize=(4.1 * ng, 8.6), squeeze=False)
    for j, g in enumerate(genes):
        vals = pd.to_numeric(df[g], errors="coerce").to_numpy()
        ssum, _ = _bin_sum_count(nt, dv, vals, nt_edges, dv_edges)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = ssum / cnt
        img = _smooth_masked(mean, mask, smooth_sigma)
        fin = np.isfinite(img)
        vmax = np.nanpercentile(img[fin], 97) if fin.any() else 1.0
        vmax = vmax if vmax > 0 else 1.0
        norm = mpl.colors.Normalize(0, vmax); cmap = mpl.cm.magma

        ax = axes[0][j]
        ax.pcolormesh(nt_edges, dv_edges, np.ma.masked_invalid(img), cmap=cmap, norm=norm, shading="flat")
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{g}\n(NT, DV) square", fontsize=10)
        if j == 0:
            ax.set_xlabel("NT  (T <- . -> N)", fontsize=8); ax.set_ylabel("DV  (V <- . -> D)", fontsize=8)

        ax = axes[1][j]
        polys, cols = [], []
        for iy in range(nbins):
            for ix in range(nbins):
                if mask[iy, ix] or not np.isfinite(img[iy, ix]):
                    continue
                gs = gore_v[iy:iy + 2, ix:ix + 2]
                if gs.min() != gs.max() and rho_v[iy:iy + 2, ix:ix + 2].min() > rho_join:
                    continue   # quad crosses a rip beyond the joined center -> drop (slit)
                polys.append(np.column_stack([
                    [Xv[iy, ix], Xv[iy, ix + 1], Xv[iy + 1, ix + 1], Xv[iy + 1, ix]],
                    [Yv[iy, ix], Yv[iy, ix + 1], Yv[iy + 1, ix + 1], Yv[iy + 1, ix]]]))
                cols.append(img[iy, ix])
        pc = PolyCollection(polys, array=np.asarray(cols), cmap=cmap, norm=norm, edgecolors="none")
        ax.add_collection(pc); ax.autoscale_view(); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title("flower / orange-peel", fontsize=10, pad=12)
        for txt, (fx, fy, ha, va) in {"N": (0.99, 0.5, "right", "center"),
                                      "T": (0.01, 0.5, "left", "center"),
                                      "D": (0.5, 0.90, "center", "top"),
                                      "V": (0.5, 0.02, "center", "bottom")}.items():
            ax.annotate(txt, (fx, fy), xycoords="axes fraction", ha=ha, va=va, fontsize=9, color="0.35")

    if show_title:
        fig.suptitle((title + "  |  " if title else "") + _pstr(params), fontsize=10)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
    else:
        fig.tight_layout(rect=[0, 0, 1, 0.99])
    fig.savefig(out_png, dpi=150)
    print("wrote", out_png)

def render_grid(df, genes, params, out_png, ncols=4, nbins=64, min_cells=6, smooth_sigma=1.3,
                title="", ref_density=False, gene_inset=True, show_title=True):
    """Flower-only grid for many genes. Warp geometry is gene-independent -> computed once."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    dv = pd.to_numeric(df["DV.Score"], errors="coerce").to_numpy()
    nt = pd.to_numeric(df["NT.Score"], errors="coerce").to_numpy()
    frac = {**DEFAULT_PARAMS, **params}.get("missing_nasal_frac", 0.0)
    nt_edges = np.linspace(np.nanmin(nt), np.nanmax(nt), nbins + 1)
    dv_edges = np.linspace(np.nanmin(dv), np.nanmax(dv), nbins + 1)
    NTv, DVv = np.meshgrid(nt_edges, dv_edges)
    Xv, Yv, diagv = flower_transform(DVv.ravel(), NTv.ravel(), **params)
    Xv = Xv.reshape(NTv.shape); Yv = Yv.reshape(NTv.shape)
    gore_v = diagv["gore"].reshape(NTv.shape); rho_v = diagv["rho"].reshape(NTv.shape)
    rho_join = np.deg2rad({**DEFAULT_PARAMS, **params}["rho_join_deg"])
    _, cnt = _bin_sum_count(nt, dv, np.ones_like(dv), nt_edges, dv_edges)
    mask = cnt < min_cells

    polys, ij = [], []
    for iy in range(nbins):
        for ix in range(nbins):
            if mask[iy, ix]:
                continue
            gs = gore_v[iy:iy + 2, ix:ix + 2]
            if gs.min() != gs.max() and rho_v[iy:iy + 2, ix:ix + 2].min() > rho_join:
                continue
            polys.append(np.column_stack([
                [Xv[iy, ix], Xv[iy, ix + 1], Xv[iy + 1, ix + 1], Xv[iy + 1, ix]],
                [Yv[iy, ix], Yv[iy, ix + 1], Yv[iy + 1, ix + 1], Yv[iy + 1, ix]]]))
            ij.append((iy, ix))
    ij = np.asarray(ij); iy_a, ix_a = ij[:, 0], ij[:, 1]

    ghost_polys = []                                # empty "missing/uncaptured" most-nasal bins
    if frac > 0:                                    # extrude nasal-facing outer-edge data bins outward
        nlayers = int({**DEFAULT_PARAMS, **params}.get("missing_layers", 2))
        dataset = set(map(tuple, ij.tolist()))
        for (iy, ix) in dataset:
            if (iy, ix + 1) in dataset:             # has a nasal neighbour -> not the nasal outer edge
                continue
            cx = 0.25 * (Xv[iy, ix] + Xv[iy, ix + 1] + Xv[iy + 1, ix] + Xv[iy + 1, ix + 1])
            cy = 0.25 * (Yv[iy, ix] + Yv[iy, ix + 1] + Yv[iy + 1, ix] + Yv[iy + 1, ix + 1])
            if cx <= 0 or cy > cx or cy < -1.4 * cx:  # nasal side; dorso<=45deg, ventro<=~54deg (wider VN)
                continue
            L = nlayers + (1 if cy < 0 else 0)      # one extra layer ventro-nasally
            p1 = np.array([Xv[iy, ix + 1], Yv[iy, ix + 1]])
            p2 = np.array([Xv[iy + 1, ix + 1], Yv[iy + 1, ix + 1]])
            for k in range(L):
                f0, f1 = 1.0 + frac * k, 1.0 + frac * (k + 1)
                ghost_polys.append(np.array([p1 * f0, p1 * f1, p2 * f1, p2 * f0]))

    ng = len(genes); nrows = int(np.ceil(ng / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.7 * ncols, 2.9 * nrows), squeeze=False)
    for k, g in enumerate(genes):
        ax = axes[k // ncols][k % ncols]
        vals = pd.to_numeric(df[g], errors="coerce").to_numpy()
        ssum, _ = _bin_sum_count(nt, dv, vals, nt_edges, dv_edges)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = ssum / cnt
        img = _smooth_masked(mean, mask, smooth_sigma)
        fin = np.isfinite(img)
        vmax = np.nanpercentile(img[fin], 97) if fin.any() else 1.0
        vmax = vmax if vmax > 0 else 1.0
        norm = mpl.colors.Normalize(0, vmax)
        pc = PolyCollection(polys, array=np.ma.masked_invalid(img[iy_a, ix_a]),
                            cmap=mpl.cm.magma, norm=norm, edgecolors="none")
        ax.add_collection(pc)
        if ghost_polys:  # empty outlined bins = anatomically present but uncaptured most-nasal retina
            gpc = PolyCollection(ghost_polys, facecolors="none", edgecolors="0.4", linewidths=0.6)
            gpc.set_hatch("///")
            ax.add_collection(gpc)
        ax.autoscale_view(); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title(g, fontsize=11, style="italic")
        if gene_inset:  # small non-stretched (original square NT/DV) map of the same gene
            axin = ax.inset_axes([0.0, 0.72, 0.28, 0.28])
            axin.pcolormesh(nt_edges, dv_edges, np.ma.masked_invalid(img),
                            cmap=mpl.cm.magma, norm=norm, shading="flat")
            axin.set_aspect("equal"); axin.set_xticks([]); axin.set_yticks([])
            for s in axin.spines.values():
                s.set_linewidth(0.5); s.set_edgecolor("0.55")
            if k == 0:
                axin.set_title("orig.", fontsize=6, pad=1)
    for k in range(ng, nrows * ncols):
        axes[k // ncols][k % ncols].set_visible(False)
    for txt, (fx, fy, ha, va) in {"N": (0.99, 0.5, "right", "center"),
                                  "T": (0.01, 0.5, "left", "center"),
                                  "D": (0.5, 0.94, "center", "top"),
                                  "V": (0.5, 0.04, "center", "bottom")}.items():
        axes[0][0].annotate(txt, (fx, fy), xycoords="axes fraction", ha=ha, va=va, fontsize=9, color="0.4")
    if show_title:
        fig.suptitle((title + "  |  " if title else "") + _pstr(params), fontsize=9, y=0.995)
        fig.tight_layout(rect=[0, 0, 1, 0.90])
    else:
        fig.tight_layout(rect=[0, 0, 1, 0.985])
    if ref_density:
        axr = fig.add_axes([0.015, 0.905, 0.13, 0.085])  # small reference: un-warped square density
        dens = np.where(mask, np.nan, cnt)
        with np.errstate(invalid="ignore"):
            axr.pcolormesh(nt_edges, dv_edges, np.ma.masked_invalid(np.log10(dens)),
                           cmap="viridis", shading="flat")
        axr.set_aspect("equal"); axr.set_xticks([]); axr.set_yticks([])
        for s in axr.spines.values():
            s.set_linewidth(0.4)
        axr.set_title("original (NT,DV)\ncell density", fontsize=6.5)
    fig.savefig(out_png, dpi=150)
    print("wrote", out_png)

def synthetic_df(n=14000, seed=0):
    rng = np.random.default_rng(seed)
    ang = rng.uniform(0, 2 * np.pi, n)
    rad = np.sin(np.sqrt(rng.uniform(0, 1, n)) * (np.pi / 2))   # crowd the rim (saturation)
    nt = rad * np.cos(ang); dv = rad * np.sin(ang)
    df = pd.DataFrame({"DV.Score": dv, "NT.Score": nt})
    df["NASAL_grad"] = np.clip(nt, 0, None)
    df["DORSAL_grad"] = np.clip(dv, 0, None)
    df["HAA_center"] = np.exp(-(nt ** 2 + dv ** 2) / 0.05)
    return df, ["HAA_center", "NASAL_grad", "DORSAL_grad"]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", type=Path, default=None, help="optional pre-extracted DV/NT and gene table")
    ap.add_argument("--h5ad", type=Path, default=REPO / "data/20250604_chick_RPC.h5ad")
    ap.add_argument("--genes", default=PUBLISHED_GENES)
    ap.add_argument("--out", type=Path, default=REPO / "figures/Figure_SF15/SF15C_flower_reprojection.png")
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--nbins", type=int, default=64)
    ap.add_argument("--min_cells", type=int, default=6)
    ap.add_argument("--smooth", type=float, default=1.3)
    ap.add_argument("--rho_nt", type=float, default=DEFAULT_PARAMS["rho_max_nt_deg"])
    ap.add_argument("--rho_dv", type=float, default=DEFAULT_PARAMS["rho_max_dv_deg"])
    ap.add_argument("--gap_gain", type=float, default=DEFAULT_PARAMS["gap_gain"])
    ap.add_argument("--rho_join", type=float, default=DEFAULT_PARAMS["rho_join_deg"])
    ap.add_argument("--cuts", default="0,90,180,270")
    ap.add_argument("--dewarp", default="arcsin")
    ap.add_argument("--grid", action=argparse.BooleanOptionalAction, default=True, help="flower-only grid (published layout)")
    ap.add_argument("--ncols", type=int, default=4)
    ap.add_argument("--gap_mode", default=DEFAULT_PARAMS["gap_mode"], choices=["linear", "deficit"])
    ap.add_argument("--gap_frac", type=float, default=DEFAULT_PARAMS["gap_frac"])
    ap.add_argument("--stretch_t", type=float, default=0.60, help="temporal (-NT, 180deg) radial stretch amplitude")
    ap.add_argument("--stretch_vt", type=float, default=0.85, help="ventro-temporal (225deg) radial stretch amplitude")
    ap.add_argument("--stretch_dn", type=float, default=0.40, help="dorso-nasal (45deg) radial stretch amplitude")
    ap.add_argument("--stretch_width", type=float, default=50.0, help="angular width (deg) of stretch bumps")
    ap.add_argument("--missing_nasal", type=float, default=0.13, help="per-layer depth of empty 'missing' nasal cap (e.g. 0.13)")
    ap.add_argument("--missing_layers", type=int, default=2, help="radial layers of missing nasal tiles (ventro-nasal +1)")
    ap.add_argument("--no_title", action=argparse.BooleanOptionalAction, default=True, help="suppress the parameter suptitle (published layout)")
    args = ap.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    _bumps = []
    if args.stretch_t:
        _bumps.append((180.0, args.stretch_t, args.stretch_width))
    if args.stretch_vt:
        _bumps.append((225.0, args.stretch_vt, args.stretch_width))
    if args.stretch_dn:
        _bumps.append((45.0, args.stretch_dn, args.stretch_width))
    params = dict(rho_max_nt_deg=args.rho_nt, rho_max_dv_deg=args.rho_dv,
                  gap_gain=args.gap_gain, gap_mode=args.gap_mode, gap_frac=args.gap_frac,
                  rho_join_deg=args.rho_join, dewarp=args.dewarp, stretch_bumps=tuple(_bumps),
                  missing_nasal_frac=args.missing_nasal, missing_layers=args.missing_layers,
                  cut_angles_deg=tuple(float(c) for c in args.cuts.split(",") if c.strip()))
    if args.synthetic:
        df, genes = synthetic_df()
        render_binned(df, genes, params, args.out, nbins=args.nbins,
                      min_cells=args.min_cells, smooth_sigma=args.smooth, title="SYNTHETIC",
                      show_title=not args.no_title)
        return
    if args.tsv is None:
        from sfig15c_extract_dvnt import extract_dataframe
        df = extract_dataframe(args.h5ad)
    else:
        df = pd.read_csv(args.tsv, sep="\t", index_col=0)
    genes = [g for g in args.genes.split(",") if g in df.columns]
    missing = [g for g in args.genes.split(",") if g not in df.columns]
    if missing:
        print(f"WARNING: requested genes absent from TSV (skipped): {missing}")
    if not genes:
        raise SystemExit("none of the requested genes present in TSV")
    if args.grid or len(genes) > 3:
        render_grid(df, genes, params, args.out, ncols=args.ncols, nbins=args.nbins,
                    min_cells=args.min_cells, smooth_sigma=args.smooth, title="chick RPC",
                    show_title=not args.no_title)
    else:
        render_binned(df, genes, params, args.out, nbins=args.nbins,
                      min_cells=args.min_cells, smooth_sigma=args.smooth, title="chick RPC",
                      show_title=not args.no_title)

if __name__ == "__main__":
    main()
