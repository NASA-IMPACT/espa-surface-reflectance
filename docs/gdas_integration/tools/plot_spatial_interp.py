#!/usr/bin/env python3
"""How LaSRC turns four 1-degree GRIB values into the value of one 0.05-degree CMG cell (step 3), and what would change
if the bilinear step were replaced by the nearest GRIB value.

Left: the 4 GRIB grid points around one scene, the 0.05-degree CMG cells LaSRC fills by bilinear interpolation,
the CMG cell at the scene centre and the bilinear weights (same code path as scene_centre_values.py / lut_subr.c).
Right: over all land CMG cells (every 4th row/column), % difference between the nearest GRIB value and the bilinear
value, and how often the stored integer (wv = int(PWAT*10), oz = int(TOZNE*0.4)) would change.

Usage: plot_spatial_interp.py --record REC.json --pair F1 F2 --anl F [F ...] --name NAME --out DIR
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import scene_centre_values as S  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.path import Path  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

INK, INK2, GRID, RED = "#0b0b0b", "#52514e", "#d9d8d4", "#d62728"
SEQ = LinearSegmentedColormap.from_list("seq", ["#d6e6f8", "#86b3e6", "#2a78d6", "#0d3f86"])
plt.rcParams.update({"font.size": 6.5, "text.color": INK})
SHIFT = 1


def lasrc_grid(fields, var):
    return np.roll(fields[var].astype(np.float32), -(180 + SHIFT), axis=1)   # row i: lat 90-i, col j: lon j-180


def linesamp(lat, lon):
    lat = np.asarray(lat, dtype=np.float32); lon = np.asarray(lon, dtype=np.float32)
    il = ((90 - lat) + 0.5).astype(int); rf = (90 - il) - lat
    neg = rf < 0; il = np.where(neg, il - 1, il); rf = np.where(neg, 1 - np.abs(rf), rf)
    lon = np.where(lon > 179.5, -180 - (180 - lon), lon)
    js = ((lon + 180) + 0.5).astype(int); cf = lon + (180 - js)
    neg = cf < 0; js = np.where(neg, js - 1, js) % 360; cf = np.where(neg, 1 - np.abs(cf), cf)
    return il, js, rf, cf


def bilinear(g, il, js, rf, cf):
    j1 = (js + 1) % 360; i1 = np.minimum(il + 1, 180)
    return ((1 - rf) * (1 - cf) * g[il, js] + rf * (1 - cf) * g[i1, js] + rf * cf * g[i1, j1]
            + cf * (1 - rf) * g[il, j1])


def nearest(g, il, js, rf, cf):
    return g[np.where(rf < 0.5, il, np.minimum(il + 1, 180)), np.where(cf < 0.5, js, (js + 1) % 360)]


def land_mask(lat, lon):
    rings = json.load(open(os.path.join(HERE, "world_land_110m.json")))["rings"]
    pts = np.column_stack([lon.ravel(), lat.ravel()]); m = np.zeros(len(pts), bool)
    for r in rings:
        r = np.asarray(r)
        if len(r) < 3:
            continue
        bb = (pts[:, 0] >= r[:, 0].min()) & (pts[:, 0] <= r[:, 0].max()) & (pts[:, 1] >= r[:, 1].min()) & \
             (pts[:, 1] <= r[:, 1].max())
        if bb.any():
            idx = np.where(bb)[0]
            m[idx[Path(r).contains_points(pts[idx])]] = True
    return m.reshape(lat.shape)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--record", required=True)
    ap.add_argument("--pair", nargs=2, required=True, help="the two V0 files of the example scene")
    ap.add_argument("--anl", nargs="+", required=True, help="analysis files for the global comparison")
    ap.add_argument("--name", default="")
    ap.add_argument("--tile", default="", help="MGRS tile to outline")
    ap.add_argument("--wv-range", nargs=2, type=float, default=[0, 20])
    ap.add_argument("--oz-range", nargs=2, type=float, default=[250, 470])
    ap.add_argument("--no-hist", action="store_true", help="only the interpolation panels (no global comparison)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    r = json.load(open(a.record))
    w1, w2, lat, lon = r["first_weight"], r["second_weight"], r["scene_lat"], r["scene_lon"]
    A, B = S.read_fields(a.pair[0]), S.read_fields(a.pair[1])
    lcmg = round((89.975 - lat) * 20); scmg = round((179.975 + lon) * 20)
    clat, clon = 90 - lcmg * 0.05, -180 + scmg * 0.05            # CMG cell LaSRC reads for this scene
    il, js, rf, cf = [x.item() for x in linesamp(clat, clon)]

    poly = None
    if a.tile:
        import mgrs
        from pyproj import Transformer
        z, h, e, n = mgrs.MGRS().MGRSToUTM(a.tile + "0000000000")
        tr = Transformer.from_crs((32600 if h == "N" else 32700) + int(z), 4326, always_xy=True)
        cs = [tr.transform(e + de, n + dn) for de, dn in [(-40, 100020), (109760, 100020), (109760, -9780), (-40, -9780), (-40, 100020)]]
        poly = (np.array([c[0] for c in cs]), np.array([c[1] for c in cs]))
    CELL = "#ff8c00"
    # ---- global comparison: nearest GRIB value vs bilinear, land CMG cells, every 4th row / column
    glat = (90 - np.arange(0, 3600, 4) * 0.05)[:, None] * np.ones((1, 1800))
    glon = np.ones((900, 1)) * (-180 + np.arange(0, 7200, 4) * 0.05)[None, :]
    keep = land_mask(glat, glon) & (glat > -60) & (glat < 75)
    gi, gj, grf, gcf = linesamp(glat[keep], glon[keep])
    res = {"PWAT": [], "TOZNE": []}
    for f in a.anl:
        F = S.read_fields(f)
        for var, mult in (("PWAT", 10), ("TOZNE", 0.4)):
            g = lasrc_grid(F, var)
            b = bilinear(g, gi, gj, grf, gcf); n = nearest(g, gi, gj, grf, gcf)
            ok = b > (1.0 if var == "PWAT" else 50)
            res[var].append((100 * (n[ok] - b[ok]) / b[ok], (n[ok] * mult).astype(int) != (b[ok] * mult).astype(int)))

    nh = a.no_hist
    if nh:
        fig = plt.figure(figsize=(7.3, 3.7), facecolor="white")
        gs = fig.add_gridspec(1, 4, width_ratios=[1.0, 0.72, 1.0, 0.72], wspace=0.08, left=0.04, right=0.99,
                              top=0.80, bottom=0.12)
    else:
        fig = plt.figure(figsize=(7.3, 4.3), facecolor="white")
        gs = fig.add_gridspec(2, 3, width_ratios=[0.78, 0.62, 1.15], wspace=0.12, hspace=0.62, left=0.04,
                              right=0.985, top=0.80, bottom=0.12)
    meta = {"PWAT": ("Water vapour (kg/m²)", "%.2f", 10, 0.1, "kg/m²"), "TOZNE": ("Total ozone (DU)", "%.1f", 0.4, 2.5, "DU")}
    summary = {}
    for r_, var in enumerate(("PWAT", "TOZNE")):
        label, fmt, mult, step, unit = meta[var]
        g = lasrc_grid(A, var) * np.float32(w1) + lasrc_grid(B, var) * np.float32(w2)
        # zoom: the GRIB cell containing the scene, CMG cells at 0.05 deg
        la0, lo0 = 90 - il, js - 180                                 # upper-left grid point
        x0, x1, y0, y1 = lo0 - 0.62, lo0 + 1.62, la0 - 1.62, la0 + 0.62
        if poly is not None:
            x0, x1 = min(x0, poly[0].min() - 0.08), max(x1, poly[0].max() + 0.08)
            y0, y1 = min(y0, poly[1].min() - 0.08), max(y1, poly[1].max() + 0.08)
        ys = la0 + 0.05 * np.arange(int(np.floor((y0 - la0) / 0.05)), int(np.ceil((y1 - la0) / 0.05)) + 1)[::-1]
        xs = lo0 + 0.05 * np.arange(int(np.floor((x0 - lo0) / 0.05)), int(np.ceil((x1 - lo0) / 0.05)) + 1)
        YY, XX = np.meshgrid(ys, xs, indexing="ij")
        zi, zj, zr, zc = linesamp(YY, XX)
        Z = bilinear(g, zi, zj, zr, zc)
        ax = fig.add_subplot(gs[0, 2 * r_] if nh else gs[r_, 0])
        vr = a.wv_range if var == "PWAT" else a.oz_range
        im = ax.pcolormesh(xs, ys, Z, cmap=SEQ, shading="nearest", edgecolors=(1, 1, 1, 0.35), linewidth=0.15,
                           vmin=vr[0], vmax=vr[1])
        corners = [(il, js, lo0, la0), (il + 1, js, lo0, la0 - 1), (il + 1, (js + 1) % 360, lo0 + 1, la0 - 1),
                   (il, (js + 1) % 360, lo0 + 1, la0)]
        wts = [(1 - rf) * (1 - cf), rf * (1 - cf), rf * cf, cf * (1 - rf)]
        ax.plot([lo0, lo0 + 1, lo0 + 1, lo0, lo0], [la0, la0, la0 - 1, la0 - 1, la0], color=INK, lw=0.6,
                ls=(0, (3, 2)), zorder=5)
        for k in (lo0 - 0.5, lo0 + 0.5, lo0 + 1.5):
            ax.axvline(k, color="white", lw=0.9, zorder=4)
        for k in (la0 + 0.5, la0 - 0.5, la0 - 1.5):
            ax.axhline(k, color="white", lw=0.9, zorder=4)
        lines = []
        for k, ((ii, jj, x, y), w) in enumerate(zip(corners, wts), 1):
            ax.plot(x, y, "s", ms=5, mfc="white", mec=INK, mew=0.9, clip_on=False, zorder=9)
            ax.plot([x, clon], [y, clat], color=INK, lw=0.5, ls=(0, (1.5, 1.5)), zorder=5)
            ax.text(x + (-0.07 if x == lo0 else 0.07), y + (0.07 if y == la0 else -0.07), "v%d" % k, fontsize=5.6,
                    ha="right" if x == lo0 else "left", va="bottom" if y == la0 else "top", fontweight="bold",
                    zorder=9, bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.85))
            lines.append("v%d = %s   w = %.3f" % (k, fmt % g[ii, jj], w))
        ax.add_patch(plt.Rectangle((clon - 0.025, clat - 0.025), 0.05, 0.05, fill=False, ec=CELL, lw=1.3, zorder=7))
        if poly is not None:
            ax.plot(poly[0], poly[1], color=RED, lw=1.0, zorder=6)
        ax.plot(lon, lat, "o", ms=2.4, mfc="white", mec=INK, mew=0.6, zorder=8)
        val = bilinear(g, *[np.array(v) for v in (il, js, rf, cf)]).item()
        tx = fig.add_subplot(gs[0, 2 * r_ + 1] if nh else gs[r_, 1]); tx.axis("off")
        tx.text(0.02, 0.84, "\n".join(lines), fontsize=6, va="top", family="DejaVu Sans Mono", linespacing=1.5)
        tx.text(0.02, 0.36, "Σ w·v = %s\n(orange CMG cell)\nstored after truncation:\n%s %s" %
                (fmt % val, fmt % (int(val * mult) / mult), unit), fontsize=6.2, va="top", color="#b35900",
                fontweight="bold", linespacing=1.4)
        tx.text(0.02, 0.97, "row r = %.2f, column c = %.2f" % (rf, cf), fontsize=5.8, color=INK2, va="top")
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect(1 / np.cos(np.radians(lat)))
        ax.set_xticks([lo0, lo0 + 1]); ax.set_yticks([la0 - 1, la0])
        ax.set_anchor("W")
        ax.set_xticklabels(["%d°" % lo0, "%d°" % (lo0 + 1)], fontsize=5.5)
        ax.set_yticklabels(["%d°" % (la0 - 1), "%d°" % la0], fontsize=5.5)
        ax.tick_params(length=1.5, pad=1)
        for sp in ax.spines.values():
            sp.set_linewidth(0.5)
        ax.set_title("%s: 4 GRIB points (□) → 0.05° CMG cells" % label.split(" (")[0], fontsize=6.6,
                     fontweight="bold", pad=9, loc="left")
        cb = (fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.05, pad=0.08, shrink=0.8, extend="both") if nh else
              fig.colorbar(im, ax=ax, fraction=0.06, pad=0.03, shrink=0.9)); cb.ax.tick_params(labelsize=5, length=1.5)
        cb.outline.set_linewidth(0.3); cb.set_label(unit, fontsize=5.5, labelpad=1)
        if nh:
            continue
        # global histogram
        d = np.concatenate([x[0] for x in res[var]]); ch = np.concatenate([x[1] for x in res[var]])
        summary[var] = (np.median(np.abs(d)), np.percentile(np.abs(d), 95), 100 * ch.mean())
        ax = fig.add_subplot(gs[r_, 2])
        lim = np.percentile(np.abs(d), 99.5)
        ax.hist(d, bins=81, range=(-lim, lim), color="#2a78d6", alpha=0.85, lw=0)
        ax.axvline(0, color=INK, lw=0.6)
        ax.set_yticks([]); ax.tick_params(axis="x", labelsize=5.5, length=2)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.set_xlabel("nearest GRIB value − bilinear value (% of bilinear; central 99.5% shown)", fontsize=5.8, labelpad=2)
        ax.set_title("Nearest GRIB point instead of bilinear (%s)" % label.split(" (")[0].lower(), fontsize=6.6,
                     fontweight="bold", pad=3)
        ax.text(0.99, 0.95, "median |Δ| %.1f%%, 95th pct %.1f%%\nstored value changes in %.0f%% of cells" %
                summary[var], transform=ax.transAxes, ha="right", va="top", fontsize=6, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=GRID))
    fig.text(0.05 if not nh else 0.02, 0.975, ("Step 3: LaSRC fills a 0.05° grid (CMG) from the 1° GRIB grid by bilinear "
             "interpolation,%sthen reads the one cell at the scene centre") % " ", fontsize=7.6,
             fontweight="bold", va="top")
    sub = ("Left: %s (V0, time-interpolated); white lines = 1° GRIB cells (one per point v1–v4), red box = tile, dot = scene centre, orange = CMG cell LaSRC reads.\nRight: all land CMG cells 60°S–75°N (every 4th), %d GDAS "
             "analyses (Jan and Jul 2025)." % (a.name, len(a.anl)))
    if nh:
        sub = ("%s (V0, time-interpolated). White lines = 1° GRIB cells (one per point v1–v4); red box = tile;\n"
               "dot = scene centre; orange = the CMG cell LaSRC reads." % a.name)
    fig.text(0.05 if not nh else 0.02, 0.935 if not nh else 0.925, sub, fontsize=6, color=INK2, va="top")
    fig.text(0.05 if not nh else 0.02, 0.0, "value = (1−r)(1−c)·v₁ + r(1−c)·v₂ + r·c·v₃ + c(1−r)·v₄,  r, c = fractional row / "
             "column position of the CMG cell between the 4 GRIB points (lut_subr.c: get_cart_interp)",
             fontsize=6.2, color=INK, va="bottom")
    os.makedirs(a.out, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(a.out, "slide_spatial_interp%s.%s" % ("" if not nh else "_nohist", ext)), dpi=300, bbox_inches="tight", pad_inches=0.04)
    print("weights", [round(w, 3) for w in [(1 - rf) * (1 - cf), rf * (1 - cf), rf * cf, cf * (1 - rf)]], summary)


if __name__ == "__main__":
    main()
