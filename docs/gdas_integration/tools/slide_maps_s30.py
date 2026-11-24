#!/usr/bin/env python3
"""S30 versions of the two L30 map slides (same layout and style as slide_figures.py):
  slide_vs_viirs_B03_S30   each run vs the VIIRS-based run, signed B03 mean difference (DN), pass/fail pies
  slide_change_vs_v0_B03_S30   S1-S3 change vs V0 in % of V0 B03 reflectance, plus a count table
Granule locations are the MGRS tile centres (= LaSRC's scene centre for S30).
Hollow grey on the change map = no change at all vs V0 (|mean difference| < 0.5e-6 reflectance).

Usage: slide_maps_s30.py --root ROOT [--stage immediate] [--band B03] --out DIR
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import slide_figures as F  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402
import mgrs  # noqa: E402
from pyproj import Transformer  # noqa: E402

_TR = {}


def tile_centre(tile):
    z, h, e, n = mgrs.MGRS().MGRSToUTM(tile + "0000000000")
    epsg = (32600 if h == "N" else 32700) + int(z)
    if epsg not in _TR:
        _TR[epsg] = Transformer.from_crs(epsg, 4326, always_xy=True)
    lon, lat = _TR[epsg].transform(e - 40 + 54900, n + 100020 - 54900)
    return lat, lon


def newest(d):
    f = sorted(x for x in glob.glob(os.path.join(d, "reports", "lasrc_validation_metrics_*.csv"))
               if not os.path.basename(x).startswith("._"))
    return pd.read_csv(f[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/rtmp/tvo/LaSRC_ancillary_latency_20261120")
    ap.add_argument("--stage", default="immediate")
    ap.add_argument("--band", default="B03")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    vf = os.path.join(a.root, "validation_full")
    viirs = {r: newest(os.path.join(vf, "S30_%s_vs_VIIRS_%s" % (r, a.stage))) for r in ("V0", "S1", "S2", "S3")}
    v0 = {s: newest(os.path.join(vf, "S30_%s_vs_V0_%s" % (s, a.stage))) for s in ("S1", "S2", "S3")}
    grans = set().union(*[set(m["granule"]) for m in list(viirs.values()) + list(v0.values())])
    loc = pd.DataFrame([dict(granule=g, tile=re.search(r"T(\d{2}[A-Z]{3})", g).group(1)) for g in grans]).set_index("granule")
    cent = {t: tile_centre(t) for t in loc["tile"].unique()}
    loc["scene_lat"] = [cent[t][0] for t in loc["tile"]]
    loc["scene_lon"] = [cent[t][1] for t in loc["tile"]]
    rings = json.load(open(os.path.join(HERE, "world_land_110m.json")))["rings"]
    b = a.band
    lab = {"B03": "Green", "B04": "Red", "B12": "SWIR2"}.get(b, b)

    # ---- vs VIIRS
    vmax = 40.0
    fig, axes, cb_y = F.map_figure("S30 %s %s: each run compared with the VIIRS-based run (signed mean difference, HLS DN)"
                                   % (b, lab))
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
    thr = None
    for ax, r in zip(axes, ("V0", "S1", "S2", "S3")):
        F.land(ax, rings)
        m = viirs[r][viirs[r]["band"] == b].set_index("granule")
        thr = float(m["threshold_refl"].iloc[0]) * 1e4
        dn = m["mean_diff_refl"] * 1e4
        fail = m["decision"].astype(str).str.upper().str.startswith("FAIL")
        ll = loc.loc[dn.index]
        o = np.argsort(np.abs(dn.to_numpy()))
        ax.scatter(ll["scene_lon"].to_numpy()[o], ll["scene_lat"].to_numpy()[o], s=F.size_of(dn.to_numpy()[o], vmax),
                   c=dn.clip(-vmax, vmax).to_numpy()[o], cmap=F.M.BRTE, norm=norm,
                   edgecolors=np.where(fail.to_numpy()[o], F.RED, F.INK2),
                   linewidths=np.where(fail.to_numpy()[o], 1.1, 0.3), zorder=3)
        ax.set_title("%s  |  median |diff| %.1f DN  |  fail %d of %d" % (F.LAB[r], dn.abs().median(), fail.sum(), len(dn)),
                     loc="left", fontsize=6.6, pad=2)
        F.pie(ax, int((~fail).sum()), int(fail.sum()))
    F.cbar(fig, norm, F.M.BRTE, 0.2, 0.8, cb_y,
           "difference from VIIRS run (HLS DN): teal = brighter, brown = darker; red outline = above the %.0f DN "
           "threshold" % thr)
    F.save(fig, a.out, "slide_vs_viirs_%s_S30" % b)

    # ---- change vs V0 (%)
    pan = {}
    for s, m in v0.items():
        m = m[m["band"] == b].set_index("granule")
        pan[s] = (100 * m["mean_diff_refl"] / m["mean_pre_refl"], m["mean_diff_refl"].abs() < 0.5e-6,
                  m["decision"].astype(str).str.upper().str.startswith("FAIL"))
    vmax = max(float(p.abs().max()) for p, _, _ in pan.values())
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
    fig, axes, cb_y = F.map_figure("S30 %s %s: change in surface reflectance vs V0 (%% of V0 reflectance)" % (b, lab))
    n = None
    for ax, s in zip(axes, ("S1", "S2", "S3")):
        F.land(ax, rings)
        p, same, fl = pan[s]
        n = len(p)
        ll = loc.loc[p.index]
        ax.scatter(ll.loc[same, "scene_lon"], ll.loc[same, "scene_lat"], s=6, facecolors="none", edgecolors="#9a9893",
                   linewidths=0.5, zorder=3)
        d = p[~same]
        o = np.argsort(np.abs(d.to_numpy()))
        ax.scatter(ll.loc[d.index, "scene_lon"].to_numpy()[o], ll.loc[d.index, "scene_lat"].to_numpy()[o],
                   s=F.size_of(d.to_numpy()[o], vmax), c=d.to_numpy()[o], cmap=F.M.BRTE, norm=norm,
                   edgecolors=np.where(fl.reindex(d.index).to_numpy()[o], F.RED, F.INK2),
                   linewidths=np.where(fl.reindex(d.index).to_numpy()[o], 1.1, 0.3), zorder=4)
        ax.set_title("%s vs V0  |  median %.2f%%, largest %.2f%%" % (F.LAB[s], p.abs().median(), p.abs().max()),
                     loc="left", fontsize=6.6, pad=2)
    ax = axes[3]
    ax.axis("off")
    rows = [["|change| above", "S1", "S2", "S3"]]
    for t_ in (0.10, 0.25, 0.50, 1.00):
        rows.append(["%.2f%%" % t_] + [str(int((pan[s][0].abs() > t_).sum())) for s in ("S1", "S2", "S3")])
    rows.append(["no change"] + [str(int(pan[s][1].sum())) for s in ("S1", "S2", "S3")])
    tb = ax.table(cellText=rows, loc="center", cellLoc="center", colWidths=[0.34, 0.16, 0.16, 0.16],
                  bbox=[0.12, 0.3, 0.76, 0.65])
    tb.auto_set_font_size(False); tb.set_fontsize(6.8)
    for (i, j), c in tb.get_celld().items():
        c.set_edgecolor(F.GRID); c.set_linewidth(0.5)
        if i == 0:
            c.set_text_props(fontweight="bold")
    ax.text(0.12, 0.2, "granules out of %d.  %s threshold (%.0f DN).\nHollow grey = no change vs V0; red outline = "
            "fails\nvs V0 in this band" % (n, b, float(v0["S1"][v0["S1"]["band"] == b]["threshold_refl"].iloc[0]) * 1e4),
            transform=ax.transAxes, fontsize=6.3, color=F.INK2, va="top")
    F.cbar(fig, norm, F.M.BRTE, 0.2, 0.8, cb_y, "change vs V0 (%): teal = brighter, brown = darker; bigger dot = bigger change")
    F.save(fig, a.out, "slide_change_vs_v0_%s_S30" % b)
    for s in ("S1", "S2", "S3"):
        p, same, fl = pan[s]
        print(s, "n=%d no-change=%d >0.25%%=%d >0.5%%=%d >1%%=%d max=%.2f fail=%d" %
              (len(p), same.sum(), (p.abs() > .25).sum(), (p.abs() > .5).sum(), (p.abs() > 1).sum(), p.abs().max(), fl.sum()))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
