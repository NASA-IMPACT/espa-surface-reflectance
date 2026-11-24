#!/usr/bin/env python3
"""One granule, four ancillary scenarios: the water vapour and ozone fields LaSRC reads (left) and how each
becomes the single scene-centre value LaSRC uses (right).

Fields: time-interpolated to the acquisition (LaSRC's weights from the run record), on the 1-degree grid as LaSRC
samples it (same sampling as scene_centre_values.py). Right panels: the continuous scene-centre value of each
scenario (dot) and the stored value after LaSRC's integer truncation (diamond); grey lines = storage steps
(0.1 kg/m2 PWAT, 2.5 DU ozone).

Usage: plot_fields_to_scene.py --record S1.json --v0 F1 F2 --s1 F1 F2 --s2 F1 F2 --s3 F1 F2 --name T15VVD --out DIR
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
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, Normalize  # noqa: E402

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
SEQ = LinearSegmentedColormap.from_list("seq", ["#d6e6f8", "#86b3e6", "#2a78d6", "#0d3f86"])
DIV = LinearSegmentedColormap.from_list("prgn", ["#5b2a86", "#a77fc8", "#efeeec", "#7fbf7b", "#1b7837"])
COL = {"V0": "#111111", "S1": "#2a78d6", "S2": "#eb6834", "S3": "#1baf7a"}
LAB = {"V0": "V0 GDAS analyses", "S1": "S1 GDAS f000+f006", "S2": "S2 GFS f000+f006", "S3": "S3 GFS f006+f012"}
plt.rcParams.update({"font.size": 6.5, "axes.titlesize": 6.6, "text.color": INK})


def cmg_field(g, lat_top, lat_bot, lon_left, lon_right):
    """Bilinear values on LaSRC's 0.05-degree CMG cells (same formula as lut_subr.c) covering the box."""
    from plot_spatial_interp import linesamp, bilinear
    ls = np.arange(int(np.ceil((90 - lat_top) / 0.05)), int(np.floor((90 - lat_bot) / 0.05)) + 1)
    ss = np.arange(int(np.ceil((lon_left + 180) / 0.05)), int(np.floor((lon_right + 180) / 0.05)) + 1)
    YY, XX = np.meshgrid(90 - ls * 0.05, -180 + ss * 0.05, indexing="ij")
    return bilinear(g, *linesamp(YY, XX))


def field(f1, f2, w1, w2, var, shift=1):
    A, B = S.read_fields(f1), S.read_fields(f2)
    g = A[var].astype(np.float32) * np.float32(w1) + B[var].astype(np.float32) * np.float32(w2)
    return np.roll(g, -(180 + shift), axis=1)            # columns: lon = j - 180, rows: lat = 90 - i


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--record", required=True)
    for s in ("v0", "s1", "s2", "s3"):
        ap.add_argument("--" + s, nargs=2, required=True)
    ap.add_argument("--name", default="")
    ap.add_argument("--half", type=int, default=4)
    ap.add_argument("--tile", default="", help="MGRS tile id (e.g. 15VVD) to outline in red")
    ap.add_argument("--landsat", action="store_true", help="outline an approximate WRS-2 scene (185 km) instead")
    ap.add_argument("--label", default="", help="text after the name in the title, e.g. 'S30, 27 Jan 2025'")
    ap.add_argument("--example", default="S1", help="scenario used in the worked example under the figure")
    ap.add_argument("--wv-range", nargs=2, type=float, default=[0, 20], help="V0 water vapour colour range (kg/m2)")
    ap.add_argument("--oz-range", nargs=2, type=float, default=[250, 470], help="V0 ozone colour range (DU)")
    ap.add_argument("--wv-pct", type=float, default=30, help="water vapour % difference colour range (+/-)")
    ap.add_argument("--oz-pct", type=float, default=3, help="ozone % difference colour range (+/-)")
    ap.add_argument("--grib", action="store_true", help="show the raw 1-degree GRIB cells instead of the CMG grid")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    r = json.load(open(a.record))
    w1, w2, lat, lon = r["first_weight"], r["second_weight"], r["scene_lat"], r["scene_lon"]
    files = {"V0": a.v0, "S1": a.s1, "S2": a.s2, "S3": a.s3}
    sc = {s: S.lasrc_scene_centre(f[0], f[1], w1, w2, lat, lon) for s, f in files.items()}
    ic, jc = int(round(90 - lat)), int(round(lon + 180))
    rows = slice(ic - a.half, ic + a.half + 1)
    cols = np.arange(jc - a.half, jc + a.half + 1) % 360
    lats = 90 - np.arange(ic - a.half, ic + a.half + 1)
    lons = np.arange(jc - a.half, jc + a.half + 1) - 180

    poly = None
    if a.tile:
        import mgrs
        from pyproj import Transformer
        z, h, e, n = mgrs.MGRS().MGRSToUTM(a.tile + "0000000000")
        epsg = (32600 if h == "N" else 32700) + int(z)
        tr = Transformer.from_crs(epsg, 4326, always_xy=True)
        cs = [tr.transform(e + de, n + dn) for de, dn in
              [(-40, 100020), (109760, 100020), (109760, -9780), (-40, -9780), (-40, 100020)]]
        poly = ([c[0] for c in cs], [c[1] for c in cs])
    box = "the granule's tile (~110 km)"
    if a.landsat:
        tilt = np.arcsin(np.cos(np.radians(98.2)) / np.cos(np.radians(lat)))   # descending track, top leans east
        cs = []
        for x, y in [(-92.5, -90), (92.5, -90), (92.5, 90), (-92.5, 90), (-92.5, -90)]:
            xe = x * np.cos(tilt) - y * np.sin(tilt); yn = x * np.sin(tilt) + y * np.cos(tilt)
            cs.append((lon + xe / (111.2 * np.cos(np.radians(lat))), lat + yn / 111.2))
        poly = ([c[0] for c in cs], [c[1] for c in cs])
        box = "approx. Landsat scene (~185 km)"
    fig = plt.figure(figsize=(7.3, 4.95), facecolor="white")
    gs = fig.add_gridspec(2, 6, width_ratios=[1, 1, 1, 1, 0.08, 1.65], wspace=0.12, hspace=0.62,
                          left=0.06, right=0.985, top=0.855, bottom=0.225)
    meta = {"PWAT": ("Water vapour (PWAT, kg/m²)", 10, "uwv", 0.1, "kg/m²", "%.2f"),
            "TOZNE": ("Total ozone (TOZNE, DU)", 0.4, "uoz", 2.5, "DU", "%.1f")}
    for r_, var in enumerate(("PWAT", "TOZNE")):
        label, mult, key, step, unit, fmt = meta[var]
        ext = [lons[0] - 0.5, lons[-1] + 0.5, lats[-1] - 0.5, lats[0] + 0.5]
        if a.grib:
            F = {s: field(f[0], f[1], w1, w2, var)[rows][:, cols] for s, f in files.items()}
        else:
            F = {s: cmg_field(field(f[0], f[1], w1, w2, var), ext[3], ext[2], ext[0], ext[1]) for s, f in files.items()}
        pct = {s: 100 * (F[s] - F["V0"]) / F["V0"] for s in ("S1", "S2", "S3")}
        pm = a.wv_pct if var == "PWAT" else a.oz_pct
        dn = TwoSlopeNorm(vcenter=0, vmin=-pm, vmax=pm)
        vr = a.wv_range if var == "PWAT" else a.oz_range
        an = Normalize(vmin=vr[0], vmax=vr[1])
        print(var, "V0 field %.2f-%.2f" % (F["V0"].min(), F["V0"].max()),
              "max |pct| %.2f" % max(float(np.abs(p).max()) for p in pct.values()))
        for c, s in enumerate(("V0", "S1", "S2", "S3")):
            ax = fig.add_subplot(gs[r_, c])
            if s == "V0":
                im0 = ax.imshow(F["V0"], cmap=SEQ, norm=an, extent=ext, aspect="auto", interpolation="nearest")
                ttl = "V0 field"
            else:
                im1 = ax.imshow(pct[s], cmap=DIV, norm=dn, extent=ext, aspect="auto", interpolation="nearest")
                ttl = "%s − V0" % s
            for k in (np.arange(lons[0] - 0.5, lons[-1] + 0.51, 1.0) if not a.grib else []):
                ax.axvline(k, color="white", lw=0.35, alpha=0.7, zorder=3)
            for k in (np.arange(lats[-1] - 0.5, lats[0] + 0.51, 1.0) if not a.grib else []):
                ax.axhline(k, color="white", lw=0.35, alpha=0.7, zorder=3)
            if poly:
                ax.plot(poly[0], poly[1], color="#d62728", lw=1.1, zorder=4)
                ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
            ax.plot(lon, lat, "o", ms=3.2, mfc="white", mec=INK, mew=0.8)
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_linewidth(0.4); sp.set_color(INK2)
            ax.set_title(ttl, fontsize=6.6, pad=2, color=COL[s] if s != "V0" else INK, fontweight="bold")
            v = sc[s][key]
            cont = sc[s]["pwat_interp" if var == "PWAT" else "tozne_interp"]
            if s == "V0":
                ax.text(0.5, 0.03, "centre " + (fmt % cont), transform=ax.transAxes, ha="center", va="bottom", fontsize=5.6,
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))
            else:
                d0 = 100 * (cont - sc["V0"]["pwat_interp" if var == "PWAT" else "tozne_interp"]) / \
                    sc["V0"]["pwat_interp" if var == "PWAT" else "tozne_interp"]
                d1 = 100 * (v - sc["V0"][key]) / sc["V0"][key]
                ax.text(0.5, 0.03, "centre %+.2f%%\nstored %+.1f%%" % (d0, d1), transform=ax.transAxes, ha="center",
                        va="bottom", fontsize=5.6, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))
            if c == 0:
                ax.set_ylabel(label.split(" (")[0], fontsize=6.8, fontweight="bold")
        # colourbars under the map columns
        cax = fig.add_axes([0.06, (0.225 if r_ else 0.6) - 0.045, 0.12, 0.012])
        cb = fig.colorbar(im0, cax=cax, orientation="horizontal", extend="both"); cb.ax.tick_params(labelsize=5, length=1.5)
        cb.outline.set_linewidth(0.3)
        cax = fig.add_axes([0.22, (0.225 if r_ else 0.6) - 0.045, 0.33, 0.012])
        cb = fig.colorbar(im1, cax=cax, orientation="horizontal", extend="both"); cb.ax.tick_params(labelsize=5, length=1.5)
        cb.outline.set_linewidth(0.3)
        cb.set_label("% difference from V0 (fixed scale per gas, same on all granules)", fontsize=5.3, labelpad=1)
        # right: continuous value -> stored value
        ax = fig.add_subplot(gs[r_, 5])
        conts = {s: sc[s]["pwat_interp" if var == "PWAT" else "tozne_interp"] for s in sc}
        stored = {s: (sc[s][key] * 200 / 10 if var == "PWAT" else sc[s][key] * 1000) for s in sc}
        lo = min(min(conts.values()), min(stored.values())) - step * 0.6
        hi = max(conts.values()) + step * 0.6
        k0, k1 = int(np.floor(lo / step)), int(np.ceil(hi / step))
        for k in range(k0, k1 + 1):
            ax.axvline(k * step, color=GRID, lw=0.8, zorder=0)
        for i, s in enumerate(("V0", "S1", "S2", "S3")):
            y = 3 - i
            ax.annotate("", xy=(stored[s], y), xytext=(conts[s], y),
                        arrowprops=dict(arrowstyle="->", color=COL[s], lw=0.8, shrinkA=2, shrinkB=3))
            ax.plot(conts[s], y, "o", ms=4, color=COL[s])
            ax.plot(stored[s], y, "D", ms=4.5, mfc="white", mec=COL[s], mew=1)
            ax.text(conts[s], y + 0.28, fmt % conts[s], ha="center", fontsize=5.4, color=COL[s])
        ax.set_yticks(range(4)); ax.set_yticklabels(["S3", "S2", "S1", "V0"], fontsize=6.2, fontweight="bold")
        for t, s in zip(ax.get_yticklabels(), ("S3", "S2", "S1", "V0")):
            t.set_color(COL[s])
        ax.set_ylim(-1.25, 3.75); ax.set_xlim(lo, hi)
        ax.tick_params(axis="x", labelsize=5.5, length=2); ax.tick_params(axis="y", length=0)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.set_xlabel("%s; grey lines = storage steps (%s %s)" % (unit, ("%g" % step), unit), fontsize=5.6,
                      labelpad=2)
        same = len(set(round(v, 4) for v in stored.values())) == 1
        sv = stored["V0"]
        if same:
            ax.text(sv, -0.95, ("all stored as %s (%s %s)" % (fmt % sv, key, ("%.2f" % sc["V0"][key]) if var == "PWAT" else ("%.4f" % sc["V0"][key]))),
                    ha="left", va="center", fontsize=5.6, color=INK, fontweight="bold")
        else:
            for i, s_ in enumerate(("V0", "S1", "S2", "S3")):
                ax.text(stored[s_], 3 - i - 0.33, (fmt % stored[s_]), ha="center", va="center", fontsize=5.2,
                        color=COL[s_])
        ax.set_title("scene centre ● → value LaSRC stores ◆", fontsize=6.4, fontweight="bold", pad=2)
    fig.text(0.06, 0.975, "%s (%s): the four sources give different fields; LaSRC reduces each to one scene-centre value and "
             "truncates it" % (a.name, a.label), fontsize=7.6, fontweight="bold", va="top")
    used = " · ".join("%s %.3f / %.4f" % (s_, sc[s_]["uwv"], sc[s_]["uoz"]) for s_ in ("V0", "S1", "S2", "S3"))
    if len(set((round(d["uwv"], 4), round(d["uoz"], 4)) for d in sc.values())) == 1:
        used = "water vapour %.3f, ozone %.4f in all four runs" % (sc["V0"]["uwv"], sc["V0"]["uoz"])
    else:
        used = "uwv / uoz: " + used
    lead = ("Each square = one original 1° GRIB cell (before LaSRC's bilinear step)" if a.grib else
            "Colours = LaSRC's 0.05° CMG cells after bilinear interpolation; white lines = 1° GRIB cells")
    fig.text(0.06, 0.94, lead + "; red box = %s; dot = scene centre.\nFields time-interpolated to the "
             "acquisition (%s UTC). Values LaSRC used: %s." % (box, r["acquisition"].replace("T", " "), used),
             fontsize=6, color=INK2, va="top")
    ex = a.example

    def wv_line(s_):
        p_ = sc[s_]["pwat_interp"]; n = int(p_ * 10)
        return "%s %.3f × 10.0 = %.2f → %d → %.2f kg/m²" % (s_, p_, p_ * 10, n, n / 10)

    def oz_line(s_):
        t_ = sc[s_]["tozne_interp"]; n = int(t_ * 0.4)
        return "%s %.2f × 0.40 = %.2f → %d → %.1f DU" % (s_, t_, t_ * 0.4, n, n * 2.5)

    def pc(s_, k):
        return 100 * (sc[s_][k] - sc["V0"][k]) / sc["V0"][k]
    fig.text(0.06, 0.0, "LaSRC stores wv = (uint16)(PWAT × 10.0) and oz = (uint8)(TOZNE × 0.40); the C cast drops the decimals "
             "(truncation, not rounding).\n"
             "Water vapour: %s;  %s  (stored %+.1f%%).\n"
             "Ozone: %s;  %s  (stored %+.1f%%)."
             % (wv_line("V0"), wv_line(ex), pc(ex, "uwv"), oz_line("V0"), oz_line(ex), pc(ex, "uoz")),
             fontsize=6.2, color=INK, va="bottom", linespacing=1.5)
    os.makedirs(a.out, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(a.out, "slide_fields_to_scene_%s%s.%s" % (a.name, "_grib1deg" if a.grib else "", ext)), dpi=300,
                    bbox_inches="tight", pad_inches=0.04)
    print({s: {k: round(v, 4) for k, v in d.items()} for s, d in sc.items()})


if __name__ == "__main__":
    main()
