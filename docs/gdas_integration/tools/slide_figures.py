#!/usr/bin/env python3
"""Slide-sized versions of the L30 figures used in the team deck.

Each figure is drawn at the size it occupies on a 10 x 5.625 in slide (about 7.3 x 4.45 in) with fonts sized for
that, and saved at 300 dpi, so text stays sharp when shown at full size (the analysis versions are drawn ~13 in
wide and their text shrinks to ~5 pt on a slide).

Inputs:
  --table     reports_full/scene_centre_values_l30.csv
  --vs-viirs  validation_full/L30_<R>_vs_VIIRS_immediate/reports/lasrc_validation_metrics_*.csv  (R = V0,S1,S2,S3)
  --vs-v0     validation_full/L30_<S>_vs_V0_immediate/reports/lasrc_validation_metrics_*.csv     (S = S1,S2,S3)
Usage:
  slide_figures.py --table T --vs-viirs V0.csv --vs-viirs S1.csv ... --vs-v0 S1.csv ... --out DIR
Writes: slide_vs_viirs_B03.png, slide_change_vs_v0_B03.png, slide_uwv_map.png, slide_uoz_map.png,
        slide_boxplots.png, slide_sr_vs_inputs.png
"""
import argparse
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plot_pilot_difference_map as M  # noqa: E402  palettes, land outlines
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import Normalize, TwoSlopeNorm, to_rgba  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

W, H, DPI = 7.3, 4.45, 300
INK, INK2, GRID, RED = "#0b0b0b", "#52514e", "#e4e3df", "#d03b3b"
SC = {"S1": "#2a78d6", "S2": "#eb6834", "S3": "#1baf7a"}
LAB = {"V0": "V0  GDAS analyses (current)", "S1": "S1  GDAS f000 + f006", "S2": "S2  GFS f000 + f006",
       "S3": "S3  GFS f006 + f012 (12-h lead)"}
SMIN, SMAX = 3, 26
plt.rcParams.update({"font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5,
                     "ytick.labelsize": 6.5, "text.color": INK, "axes.titlecolor": INK})



def save(fig, out, name):
    """PNG at 300 dpi plus a vector PDF (never blurry at any zoom)."""
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, "%s.%s" % (name, ext)), dpi=DPI, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


# ------------------------------------------------------------------ helpers
def map_figure(title, nrow=2):
    """2 x nrow world panels filling a slide-sized figure; returns fig, axes, colourbar y."""
    fig = plt.figure(figsize=(W, H), facecolor="white")
    pw, gx, lm = 3.5, 0.12, 0.09
    ph = pw * 134 / 360
    top, ttl = 0.32, 0.17
    axes = []
    for r in range(nrow):
        y0 = H - top - ttl - r * (ph + ttl + 0.1) - ph
        for c in range(2):
            axes.append(fig.add_axes([(lm + c * (pw + gx)) / W, y0 / H, pw / W, ph / H]))
    fig.text(lm / W, 1 - 0.08 / H, title, ha="left", va="top", fontsize=8, fontweight="bold")
    cb_y = (H - top - nrow * (ph + ttl + 0.1) - 0.12) / H
    return fig, axes, cb_y


def land(ax, rings):
    ax.add_collection(PolyCollection(rings, facecolor=M.LAND, edgecolor=M.EDGE, linewidth=0.2, zorder=0))
    ax.set_xlim(-180, 180); ax.set_ylim(-56, 78); ax.set_aspect("equal")
    ax.set_facecolor(M.SEA); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color(M.EDGE); sp.set_linewidth(0.5)


def cbar(fig, norm, cmap, x0, x1, y, label, ticks=None):
    cax = fig.add_axes([x0, y, x1 - x0, 0.022])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.ax.tick_params(labelsize=6, length=2, width=0.5)
    cb.outline.set_linewidth(0.5)
    if ticks is not None:
        cb.set_ticks(ticks)
    cb.set_label(label, fontsize=6.3)
    return cb


def size_of(v, vmax):
    return SMIN + (SMAX - SMIN) * (np.abs(v).clip(0, vmax) / vmax) ** 0.5


def load(a):
    t = pd.read_csv(a.table, dtype={"granule": str})
    t = t[t["product"] == "L30"]
    loc = t.groupby("granule")[["scene_lat", "scene_lon"]].first()
    used = {}
    for var, rep, log in (("uwv", "uwv_lasrc", "uwv_logged"), ("uoz", "uoz_lasrc_du", "uoz_logged_du")):
        v = t[log].where(t["scenario"] != "V0").fillna(t[rep])
        used[var] = t.assign(v=v).pivot_table(index="granule", columns="scenario", values="v", aggfunc="first")
    viirs = {re.search(r"L30_(\w+?)_vs_VIIRS", f).group(1): pd.read_csv(f) for f in a.vs_viirs}
    v0 = {re.search(r"L30_(S\d)_vs_V0", f).group(1): pd.read_csv(f) for f in a.vs_v0}
    return loc, used, viirs, v0


def pct_inputs(used, s):
    return {var: 100 * (used[var][s] - used[var]["V0"]) / used[var]["V0"] for var in used}



def pie(ax, n_pass, n_fail):
    """Small pass/fail pie in the South-Pacific corner of a map panel (same as the analysis maps)."""
    n = n_pass + n_fail
    ia = ax.inset_axes([0.012, 0.2, 0.11, 0.27])
    ia.pie([n_pass, n_fail], colors=[M.PASS, RED], startangle=90, counterclock=False,
           wedgeprops=dict(edgecolor="white", linewidth=0.6))
    ia.set_aspect("equal")
    pct = lambda k: ("%.1f" % (100 * k / n)) if 0 < 100 * k / n < 1 or 99 < 100 * k / n < 100 else "%.0f" % (100 * k / n)
    ia.text(0, -1.25, "pass %d (%s%%)\nfail %d (%s%%)" % (n_pass, pct(n_pass), n_fail, pct(n_fail)),
            ha="center", va="top", fontsize=5, color=INK, linespacing=1.2)


# ------------------------------------------------------------------ figures
def fig_vs_viirs(loc, viirs, rings, out, band="B03"):
    vmax, thr = 40.0, None
    fig, axes, cb_y = map_figure("L30 B03 Green: each run compared with the VIIRS-based run (signed mean difference, HLS DN)")
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
    for ax, r in zip(axes, ("V0", "S1", "S2", "S3")):
        land(ax, rings)
        m = viirs[r][viirs[r]["band"] == band].set_index("granule")
        thr = float(m["threshold_refl"].iloc[0]) * 1e4
        dn = m["mean_diff_refl"] * 1e4
        fail = m["decision"].astype(str).str.upper().str.startswith("FAIL")
        ll = loc.loc[dn.index]
        o = np.argsort(np.abs(dn.to_numpy()))
        ax.scatter(ll["scene_lon"].to_numpy()[o], ll["scene_lat"].to_numpy()[o], s=size_of(dn.to_numpy()[o], vmax),
                   c=dn.clip(-vmax, vmax).to_numpy()[o], cmap=M.BRTE, norm=norm,
                   edgecolors=np.where(fail.to_numpy()[o], RED, INK2), linewidths=np.where(fail.to_numpy()[o], 1.1, 0.3),
                   zorder=3)
        ax.set_title("%s  |  median |diff| %.1f DN  |  fail %d of %d" % (LAB[r], dn.abs().median(), fail.sum(), len(dn)),
                     loc="left", fontsize=6.6, pad=2)
        pie(ax, int((~fail).sum()), int(fail.sum()))
    cb = cbar(fig, norm, M.BRTE, 0.2, 0.8, cb_y,
              "difference from VIIRS run (HLS DN): teal = brighter, brown = darker; red outline = above the %.0f DN "
              "threshold" % thr)
    save(fig, out, "slide_vs_viirs_B03")


def fig_change_v0(loc, used, v0, rings, out, band="B03"):
    pan = {}
    for s, m in v0.items():
        m = m[m["band"] == band].set_index("granule")
        pin = pct_inputs(used, s)
        pan[s] = (100 * m["mean_diff_refl"] / m["mean_pre_refl"],
                  (pin["uwv"].abs() < 1e-9) & (pin["uoz"].abs() < 1e-9),
                  m["decision"].astype(str).str.upper().str.startswith("FAIL"))
    vmax = max(float(p.abs().max()) for p, _, _ in pan.values())
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
    fig, axes, cb_y = map_figure("L30 B03 Green: change in surface reflectance vs V0 (% of V0 reflectance)")
    for ax, s in zip(axes, ("S1", "S2", "S3")):
        land(ax, rings)
        p, same, fl = pan[s]
        same = same.reindex(p.index).fillna(False)
        ll = loc.loc[p.index]
        ax.scatter(ll.loc[same, "scene_lon"], ll.loc[same, "scene_lat"], s=6, facecolors="none", edgecolors="#9a9893",
                   linewidths=0.5, zorder=3)
        d = p[~same]
        o = np.argsort(np.abs(d.to_numpy()))
        ax.scatter(ll.loc[d.index, "scene_lon"].to_numpy()[o], ll.loc[d.index, "scene_lat"].to_numpy()[o],
                   s=size_of(d.to_numpy()[o], vmax), c=d.to_numpy()[o], cmap=M.BRTE, norm=norm, edgecolors=INK2,
                   linewidths=0.3, zorder=4)
        ax.set_title("%s vs V0  |  median %.2f%%, largest %.2f%%" % (LAB[s], p.abs().median(), p.abs().max()),
                     loc="left", fontsize=6.6, pad=2)
    ax = axes[3]
    ax.axis("off")
    rows = [["|change| above", "S1", "S2", "S3"]]
    for t_ in (0.10, 0.25, 0.50):
        rows.append(["%.2f%%" % t_] + [str(int((pan[s][0].abs() > t_).sum())) for s in ("S1", "S2", "S3")])
    tb = ax.table(cellText=rows, loc="center", cellLoc="center", colWidths=[0.34, 0.16, 0.16, 0.16],
                  bbox=[0.12, 0.38, 0.76, 0.55])
    tb.auto_set_font_size(False); tb.set_fontsize(6.8)
    for (i, j), c in tb.get_celld().items():
        c.set_edgecolor(GRID); c.set_linewidth(0.5)
        if i == 0:
            c.set_text_props(fontweight="bold")
    ax.text(0.12, 0.2, "granules out of 614.  B03 threshold (54 DN) is about\n4.7% of a typical green reflectance; "
            "hollow grey = same\nwater vapour and ozone as V0", transform=ax.transAxes, fontsize=6.3, color=INK2,
            va="top")
    cbar(fig, norm, M.BRTE, 0.2, 0.8, cb_y, "change vs V0 (%): teal = brighter, brown = darker; bigger dot = bigger change")
    save(fig, out, "slide_change_vs_v0_B03")


def fig_inputs(loc, used, rings, out):
    pmax = max(float(pct_inputs(used, s)[v].abs().max()) for s in ("S1", "S2", "S3") for v in used)
    meta = {"uwv": ("water vapour", "g/cm²", "%.2f"), "uoz": ("total ozone", "DU", "%.0f")}
    for var, (label, unit, fmt) in meta.items():
        w = used[var]
        fig, axes, cb_y = map_figure("L30: %s LaSRC used (after its rounding). First panel: V0 value; others: change "
                                     "vs V0 in %%" % label)
        v0 = w["V0"].dropna()
        anorm = Normalize(vmin=float(v0.min()), vmax=float(v0.max()))
        dnorm = TwoSlopeNorm(vcenter=0, vmin=-pmax, vmax=pmax)
        for ax, s in zip(axes, ("V0", "S1", "S2", "S3")):
            land(ax, rings)
            if s == "V0":
                ll = loc.loc[v0.index]
                ax.scatter(ll["scene_lon"], ll["scene_lat"], s=9, c=v0, cmap=M.SEQ, norm=anorm, edgecolors=INK2,
                           linewidths=0.3, zorder=3)
                ax.set_title("V0  GDAS analyses  |  %s to %s %s" % (fmt % v0.min(), fmt % v0.max(), unit), loc="left",
                             fontsize=6.6, pad=2)
                continue
            p = (100 * (w[s] - w["V0"]) / w["V0"]).dropna()
            same, diff = p[p.abs() < 1e-9], p[p.abs() >= 1e-9]
            ax.scatter(loc.loc[same.index, "scene_lon"], loc.loc[same.index, "scene_lat"], s=6, facecolors="none",
                       edgecolors="#9a9893", linewidths=0.5, zorder=3)
            o = np.argsort(np.abs(diff.to_numpy()))
            ax.scatter(loc.loc[diff.index, "scene_lon"].to_numpy()[o], loc.loc[diff.index, "scene_lat"].to_numpy()[o],
                       s=size_of(diff.to_numpy()[o], pmax), c=diff.to_numpy()[o], cmap=M.DIV2, norm=dnorm,
                       edgecolors=INK2, linewidths=0.3, zorder=4)
            ax.set_title("%s  |  %d of %d same as V0, largest %.1f%%" % (LAB[s], len(same), len(p), p.abs().max()),
                         loc="left", fontsize=6.6, pad=2)
        cbar(fig, anorm, M.SEQ, 0.03, 0.45, cb_y, "V0 %s (%s)" % (label, unit))
        cbar(fig, dnorm, M.DIV2, 0.53, 0.97, cb_y,
             "change vs V0 (%), same scale for water vapour and ozone; hollow grey = same as V0")
        save(fig, out, "slide_%s_map" % var)


def fig_boxplots(viirs, out):
    bands = ["B01", "B02", "B03", "B04", "B06", "B07"]
    fig, axes = plt.subplots(2, 2, figsize=(W, H + 0.05), facecolor="white", sharey=True)
    top = 0
    for r in ("V0", "S1", "S2", "S3"):
        m = viirs[r]
        top = max(top, float((m[m["band"].isin(bands)]["abs_mean_diff_refl"] * 1e4).max()),
                  float((m[m["band"].isin(bands)]["threshold_refl"] * 1e4).max()))
    for ax, r in zip(axes.flat, ("V0", "S1", "S2", "S3")):
        m = viirs[r]
        data = [m[m["band"] == b]["abs_mean_diff_refl"].dropna().to_numpy() * 1e4 for b in bands]
        thr = [float(m[m["band"] == b]["threshold_refl"].iloc[0]) * 1e4 for b in bands]
        lbl = [str(m[m["band"] == b]["band_label"].iloc[0]) for b in bands]
        bp = ax.boxplot(data, patch_artist=True, widths=0.55, showfliers=True,
                        flierprops=dict(marker="o", markersize=2.2, markerfacecolor="none", markeredgecolor=INK2,
                                        markeredgewidth=0.4),
                        medianprops=dict(color="#0d3f86", lw=1.1), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6),
                        boxprops=dict(lw=0.6))
        for pch in bp["boxes"]:
            pch.set_facecolor("#cfe0f5")
        for i, t_ in enumerate(thr, start=1):
            ax.hlines(t_, i - 0.32, i + 0.32, colors=RED, linestyles="--", lw=1.1, zorder=3)
            ax.text(i, t_ + 2, "%.0f" % t_, ha="center", va="bottom", fontsize=5.8, color=RED)
        nfail = int(m[m["band"].isin(bands)].groupby("granule")["decision"].apply(
            lambda d: d.astype(str).str.upper().str.startswith("FAIL").any()).sum())
        ax.set_title("%s vs VIIRS  |  %d granules, %d fail" % (LAB[r], m["granule"].nunique(), nfail), loc="left",
                     fontsize=6.8, fontweight="bold")
        ax.set_xticks(range(1, len(bands) + 1))
        ax.set_xticklabels(["%s\n%s" % (b, l) for b, l in zip(bands, lbl)], fontsize=6)
        ax.set_ylim(-3, top * 1.12)
        ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel("|mean difference| (HLS DN)", fontsize=6.6)
    fig.legend(handles=[Line2D([0], [0], color=RED, ls="--", lw=1.1, label="band threshold (Ju et al. 2025)")],
               loc="lower center", frameon=False, fontsize=6.6)
    fig.tight_layout(rect=(0, 0.04, 1, 1), h_pad=0.6, w_pad=0.6)
    save(fig, out, "slide_boxplots")


def fig_scatter(used, v0, out):
    rows = []
    for s, m in v0.items():
        m = m.assign(dn=100 * m["mean_diff_refl"] / m["mean_pre_refl"])[["granule", "band", "dn"]].copy()  # % of V0 mean SR
        pin = pct_inputs(used, s)
        for var in pin:
            m["p_" + var] = m["granule"].map(pin[var])
        m["same"] = (m["p_uwv"].abs() < 1e-9) & (m["p_uoz"].abs() < 1e-9)
        rows.append(m.assign(scenario=s))
    d = pd.concat(rows, ignore_index=True)
    fig, axes = plt.subplots(2, 3, figsize=(W, H + 0.05), facecolor="white", sharex=True, sharey=True)
    bands = [("B03", "Green"), ("B04", "Red"), ("B07", "SWIR2")]
    names = {"uwv": "water vapour", "uoz": "ozone"}
    for r, var in enumerate(("uwv", "uoz")):
        for c, (b, nm) in enumerate(bands):
            ax = axes[r, c]
            x = d[d["band"] == b]
            y = x[x["same"]]
            ax.scatter(y["p_" + var], y["dn"], s=12, facecolor="none", edgecolor="#6b6a66", linewidth=0.5, zorder=2)
            for s, col in SC.items():
                y = x[(x["scenario"] == s) & ~x["same"]]
                ax.scatter(y["p_" + var], y["dn"], s=12, facecolor=to_rgba(col, 0.4), edgecolor=col, linewidth=0.5,
                           zorder=3)
            ax.axhline(0, color="#c9c7c0", lw=0.6, zorder=0); ax.axvline(0, color="#c9c7c0", lw=0.6, zorder=0)
            ax.set_title("%s %s" % (b, nm), loc="left", fontsize=7, fontweight="bold", pad=2)
            ax.set_xlabel("%s change vs V0 (%%)" % names[var], fontsize=6.5)
            if c == 0:
                ax.set_ylabel("SR change vs V0 (%)", fontsize=6.5)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            for sp in ("left", "bottom"):
                sp_ = ax.spines[sp]; sp_.set_linewidth(0.6)
            ax.tick_params(labelbottom=True, labelleft=True, width=0.5, length=2)
    h = [Line2D([0], [0], marker="o", ls="none", markersize=4.5, markerfacecolor=to_rgba(col, 0.4),
                markeredgecolor=col, label={"S1": "S1 GDAS f000+f006", "S2": "S2 GFS f000+f006",
                                            "S3": "S3 GFS f006+f012"}[s]) for s, col in SC.items()]
    h.append(Line2D([0], [0], marker="o", ls="none", markersize=4.5, markerfacecolor="none", markeredgecolor="#6b6a66",
                    label="same water vapour and ozone as V0"))
    fig.legend(handles=h, loc="lower center", ncol=4, frameon=False, fontsize=6.5)
    fig.tight_layout(rect=(0, 0.05, 1, 1), h_pad=0.5, w_pad=0.5)
    save(fig, out, "slide_sr_vs_inputs")


def fig_wait(out):
    """Slide 3: how long HLS waits for ancillary data after the overpass, per option (hours)."""
    rows = [("V0  GDAS analyses (current)", 7, 13, "7–13 h", "#2a78d6"),
            ("S1  GDAS forecast", 1, 7, "1–7 h", "#1e7b45"),
            ("S2  GFS forecast", 0, 3.6, "0–3.6 h", "#1e7b45"),
            ("S3  GFS, 12-h lead", 0, 0, "0 h: files exist before the overpass", "#1e7b45")]
    fig, ax = plt.subplots(figsize=(5.3, 2.5), facecolor="white")
    for i, (lab, lo, hi, txt, col) in enumerate(rows):
        y = len(rows) - 1 - i
        if hi > lo:
            ax.barh(y, hi - lo, left=lo, height=0.55, color=col)
            ax.text(lo + 0.15, y, txt, va="center", ha="left", fontsize=7.5, fontweight="bold",
                    color="white" if hi - lo >= 3 else INK)
        else:
            ax.plot([0], [y], marker="|", markersize=10, color=col, mew=2.5)
            ax.text(0.25, y, txt, va="center", ha="left", fontsize=7.5, fontweight="bold", color=INK)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=7.5)
    ax.set_xlim(0, 14); ax.set_xticks(range(0, 15, 2))
    ax.set_xlabel("hours after the overpass", fontsize=7)
    ax.xaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_title("Wait for ancillary data after the overpass", loc="left", fontsize=8.5, fontweight="bold",
                 color="#0b5394")
    ax.text(0, -0.42, "VIIRS today: 48–72 h (off this scale)", transform=ax.transAxes, fontsize=6.5, color=INK2)
    fig.tight_layout()
    save(fig, out, "slide_wait_times")


def fig_storage(table, used, out):
    """Before vs after LaSRC's integer storage: % difference from V0 of the continuous value (x) and of the stored
    value LaSRC uses (y). Dots on y = 0 (hollow) had a real input difference that the storage step erased."""
    t = pd.read_csv(table, dtype={"granule": str})
    t = t[t["product"] == "L30"]
    cont = {"uwv": t.pivot_table(index="granule", columns="scenario", values="pwat_continuous_kgm2"),
            "uoz": t.pivot_table(index="granule", columns="scenario", values="tozne_continuous_du")}
    names = {"uwv": ("Water vapour", "0.1 kg/m² steps"), "uoz": ("Total ozone", "2.5 DU steps")}
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.25), facecolor="white")
    for ax, var in zip(axes, ("uwv", "uoz")):
        lim = 0
        for sc, col in SC.items():
            x = 100 * (cont[var][sc] - cont[var]["V0"]) / cont[var]["V0"]
            y = 100 * (used[var][sc] - used[var]["V0"]) / used[var]["V0"]
            k = x.index.intersection(y.index)
            x, y = x[k], y[k]
            lim = max(lim, float(x.abs().max()), float(y.abs().max()))
            z = y.abs() < 1e-9
            ax.scatter(x[z], y[z], s=10, facecolor="none", edgecolor="#6b6a66", linewidth=0.5, zorder=2)
            ax.scatter(x[~z], y[~z], s=10, facecolor=to_rgba(col, 0.4), edgecolor=col, linewidth=0.5, zorder=3)
        lim *= 1.05
        ax.plot([-lim, lim], [-lim, lim], ls="--", color=INK2, lw=0.6, zorder=1)
        ax.axhline(0, color="#c9c7c0", lw=0.6, zorder=0); ax.axvline(0, color="#c9c7c0", lw=0.6, zorder=0)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
        ax.set_title("%s (LaSRC stores it in %s)" % names[var], loc="left", fontsize=7, fontweight="bold", pad=3)
        ax.set_xlabel("difference from V0 before storage (%)", fontsize=6.5)
        ax.set_ylabel("difference from V0 in the value LaSRC uses (%)", fontsize=6.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(width=0.5, length=2)
    h = [Line2D([0], [0], marker="o", ls="none", markersize=4.5, markerfacecolor=to_rgba(col, 0.4),
                markeredgecolor=col, label=sc) for sc, col in SC.items()]
    h.append(Line2D([0], [0], marker="o", ls="none", markersize=4.5, markerfacecolor="none",
                    markeredgecolor="#6b6a66", label="difference erased by storage"))
    h.append(Line2D([0], [0], ls="--", color=INK2, lw=0.6, label="no change by storage"))
    fig.legend(handles=h, loc="lower center", ncol=5, frameon=False, fontsize=6.3, handletextpad=0.3,
               columnspacing=1.0)
    fig.tight_layout(rect=(0, 0.09, 1, 1), w_pad=1.2)
    save(fig, out, "slide_storage_effect")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--vs-viirs", action="append", required=True)
    ap.add_argument("--vs-v0", action="append", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rings = json.load(open(os.path.join(HERE, "world_land_110m.json")))["rings"]
    loc, used, viirs, v0 = load(a)
    fig_vs_viirs(loc, viirs, rings, a.out)
    fig_change_v0(loc, used, v0, rings, a.out)
    fig_inputs(loc, used, rings, a.out)
    fig_boxplots(viirs, a.out)
    fig_scatter(used, v0, a.out)
    fig_wait(a.out)
    fig_storage(a.table, used, a.out)
    print("wrote 8 figures to", a.out)


if __name__ == "__main__":
    main()
