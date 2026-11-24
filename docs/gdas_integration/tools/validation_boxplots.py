#!/usr/bin/env python3
"""Boxplots and summary table from the validation CLI metrics, for one product and stage.

Reads the newest lasrc_validation_metrics_*.csv in
  <root>/validation_full/<P>_<R>_vs_VIIRS_<stage>/reports/   R = V0, S1, S2, S3
  <root>/validation_full/<P>_<S>_vs_V0_<stage>/reports/      S = S1, S2, S3
and writes to --out:
  <P>_<stage>_boxplots_vs_VIIRS.png/.pdf   2 x 2 panels (V0, S1, S2, S3): per band, per-granule |mean difference|
                                            (HLS DN) with the band threshold as a red dashed line
  <P>_<stage>_boxplots_vs_V0.png/.pdf      1 x 3 panels (S1, S2, S3), same layout
  <P>_<stage>_summary.csv                  per comparison and band: granules, median / p95 / max |diff| (DN),
                                            median / max |diff| in % of the reference mean SR, fails, pass rate
  <P>_<stage>_fails.csv                    every FAIL row (granule, band, diff, threshold)
Only bands that have a threshold are plotted (e.g. S30 B08-B10 have none).

Usage:
  python validation_boxplots.py --product S30 --stage immediate --out $R/reports_full/validation_s30
  python validation_boxplots.py --product S30 --stage final     --out $R/reports_full/validation_s30
Needs pandas + matplotlib (the lasrc_validation venv has both).
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd

INK, INK2, GRID, RED = "#0b0b0b", "#52514e", "#e4e3df", "#c0392b"
LAB = {"V0": "V0 GDAS analyses", "S1": "S1 GDAS f000+f006", "S2": "S2 GFS f000+f006",
       "S3": "S3 GFS f006+f012 (12-h lead)"}
ORDER = ["B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B09", "B10", "B11", "B12"]
plt.rcParams.update({"font.size": 6.5, "text.color": INK})


def newest(folder):
    f = sorted(x for x in glob.glob(os.path.join(folder, "reports", "lasrc_validation_metrics_*.csv"))
               if not os.path.basename(x).startswith("._"))
    return f[-1] if f else None


def load(root, product, stage):
    out = {}
    for ref, runs in (("VIIRS", ("V0", "S1", "S2", "S3")), ("V0", ("S1", "S2", "S3"))):
        for r in runs:
            d = os.path.join(root, "validation_full", "%s_%s_vs_%s_%s" % (product, r, ref, stage))
            f = newest(d)
            if f is None:
                print("missing metrics:", d)
                continue
            m = pd.read_csv(f)
            m = m[m["threshold_refl"].notna() & (m["decision"].astype(str) != "NO_THRESHOLD")]
            out[(r, ref)] = m
            print("%-28s %s  granules=%d" % ("%s vs %s" % (r, ref), os.path.basename(f), m["granule"].nunique()))
    return out


def panels(data, keys, ncols, title, out, name):
    bands = [b for b in ORDER if any(b in set(data[k]["band"]) for k in keys)]
    nrows = (len(keys) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 2.3 * nrows + 0.5), facecolor="white",
                             sharey=True, squeeze=False)
    top = max(max(float((data[k]["abs_mean_diff_refl"] * 1e4).max()), float((data[k]["threshold_refl"] * 1e4).max()))
              for k in keys)
    for ax, k in zip(axes.flat, keys):
        m = data[k]
        bs = [b for b in bands if b in set(m["band"])]
        vals = [m[m["band"] == b]["abs_mean_diff_refl"].dropna().to_numpy() * 1e4 for b in bs]
        thr = [float(m[m["band"] == b]["threshold_refl"].iloc[0]) * 1e4 for b in bs]
        lbl = [str(m[m["band"] == b]["band_label"].iloc[0]) if "band_label" in m else "" for b in bs]
        bp = ax.boxplot(vals, patch_artist=True, widths=0.55, showfliers=True,
                        flierprops=dict(marker="o", markersize=2.2, markerfacecolor="none", markeredgecolor=INK2,
                                        markeredgewidth=0.4),
                        medianprops=dict(color="#0d3f86", lw=1.1), whiskerprops=dict(lw=0.6),
                        capprops=dict(lw=0.6), boxprops=dict(lw=0.6))
        for p in bp["boxes"]:
            p.set_facecolor("#cfe0f5")
        for i, t in enumerate(thr, start=1):
            ax.hlines(t, i - 0.32, i + 0.32, colors=RED, linestyles="--", lw=1.1, zorder=3)
            ax.text(i, t + top * 0.015, "%.0f" % t, ha="center", va="bottom", fontsize=5.6, color=RED)
        nfail = int(m.groupby("granule")["decision"].apply(
            lambda d: d.astype(str).str.upper().str.startswith("FAIL").any()).sum())
        ax.set_title("%s vs %s  |  %d granules, %d fail" % (LAB[k[0]], k[1], m["granule"].nunique(), nfail),
                     loc="left", fontsize=6.8, fontweight="bold")
        ax.set_xticks(range(1, len(bs) + 1))
        ax.set_xticklabels(["%s\n%s" % (b, l) for b, l in zip(bs, lbl)], fontsize=5.8)
        ax.set_ylim(-top * 0.03, top * 1.12)
        ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    for ax in axes.flat[len(keys):]:
        ax.axis("off")
    for ax in axes[:, 0]:
        ax.set_ylabel("|mean difference| (HLS DN)", fontsize=6.6)
    fig.suptitle(title, x=0.01, ha="left", fontsize=7.6, fontweight="bold")
    fig.legend(handles=[Line2D([0], [0], color=RED, ls="--", lw=1.1, label="band threshold (Ju et al. 2025)")],
               loc="lower center", frameon=False, fontsize=6.6)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95), h_pad=0.8, w_pad=0.6)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, "%s.%s" % (name, ext)), dpi=300, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/rtmp/tvo/LaSRC_ancillary_latency_20261120")
    ap.add_argument("--product", default="S30", choices=["S30", "L30"])
    ap.add_argument("--stage", default="immediate", choices=["immediate", "final"])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    data = load(a.root, a.product, a.stage)
    tag = "%s_%s" % (a.product, a.stage)
    vk = [k for k in (("V0", "VIIRS"), ("S1", "VIIRS"), ("S2", "VIIRS"), ("S3", "VIIRS")) if k in data]
    ok = [k for k in (("S1", "V0"), ("S2", "V0"), ("S3", "V0")) if k in data]
    if vk:
        panels(data, vk, 2, "%s %s: each run vs the VIIRS-based run, per-granule |mean difference| by band"
               % (a.product, a.stage), a.out, tag + "_boxplots_vs_VIIRS")
    if ok:
        panels(data, ok, 3, "%s %s: each forecast scenario vs V0 (GDAS analyses), per-granule |mean difference| by band"
               % (a.product, a.stage), a.out, tag + "_boxplots_vs_V0")
    rows, fails = [], []
    for (r, ref), m in data.items():
        m = m.assign(dn=m["abs_mean_diff_refl"] * 1e4,
                     pct=100 * m["abs_mean_diff_refl"] / m["mean_pre_refl"].where(m["mean_pre_refl"] > 0))
        for b, g in m.groupby("band"):
            f = g["decision"].astype(str).str.upper().str.startswith("FAIL")
            rows.append(dict(comparison="%s vs %s" % (r, ref), band=b, granules=g["granule"].nunique(),
                             median_dn=round(g["dn"].median(), 2), p95_dn=round(g["dn"].quantile(0.95), 2),
                             max_dn=round(g["dn"].max(), 2), median_pct=round(g["pct"].median(), 3),
                             max_pct=round(g["pct"].max(), 3), threshold_dn=round(g["threshold_refl"].iloc[0] * 1e4, 1),
                             fails=int(f.sum()), pass_rate_pct=round(100 * (1 - f.mean()), 2)))
            for _, x in g[f].iterrows():
                fails.append(dict(comparison="%s vs %s" % (r, ref), granule=x["granule"], band=b,
                                  mean_diff_dn=round(x["mean_diff_refl"] * 1e4, 1),
                                  mean_ref_refl=round(x["mean_pre_refl"], 4),
                                  threshold_dn=round(x["threshold_refl"] * 1e4, 1)))
    s = pd.DataFrame(rows)
    s["band"] = pd.Categorical(s["band"], [b for b in ORDER if b in set(s["band"])])
    s = s.sort_values(["comparison", "band"])
    s.to_csv(os.path.join(a.out, tag + "_summary.csv"), index=False)
    pd.DataFrame(fails, columns=["comparison", "granule", "band", "mean_diff_dn", "mean_ref_refl", "threshold_dn"]) \
        .to_csv(os.path.join(a.out, tag + "_fails.csv"), index=False)
    print("\nfailing granules per comparison:")
    for k, m in data.items():
        fg = sorted(set(m[m["decision"].astype(str).str.upper().str.startswith("FAIL")]["granule"]))
        print("  %-14s %d  %s" % ("%s vs %s" % k, len(fg), ", ".join(fg)))
    print("\nwrote", a.out)


if __name__ == "__main__":
    main()
