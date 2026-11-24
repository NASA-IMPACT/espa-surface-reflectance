#!/usr/bin/env python3
"""How sensitive is LaSRC surface reflectance to the water vapour and ozone it is given?

Uses LaSRC's own gas-transmittance formulas and coefficients (lut_subr.c comptg;
compute_landsat_refl.c / compute_sentinel_refl.c):
    T_oz = exp(a_oz * m * uoz)                       m = 1/cos(SZA) + 1/cos(VZA)
    T_wv = exp(-a_wv * (m * uwv) ** b_wv)
LaSRC divides the top-of-atmosphere signal by these transmittances, so to first order a change in
uwv/uoz changes surface reflectance by  dSR/SR ~ -[d ln T_oz + d ln T_wv].

Figures (--out):
  sensitivity_curves_<P>.png       reflectance change vs % change in water vapour / ozone, per band
  predicted_vs_observed_<P>.png    each granule: change predicted from its uwv/uoz change vs change measured
Table (stdout + threshold_crossing_<P>.csv): input error needed to exceed each band's threshold.

Usage: gas_sensitivity.py --table scene_centre_values_l30.csv --change change_vs_v0_percent.csv
                          --v0-metrics <V0_vs_VIIRS metrics csv> --product L30 --out DIR
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_pilot_difference_map import key_of  # noqa: E402

# LaSRC coefficients for the HLS bands (index into LaSRC's band arrays)
COEF = {
    "L30": dict(bands=["B01", "B02", "B03", "B04", "B05", "B06", "B07"],
                names=["Ultrablue", "Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2"],
                oz=[-0.00255649, -0.0177861, -0.0969872, -0.0611428, 0.0001, 0.0001, 0.0001],
                wa=[2.29849e-27, 2.29849e-27, 0.00194772, 0.00404159, 0.000729136, 0.00067324, 0.0177533],
                wb=[0.999742, 0.999742, 0.775024, 0.774482, 0.893085, 0.939669, 0.65094]),
    # Sentinel-2 default build (bands 9 and 10 not processed): B1 B2 B3 B4 B5 B6 B7 B8 B8A B11 B12
    "S30": dict(bands=["B01", "B02", "B03", "B04", "B8A", "B11", "B12"],
                names=["Ultrablue", "Blue", "Green", "Red", "NIR narrow", "SWIR1", "SWIR2"],
                oz=[-0.00264691, -0.0272572, -0.0986512, -0.0500348, 0.0001, 0.0001, 0.0001],
                wa=[2.29849e-27, 2.29849e-27, 0.000777307, 0.00361051, 0.000390755, 0.000640155, 0.018006],
                wb=[0.999742, 0.999742, 0.891099, 0.754895, 0.900899, 0.943712, 0.647517]),
}
COL = ["#7b3294", "#2a78d6", "#1baf7a", "#d62728", "#8c564b", "#eb6834", "#1a1a1a"]
SCOL = {"S1": "#2a78d6", "S2": "#eb6834", "S3": "#1baf7a"}
INK2, EDGE = "#5a5955", "#cfccc4"


def ln_t(c, i, uwv, uoz, m):
    """ln(T_oz * T_wv) for band i (LaSRC formulas)."""
    return c["oz"][i] * m * uoz - c["wa"][i] * np.power(m * np.asarray(uwv, float), c["wb"][i])


def dsr_pct(c, i, uwv0, uoz0, uwv1, uoz1, m):
    return -(ln_t(c, i, uwv1, uoz1, m) - ln_t(c, i, uwv0, uoz0, m)) * 100


def m_of(sza_deg, vza_deg=0.0):
    return 1 / np.cos(np.radians(sza_deg)) + 1 / np.cos(np.radians(vza_deg))


def solar_zenith(lat, lon, when):
    """NOAA solar position (approximate, enough for the air-mass factor)."""
    doy = when.timetuple().tm_yday
    hour = when.hour + when.minute / 60 + when.second / 3600
    g = 2 * np.pi / 365 * (doy - 1 + (hour - 12) / 24)
    eqt = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g) - 0.014615 * np.cos(2 * g)
                    - 0.040849 * np.sin(2 * g))
    decl = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
            + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    tst = hour * 60 + eqt + 4 * lon
    ha = np.radians(tst / 4 - 180)
    la = np.radians(lat)
    return np.degrees(np.arccos(np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(ha)))


def curves(p, out, uwv0=1.0, uoz0=0.30):
    c = COEF[p]
    pct = np.linspace(-50, 50, 201)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0), facecolor="white")
    for i, (b, n) in enumerate(zip(c["bands"], c["names"])):
        for ax, var in zip(axes, ("wv", "oz")):
            lo, mid, hi = [], [], []
            for sza in (30, 45, 70):
                m = m_of(sza)
                u1 = uwv0 * (1 + pct / 100) if var == "wv" else uwv0
                o1 = uoz0 * (1 + pct / 100) if var == "oz" else uoz0
                y = dsr_pct(c, i, uwv0, uoz0, u1, o1, m)
                {30: lo, 45: mid, 70: hi}[sza].append(y)
            ax.fill_between(pct, lo[0], hi[0], color=COL[i], alpha=0.12, lw=0)
            ax.plot(pct, mid[0], color=COL[i], lw=1.8, label="%s %s" % (b, n))
    for ax, var, obs in ((axes[0], "water vapour", 7), (axes[1], "ozone", 2)):
        ax.axhline(0, color=EDGE, lw=0.8); ax.axvline(0, color=EDGE, lw=0.8)
        ax.axvspan(-obs, obs, color="#9a9893", alpha=0.18, lw=0)
        ax.text(0, ax.get_ylim()[1] * 0.92, "range of real\nforecast differences", ha="center", va="top",
                fontsize=7.5, color=INK2)
        ax.set_xlabel("change in the %s LaSRC uses (%%)" % var, fontsize=9)
        ax.set_title("Surface reflectance change for a given %s change" % var, loc="left", fontsize=10,
                     fontweight="bold")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("change in surface reflectance (%)", fontsize=9)
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.suptitle("%s: how much LaSRC's surface reflectance responds to water vapour and ozone, per band" % p,
                 x=0.01, ha="left", fontsize=11.5, fontweight="bold")
    fig.text(0.01, 0.905, "From LaSRC's own gas-absorption formulas and band coefficients. Line: sun 45° from "
             "overhead; shading: sun 30° to 70°. Reference values: water vapour %.1f g/cm² (LaSRC units), "
             "ozone %.0f DU." % (uwv0, uoz0 * 1000), ha="left", va="top", fontsize=8.5, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    f = os.path.join(out, "sensitivity_curves_%s.png" % p)
    fig.savefig(f, dpi=150); plt.close(fig); print("wrote", f)


def threshold_table(p, v0m, out, uwv0=1.0, uoz0=0.30, sza=45):
    """Input change needed for the predicted reflectance change to reach the band threshold."""
    c = COEF[p]; m = m_of(sza)
    rows = []
    for i, b in enumerate(c["bands"]):
        z = v0m[v0m["band"] == b]
        if z.empty:
            continue
        refl = float(z["mean_post_refl"].median())
        thr = float(z["threshold_refl"].iloc[0])
        need = 100 * thr / refl                                  # % change of reflectance = threshold
        res = {}
        for var in ("wv", "oz"):
            hit = None
            for k in np.arange(1, 1001):                         # search +1 ... +1000 %
                u1 = uwv0 * (1 + k / 100) if var == "wv" else uwv0
                o1 = uoz0 * (1 + k / 100) if var == "oz" else uoz0
                if abs(dsr_pct(c, i, uwv0, uoz0, u1, o1, m)) >= need:
                    hit = k; break
            res[var] = "%d%%" % hit if hit else "> 1000%"
        rows.append(dict(band=b, typical_reflectance=round(refl, 3), threshold_dn=round(thr * 1e4),
                         threshold_pct_of_reflectance=round(need, 1),
                         water_vapour_change_needed=res["wv"], ozone_change_needed=res["oz"]))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(out, "threshold_crossing_%s.csv" % p), index=False)
    print(t.to_string(index=False))
    return t


def predicted_vs_observed(p, table, change, out, stage):
    c = COEF[p]
    t = pd.read_csv(table)
    t = t[t["product"] == p].copy()
    t["wv"] = t["uwv_logged"].fillna(t["uwv_lasrc"]); t["oz"] = t["uoz_logged_du"].fillna(t["uoz_lasrc_du"]) / 1000
    t["key"] = t["granule"].map(key_of)
    v0 = {r["key"]: r for _, r in t[t["scenario"] == "V0"].iterrows()}
    ch = pd.read_csv(change, dtype={"tile_or_pathrow": str, "date": str})
    ch = ch[(ch["product"] == p) & (ch["stage"] == stage) & (ch["panel"] != "V0")]
    ch["key"] = list(zip(ch["tile_or_pathrow"], ch["date"]))
    obs = {(r["panel"], r["key"], r["band"]): r["pct"] for _, r in ch.iterrows()}
    rows = []
    for _, r in t[t["scenario"] != "V0"].iterrows():
        if r["key"] not in v0:
            continue
        b0 = v0[r["key"]]
        # Landsat and Sentinel-2 cross the equator at about 10:00-10:30 local solar time
        day = pd.Timestamp(r["key"][1])
        when = day + pd.Timedelta(hours=(10.25 - r["scene_lon"] / 15) % 24)
        sza = solar_zenith(r["scene_lat"], r["scene_lon"], when.to_pydatetime())
        m = m_of(min(sza, 80))
        for i, b in enumerate(c["bands"]):
            k = (r["scenario"], r["key"], b)
            if k in obs:
                rows.append(dict(scenario=r["scenario"], band=b, key=r["key"],
                                 pred=dsr_pct(c, i, b0["wv"], b0["oz"], r["wv"], r["oz"], m), obs=obs[k]))
    d = pd.DataFrame(rows)
    bands = [b for b in c["bands"] if b in set(d["band"])]
    ncol = 3; nrow = (len(bands) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.3 * ncol, 4.0 * nrow), facecolor="white", squeeze=False)
    stats = []
    for ax, b in zip(axes.flat, bands):
        z = d[d["band"] == b]
        lim = max(float(z[["pred", "obs"]].abs().quantile(0.995).max()) * 1.1, 0.01)
        ax.plot([-lim, lim], [-lim, lim], ls="--", color=INK2, lw=0.8)
        for s in ("S1", "S2", "S3"):
            y = z[z["scenario"] == s]
            ax.scatter(y["pred"], y["obs"], s=9, c=SCOL[s], alpha=0.6, lw=0, label=s)
        r = np.corrcoef(z["pred"], z["obs"])[0, 1] if len(z) > 2 and z["pred"].std() > 0 else np.nan
        within = float((np.abs(z["obs"] - z["pred"]) <= np.maximum(0.02, 0.5 * np.abs(z["pred"]))).mean() * 100)
        stats.append(dict(band=b, n=len(z), r=round(r, 2), pct_close=round(within)))
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_title("%s  |  r = %.2f" % (b, r), loc="left", fontsize=10, fontweight="bold")
        ax.set_xlabel("predicted from water vapour/ozone change (%)", fontsize=8)
        ax.set_ylabel("measured change vs V0 (%)", fontsize=8)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    for ax in axes.flat[len(bands):]:
        ax.axis("off")
    axes.flat[0].legend(frameon=False, fontsize=8, markerscale=2)
    fig.suptitle("%s: reflectance change predicted from LaSRC's formulas vs measured, every granule and scenario"
                 % p, x=0.01, ha="left", fontsize=11.5, fontweight="bold")
    fig.text(0.01, 0.935, "Dots on the dashed line: the measured change is exactly what the change in water vapour "
             "and ozone predicts. Measured = scenario minus V0, % of the granule's mean reflectance.",
             ha="left", va="top", fontsize=8.5, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    f = os.path.join(out, "predicted_vs_observed_%s.png" % p)
    fig.savefig(f, dpi=150); plt.close(fig); print("wrote", f)
    st = pd.DataFrame(stats); print(st.to_string(index=False))
    return d, st


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--product", default="L30")
    ap.add_argument("--table"); ap.add_argument("--change"); ap.add_argument("--v0-metrics")
    ap.add_argument("--stage", default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    plt.rcParams.update({"font.size": 9})
    curves(a.product, a.out)
    if a.v0_metrics:
        threshold_table(a.product, pd.read_csv(a.v0_metrics), a.out)
    if a.table and a.change:
        predicted_vs_observed(a.product, a.table, a.change, a.out, a.stage or ("immediate" if a.product == "L30" else "final"))


if __name__ == "__main__":
    main()
