#!/usr/bin/env python3
"""Show which ancillary (GRIB) files each scenario actually used, from the per-run records, and check every run
against the selection rule in select_pair.py.

For every run record (records/<product>/<scenario>/<granule>.json) the two staged files are decoded from their
source path into: model (GDAS/GFS), issue cycle, forecast hour, valid time. These are expressed in hours relative
to C1 (the 6-hourly cycle at/before the acquisition) and compared with the rule:

  slot  V0 (current)        S1                    S2                   S3
  1     GDAS C1 f000        GDAS C1 f000          GFS C1 f000          GFS CP f006  (valid C1)
  2     GDAS C2 f000        GDAS C1 f006          GFS C1 f006          GFS CP f012  (valid C2)
  (C2 = C1 + 6 h, CP = C1 - 6 h)

Steps (collect is stdlib-only, so it runs on compute nodes; plot needs matplotlib):
  plot_scenario_inputs.py collect --root R --products "l30 s30" --scenarios "V0 S1 S2 S3" --out inputs.csv
  plot_scenario_inputs.py plot --csv inputs.csv [--csv more.csv ...] --out DIR [--granule LC08_L1TP_001057]
Writes DIR/scenario_inputs_timeline.png (one granule, all scenarios) and DIR/scenario_inputs_check.png (all runs),
and prints a table of runs that follow / break the rule.
"""
import argparse
import csv
import glob
import json
import os
import re
import sys
from datetime import datetime, timedelta

RULE = {  # scenario -> [(model, cycle offset h, forecast h)] for slot 1, slot 2   (offsets relative to C1)
    "V0": [("GDAS", 0, 0), ("GDAS", 6, 0)],
    "S1": [("GDAS", 0, 0), ("GDAS", 0, 6)],
    "S2": [("GFS", 0, 0), ("GFS", 0, 6)],
    "S3": [("GFS", -6, 6), ("GFS", -6, 12)],
}
LABEL = {"V0": "V0  GDAS analyses (current)", "S1": "S1  GDAS f000 + f006", "S2": "S2  GFS f000 + f006",
         "S3": "S3  GFS f006 + f012 (12-h lead)"}
POOL = re.compile(r"/(gdas|gfs)/(\d{8})/(\d\d)/\1\.t(\d\d)z\.1p00\.f(\d{3})(?:\.valid(\d{10}))?\.grb2$")
MIRROR = re.compile(r"/gdas\.(\d{8})/(\d\d)/atmos/gdas\.t(\d\d)z\.pgrb2\.1p00\.f(\d{3})$")


def decode(src):
    """source path -> (model, cycle datetime, forecast hour) or None"""
    m = POOL.search(src)
    if m:
        return m.group(1).upper(), datetime.strptime(m.group(2) + m.group(3), "%Y%m%d%H"), int(m.group(5))
    m = MIRROR.search(src)
    if m:
        return "GDAS", datetime.strptime(m.group(1) + m.group(2), "%Y%m%d%H"), int(m.group(4))
    return None


def collect(a):
    rows = []
    for p in a.products.split():
        for s in a.scenarios.split():
            for f in sorted(glob.glob(os.path.join(a.root, "records", p, s, "*.json"))):
                if os.path.basename(f).startswith("._"):
                    continue
                r = json.load(open(f))
                acq = datetime.fromisoformat(r["acquisition"])
                c1 = acq.replace(minute=0, second=0, microsecond=0, hour=acq.hour // 6 * 6)
                lads = sorted(r.get("lads", []), key=lambda x: x["name"])
                row = dict(product=p.upper(), scenario=s, granule=r["granule"], acquisition=r["acquisition"],
                           acq_h=round((acq - c1).total_seconds() / 3600, 3), n_files=len(lads),
                           second_weight=r.get("second_weight"), exit_status=r.get("exit_status"))
                for k, g in enumerate(lads[:2], start=1):
                    d = decode(g["source"])
                    name_valid = datetime.strptime(g["name"][-10:-1], "%y%m%d.%H")
                    row["f%d_name_valid_h" % k] = (name_valid - c1).total_seconds() / 3600
                    if d is None:
                        row["f%d_model" % k] = "?"
                        continue
                    model, cyc, fh = d
                    row["f%d_model" % k] = model
                    row["f%d_cycle_h" % k] = (cyc - c1).total_seconds() / 3600
                    row["f%d_fh" % k] = fh
                    row["f%d_valid_h" % k] = (cyc + timedelta(hours=fh) - c1).total_seconds() / 3600
                ok = len(lads) == 2
                for k, (model, coff, fh) in enumerate(RULE[s], start=1):
                    ok &= (row.get("f%d_model" % k) == model and row.get("f%d_cycle_h" % k) == coff
                           and row.get("f%d_fh" % k) == fh
                           and row.get("f%d_valid_h" % k) == row.get("f%d_name_valid_h" % k))
                if r.get("second_weight") is not None:
                    ok &= abs(r["second_weight"] - row["acq_h"] / 6) < 0.01
                row["follows_rule"] = int(ok)
                rows.append(row)
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    bad = [r for r in rows if not r["follows_rule"]]
    print("collected %d runs to %s; %d break the rule" % (len(rows), a.out, len(bad)))
    for r in bad[:10]:
        print("  RULE BROKEN:", r["product"], r["scenario"], r["granule"])


def plot(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    x = pd.concat([pd.read_csv(c, dtype={"granule": str}) for c in a.csv], ignore_index=True)
    os.makedirs(a.out, exist_ok=True)
    INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    MODEL = {"GDAS": "#2a78d6", "GFS": "#eb6834"}
    scen = [s for s in ("V0", "S1", "S2", "S3") if s in set(x["scenario"])]

    # ---------- 1. timeline for one granule: what each scenario feeds LaSRC
    g = x[x["granule"].str.contains(a.granule)] if a.granule else x[x["scenario"] == "S3"].head(1)
    gname = g["granule"].iloc[0]
    g = x[x["granule"] == gname].set_index("scenario")
    acq_h = float(g["acq_h"].iloc[0])
    fig, ax = plt.subplots(figsize=(11, 4.6), facecolor=SURF)
    ax.set_facecolor(SURF)
    for i, s in enumerate(reversed([s for s in scen if s in g.index])):
        r = g.loc[s]
        for k, dy in ((1, 0.13), (2, -0.13)):
            if pd.isna(r.get("f%d_cycle_h" % k)):
                continue
            c, v, fh, mdl = r["f%d_cycle_h" % k], r["f%d_valid_h" % k], int(r["f%d_fh" % k]), r["f%d_model" % k]
            col = MODEL.get(mdl, INK2)
            y = i + dy
            ax.plot([c, v], [y, y], color=col, lw=2.2, solid_capstyle="round", zorder=2)
            ax.scatter([c], [y], s=46, marker="s", color=col, zorder=3, edgecolor=SURF, linewidth=1)
            ax.scatter([v], [y], s=70, marker="o", color=SURF, edgecolor=col, linewidth=2, zorder=3)
            cn = {-6: "CP", 0: "C1", 6: "C2"}.get(c, "%+gh" % c)
            vn = {-6: "CP", 0: "C1", 6: "C2"}.get(v, "%+gh" % v)
            txt = "%s %s analysis (f000)" % (mdl, cn) if fh == 0 else "%s %s f%03d, valid at %s" % (mdl, cn, fh, vn)
            ax.text(v + 0.3, y, txt, va="center", fontsize=8.5, color=INK, zorder=4,
                    bbox=dict(boxstyle="square,pad=0.15", facecolor=SURF, edgecolor="none"))
    ax.axvline(acq_h, color="#d03b3b", lw=1.6, ls="--", zorder=1)
    ax.text(acq_h, len(scen) - 0.35, " acquisition", color="#d03b3b", fontsize=9, va="bottom")
    for t in (-6, 0, 6):
        ax.axvline(t, color=GRID, lw=0.9, zorder=0)
    ax.set_yticks(range(len(scen)))
    ax.set_yticklabels([LABEL[s] for s in reversed(scen)], fontsize=10, color=INK)
    ax.set_xticks([-6, 0, 6, 12])
    ax.set_xticklabels(["CP\n(C1 − 6 h)", "C1\n(cycle before\nacquisition)", "C2\n(C1 + 6 h)", "C1 + 12 h"],
                       fontsize=9, color=INK2)
    ax.set_xlim(-7, 15.5)
    ax.set_ylim(-0.6, len(scen) - 0.2)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(axis="y", length=0)
    from matplotlib.lines import Line2D
    h = [Line2D([0], [0], marker="s", color="w", markerfacecolor=MODEL["GDAS"], markersize=8, label="GDAS file issued (cycle)"),
         Line2D([0], [0], marker="s", color="w", markerfacecolor=MODEL["GFS"], markersize=8, label="GFS file issued (cycle)"),
         Line2D([0], [0], marker="o", color="w", markeredgecolor=INK2, markerfacecolor=SURF, markersize=8,
                markeredgewidth=2, label="time the file is valid for (LaSRC interpolates between the two)")]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=3, frameon=False, fontsize=9)
    ax.set_title("Which two GRIB files each scenario gave LaSRC (read from the run records)\nexample: %s, acquired %s UTC"
                 % (gname, str(g["acquisition"].iloc[0]).replace("T", " ")),
                 fontsize=11, fontweight="bold", color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "scenario_inputs_timeline.png"), dpi=150, facecolor=SURF, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)

    # ---------- 2. every run: observed (model, cycle, forecast hour) per file slot vs the rule
    prods = [p for p in ("L30", "S30") if p in set(x["product"])]
    fig, axes = plt.subplots(1, len(prods), figsize=(5.6 * len(prods), 4.4), facecolor=SURF, squeeze=False)
    for ax, prod in zip(axes[0], prods):
        ax.set_facecolor(SURF)
        d = x[x["product"] == prod]
        rows = []
        for s in scen:
            ds = d[d["scenario"] == s]
            if ds.empty:
                continue
            for k in (1, 2):
                combo = ds.apply(lambda r: "%s %+gh f%03d" % (r["f%d_model" % k], r["f%d_cycle_h" % k],
                                                             r["f%d_fh" % k]) if pd.notna(r.get("f%d_fh" % k)) else "?",
                                 axis=1).value_counts()
                m, c, fh = RULE[s][k - 1]
                want = "%s %+gh f%03d" % (m, c, fh)
                rows.append((s, k, want, combo.get(want, 0), len(ds), {kk: v for kk, v in combo.items() if kk != want}))
        y = list(range(len(rows)))[::-1]
        for yi, (s, k, want, n_ok, n, other) in zip(y, rows):
            ax.barh(yi, n_ok / n * 100, color="#1baf7a", height=0.62)
            if other:
                ax.barh(yi, sum(other.values()) / n * 100, left=n_ok / n * 100, color="#d03b3b", height=0.62)
            ax.text(101.5, yi, "%d / %d" % (n_ok, n) + ("   ✗ %s" % ", ".join("%s ×%d" % kv for kv in other.items()) if other else "  ✓"),
                    va="center", fontsize=8.5, color=INK)
        ax.set_yticks(y)
        ax.set_yticklabels(["%s file %d: %s" % (s, k, want) for s, k, want, *_ in rows], fontsize=8.5, color=INK)
        ax.set_xlim(0, 135)
        ax.set_xticks([0, 50, 100])
        ax.set_xlabel("% of runs that used the expected file", fontsize=9, color=INK2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=INK2, labelsize=8.5)
        nb = int((d["follows_rule"] == 0).sum())
        ax.set_title("%s: %d runs, %s" % (prod, len(d), "all used the expected files" if nb == 0 else "%d used other files" % nb),
                     fontsize=10.5, fontweight="bold", color=INK, loc="left")
    fig.suptitle("Expected file = model, issue cycle relative to C1, forecast hour (green = used; red = anything else)",
                 fontsize=9.5, color=INK2, y=0.025)
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(os.path.join(a.out, "scenario_inputs_check.png"), dpi=150, facecolor=SURF, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)

    t = x.groupby(["product", "scenario"])["follows_rule"].agg(["sum", "count"])
    print("| Product | Scenario | Runs following the rule |\n|---|---|---|")
    for (p, s), r in t.iterrows():
        print("| %s | %s | %d of %d |" % (p, s, r["sum"], r["count"]))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("collect")
    c.add_argument("--root", default="/rtmp/tvo/LaSRC_ancillary_latency_20261120")
    c.add_argument("--products", default="l30 s30")
    c.add_argument("--scenarios", default="V0 S1 S2 S3")
    c.add_argument("--out", required=True)
    q = sub.add_parser("plot")
    q.add_argument("--csv", action="append", required=True)
    q.add_argument("--out", required=True)
    q.add_argument("--granule", default=None)
    a = ap.parse_args()
    collect(a) if a.cmd == "collect" else plot(a)


if __name__ == "__main__":
    main()
