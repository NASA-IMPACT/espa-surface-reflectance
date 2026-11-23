#!/usr/bin/env python3
"""Audit that every run used the ancillary (GRIB) files its scenario is supposed to use.

For each run record (records/<product>/<scenario>/<granule>.json, written by run_record.py) it checks:
  1. expected   the two files the scenario should use (select_pair.py, from the acquisition time)
  2. names      the LADS names staged for LaSRC (gdas1.PGrbF00.YYMMDD.HHz) match the expected valid times
  3. sources    each staged file links to the expected source file (pool or GDAS mirror)
  4. read       LaSRC's log reports the same pair ("AUX input file" / "AUX2 input file")
  5. checksum   the md5 of each pool file matches ancillary_pool/manifest.csv
  6. grib time  the reference time inside each pool file (manifest) equals the time in its LADS name,
                i.e. forecast files were relabelled and LaSRC reads them at the right time
  7. weights    LaSRC's time weights match (acquisition - C1) / 6 h
  8. distinct   the scenarios of one granule did not use the same pair of files
  9. values     LaSRC logged uwv/uoz and the run finished (exit status 0)

Usage:  audit_ancillary_inputs.py --root /rtmp/tvo/LaSRC_ancillary_latency_20261120 [--products "s30 l30"]
        [--scenarios "S1 S2 S3"] [--granule T16TCK] [--md5-mirror] [--out audit.csv]
--md5-mirror also recomputes the md5 of the GDAS-mirror files (slow: ~40 MB each).
Exit code 0 if no problem was found, 1 otherwise. Problems are listed and written to --out.
"""
import argparse
import csv
import glob
import hashlib
import json
import os
import sys
from datetime import datetime, timedelta

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from select_pair import pair  # noqa: E402  (the same rule the runner used)


def tail(path, marker):
    """part of a path after a marker directory (records store /nas/rtmp/..., runs use /rtmp/...)."""
    return path.split(marker, 1)[1] if marker in path else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/rtmp/tvo/LaSRC_ancillary_latency_20261120")
    ap.add_argument("--gdas-mirror", default="/nas/rstor/tvo/GFS_2025")
    ap.add_argument("--products", default="s30 l30")
    ap.add_argument("--scenarios", default="S1 S2 S3")
    ap.add_argument("--granule", default=None, help="only granules whose name contains this text")
    ap.add_argument("--md5-mirror", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    pool = os.path.join(a.root, "ancillary_pool")
    manifest = {r["path"]: r for r in csv.DictReader(open(os.path.join(pool, "manifest.csv")))}
    problems, n_runs, pairs = [], 0, {}
    mirror_md5 = {}

    def bad(rec, check, msg):
        problems.append(dict(product=rec["product"], scenario=rec["scenario"], granule=rec["granule"],
                             check=check, detail=msg))

    for p in a.products.split():
        for s in a.scenarios.split():
            for f in sorted(glob.glob(os.path.join(a.root, "records", p, s, "*.json"))):
                if os.path.basename(f).startswith("._") or (a.granule and a.granule not in f):
                    continue
                rec = json.load(open(f))
                n_runs += 1
                if rec.get("exit_status") != 0:
                    bad(rec, "values", "run did not finish (exit status %s)" % rec.get("exit_status"))
                    continue
                acq = datetime.fromisoformat(rec["acquisition"])
                exp = pair(acq, s, pool, a.gdas_mirror)                     # [(name, path), (name, path)]
                got = sorted(rec.get("lads", []), key=lambda x: x["name"])
                if len(got) != 2:
                    bad(rec, "names", "%d staged files instead of 2" % len(got)); continue
                for (en, ep), g in zip(exp, got):
                    if g["name"] != en:
                        bad(rec, "names", "staged %s, expected %s" % (g["name"], en))
                    # source: compare the part after ancillary_pool/ or GFS_2025/
                    for marker in ("ancillary_pool/", os.path.basename(a.gdas_mirror.rstrip("/")) + "/"):
                        te, tg = tail(ep, marker), tail(g["source"], marker)
                        if te is not None:
                            if tg != te:
                                bad(rec, "sources", "%s links to %s, expected %s" % (g["name"], g["source"], ep))
                            break
                    rel = tail(g["source"], "ancillary_pool/")
                    if rel is not None:                                    # pool file
                        m = manifest.get(rel)
                        if m is None:
                            bad(rec, "checksum", "%s not in manifest" % rel)
                        else:
                            if m["md5"] != g.get("md5"):
                                bad(rec, "checksum", "%s md5 %s, manifest %s" % (rel, g.get("md5"), m["md5"]))
                            want = datetime.strptime(g["name"][-10:-1], "%y%m%d.%H").strftime("%Y-%m-%dT%H")
                            if m["grib_reference"] != want or str(m["grib_forecast_h"]) not in ("0", "0.0"):
                                bad(rec, "grib time", "%s reference %s +%sh, LaSRC needs %s +0h"
                                    % (rel, m["grib_reference"], m["grib_forecast_h"], want))
                    elif a.md5_mirror:                                      # GDAS mirror file
                        src = ep
                        if src not in mirror_md5:
                            mirror_md5[src] = hashlib.md5(open(src, "rb").read()).hexdigest()
                        if mirror_md5[src] != g.get("md5"):
                            bad(rec, "checksum", "%s md5 %s, now %s" % (src, g.get("md5"), mirror_md5[src]))
                if [rec.get("lasrc_aux"), rec.get("lasrc_aux2")] != [e[0] for e in exp]:
                    bad(rec, "read", "LaSRC read %s / %s" % (rec.get("lasrc_aux"), rec.get("lasrc_aux2")))
                c1 = acq.replace(minute=0, second=0, microsecond=0, hour=acq.hour // 6 * 6)
                w2 = (acq - c1).total_seconds() / 21600
                if rec.get("second_weight") is None or abs(rec["second_weight"] - w2) > 0.01:
                    bad(rec, "weights", "second weight %s, expected %.4f" % (rec.get("second_weight"), w2))
                if rec.get("uwv") is None or rec.get("uoz") is None:
                    bad(rec, "values", "no uwv/uoz in the LaSRC log")
                pairs.setdefault((p, rec["granule"]), {})[s] = tuple(g.get("md5") for g in got)

    for (p, g), d in pairs.items():
        seen = {}
        for s, md in d.items():
            if md in seen:
                problems.append(dict(product=p.upper(), scenario=s, granule=g, check="distinct",
                                     detail="same files as %s" % seen[md]))
            seen[md] = s

    print("audited %d runs: %d problem(s)" % (n_runs, len(problems)))
    by = {}
    for x in problems:
        by.setdefault(x["check"], []).append(x)
    for c in ("values", "names", "sources", "read", "checksum", "grib time", "weights", "distinct"):
        print("  %-9s %s" % (c, "OK" if c not in by else "%d problem(s)" % len(by[c])))
        for x in by.get(c, [])[:5]:
            print("      %s %s %s: %s" % (x["product"], x["scenario"], x["granule"], x["detail"]))
    if a.out and problems:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(problems[0])); w.writeheader(); w.writerows(problems)
        print("details in", a.out)
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
