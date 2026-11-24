#!/usr/bin/env python3
"""Scene-centre water vapour (uwv) and ozone (uoz) that LaSRC uses, computed from the two GRIB2
files of a scenario WITHOUT running the container.

Used for V0 (the original GDAS-analysis run), whose outputs were produced before the per-granule
records existed. The same calculation is checked against the values LaSRC itself logged for S1,
S2 and S3 (--check), so it can be trusted for V0.

What LaSRC b271310 does with the two GRIB files (GRIB path, whichone=0), reproduced here:
  1. read PWAT (kg/m2) and TOZNE (DU) on the 1-degree global grid of each file;
  2. time-interpolate: field = w1 * file1 + w2 * file2, with w2 = (acquisition - C1) / 6 h
     (the "First/Second weight" in the LaSRC log);
  3. bilinear-resample to the 0.05-degree CMG grid; the grid is read with a 1-degree eastward
     offset (KNOWN_ISSUES #2), reproduced as a roll of the field by -(180 + 1) columns;
  4. store as integers: wv = int(PWAT * 10) (0.1 kg/m2), oz = int(TOZNE * 0.4) (2.5 DU);
  5. take the CMG cell nearest the scene centre: uwv = wv / 200 (g/cm2), uoz = oz / 400 (atm-cm).
Scene centre and weights are taken from the S1 record of the same granule (same scene, same
acquisition time, same C1/C2 valid times, so the same weights as V0).

GRIB2 decoding uses numpy only (templates 5.0 simple and 5.3 complex packing with spatial
differencing, which NCEP uses for pgrb2.1p00); no eccodes needed.

Usage:
  scene_centre_values.py --root ROOT --scenario V0 [--write]   # V0 values (needs the GDAS mirror)
  scene_centre_values.py --root ROOT --check                   # compare with LaSRC logs for S1-S3
  scene_centre_values.py --root ROOT --table OUT.csv           # V0-S3 values before and after LaSRC's
                                                               # integer storage (precision-loss figure)
--write saves records/<p>/V0/<granule>.json (fields as run_record.py, plus "source"; exit_status is
null so run_scenario.slurm would not mistake it for a finished container run).

Check on the pilot (S2 and S3, 40 runs, GFS files): 33 of 40 reproduce LaSRC's logged uwv and uoz
exactly; the other 7 differ by exactly one storage step in one of the two values (0.005 g/cm2 or
2.5 DU). In all 7 the interpolated value lies within 0.15 of an integer storage boundary, so the
difference comes from floating-point detail in LaSRC's resampling, not from a different method.
"""
import argparse
import glob
import json
import os
import struct
import sys
from datetime import datetime

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from select_pair import pair  # noqa: E402  (same file choice as the container runs)

NAMES = {(0, 1, 3): "PWAT", (0, 14, 0): "TOZNE"}


# ---------------------------------------------------------------- GRIB2 decoding (numpy only)
def _bits(buf, bitpos, nbits, count):
    """count unsigned integers of nbits each, starting at bit offset bitpos of buf."""
    if nbits == 0 or count == 0:
        return np.zeros(count, dtype=np.int64), bitpos
    b0, b1 = bitpos // 8, (bitpos + nbits * count + 7) // 8 + 1
    bits = np.unpackbits(np.frombuffer(buf[b0:b1] + b"\0", dtype=np.uint8))
    off = bitpos - 8 * b0
    bits = bits[off:off + nbits * count].reshape(count, nbits).astype(np.int64)
    return bits @ (1 << np.arange(nbits - 1, -1, -1, dtype=np.int64)), bitpos + nbits * count


def _sm(v, nbytes):
    """sign-magnitude integer (GRIB2 convention for scale factors and extra descriptors)."""
    top = 1 << (8 * nbytes - 1)
    return -(v & (top - 1)) if v & top else v


def _decode(s5, s7, npts):
    tmpl = struct.unpack(">H", s5[9:11])[0]
    R = struct.unpack(">f", s5[11:15])[0]
    E = _sm(struct.unpack(">H", s5[15:17])[0], 2)
    D = _sm(struct.unpack(">H", s5[17:19])[0], 2)
    nb = s5[19]
    data = s7[5:]
    if tmpl == 0:
        X, _ = _bits(data, 0, nb, npts)
    elif tmpl in (2, 3):
        if s5[22] != 0:
            raise SystemExit("missing-value management %d not supported" % s5[22])
        NG = struct.unpack(">I", s5[31:35])[0]
        wref, wbits = s5[35], s5[36]
        lref, linc = struct.unpack(">I", s5[37:41])[0], s5[41]
        llast, lbits = struct.unpack(">I", s5[42:46])[0], s5[46]
        order, nex = (s5[47], s5[48]) if tmpl == 3 else (0, 0)
        pos, extras = 0, []
        for _ in range(order + 1 if order else 0):
            extras.append(_sm(int.from_bytes(data[pos:pos + nex], "big"), nex)); pos += nex
        bp = pos * 8
        gref, bp = _bits(data, bp, nb, NG); bp = (bp + 7) // 8 * 8
        gw, bp = _bits(data, bp, wbits, NG); bp = (bp + 7) // 8 * 8
        gw = gw + wref
        gl, bp = _bits(data, bp, lbits, NG); bp = (bp + 7) // 8 * 8
        gl = lref + gl * linc
        gl[-1] = llast
        X = np.empty(int(gl.sum()), dtype=np.int64)
        k = 0
        for r, w, L in zip(gref, gw, gl):
            v, bp = _bits(data, bp, int(w), int(L))
            X[k:k + L] = v + r
            k += L
        if order:
            minsd = extras[-1]
            X[order:] += minsd
            if order == 1:
                X[0] = extras[0]
                X = np.cumsum(X)
            else:
                X[0], X[1] = extras[0], extras[1]
                for i in range(2, len(X)):           # second-order (rare for these fields)
                    X[i] += 2 * X[i - 1] - X[i - 2]
    else:
        raise SystemExit("data representation template 5.%d not supported" % tmpl)
    return (R + X.astype(np.float64) * 2.0 ** E) / 10.0 ** D


def read_fields(path):
    """{'PWAT': 181x360, 'TOZNE': 181x360} from the 'entire atmosphere' messages."""
    d, pos, out = open(path, "rb").read(), 0, {}
    while True:
        i = d.find(b"GRIB", pos)
        if i < 0:
            return out
        tot = struct.unpack(">Q", d[i + 8:i + 16])[0]
        m, q, disc = d[i:i + tot], 16, d[i + 6]
        sec, name, surf = {}, None, None
        while q < len(m) - 4 and m[q:q + 4] != b"7777":
            L = struct.unpack(">I", m[q:q + 4])[0]
            sec[m[q + 4]] = m[q:q + L]
            q += L
        if 4 in sec:
            name, surf = NAMES.get((disc, sec[4][9], sec[4][10])), sec[4][22]
        if name and surf == 200 and name not in out:
            nx, ny = struct.unpack(">II", sec[3][30:38])
            if sec[6][5] != 255:
                raise SystemExit("bitmap not supported in %s" % path)
            npts = struct.unpack(">I", sec[5][5:9])[0]
            out[name] = _decode(sec[5], sec[7], npts).reshape(ny, nx)
        pos = i + tot


# ---------------------------------------------------------------- LaSRC scene-centre logic
def _linesamp(lat, lon):
    il = int((90 - lat) + 0.5); rf = (90 - il) - lat
    if rf < 0:
        il -= 1; rf = 1 - abs(rf)
    isamp = int((lon + 180) + 0.5); cf = lon + (180 - isamp)
    if cf < 0:
        isamp -= 1; cf = 1 - abs(cf)
    return il, isamp, rf, cf


def lasrc_scene_centre(f1, f2, w1, w2, lat, lon, shift=1):
    A, B = read_fields(f1), read_fields(f2)
    lcmg = round((89.975 - lat) * 20); scmg = round((179.975 + lon) * 20)
    il, isamp, rf, cf = _linesamp(90 - lcmg * 0.05, -180 + scmg * 0.05)
    j1 = (isamp + 1) % 360
    val = {}
    for v in ("PWAT", "TOZNE"):
        g = np.roll(A[v].astype(np.float32) * np.float32(w1) + B[v].astype(np.float32) * np.float32(w2),
                    -(180 + shift), axis=1)
        val[v] = float((g[il, isamp] * (1 - cf) + g[il, j1] * cf) * (1 - rf)
                       + (g[il + 1, isamp] * (1 - cf) + g[il + 1, j1] * cf) * rf)
    wv, oz = int(val["PWAT"] * 10), int(val["TOZNE"] * 0.4)
    return dict(uwv=wv / 200, uoz=oz / 400, pwat_interp=val["PWAT"], tozne_interp=val["TOZNE"])


def local_path(src, root):
    """Records hold Meteor paths; map /nas/rtmp/... pool paths onto --root when needed."""
    if os.path.exists(src):
        return src
    if "ancillary_pool" in src:
        return os.path.join(root, "ancillary_pool" + src.split("ancillary_pool", 1)[1])
    return src


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/rtmp/tvo/LaSRC_ancillary_latency_20261120")
    ap.add_argument("--scenario", default="V0", choices=["V0", "S1", "S2", "S3"])
    ap.add_argument("--gdas-mirror", default="/nas/rstor/tvo/GFS_2025")
    ap.add_argument("--check", action="store_true", help="compare with the uwv/uoz LaSRC logged for S1-S3")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--granule", default=None, help="print the section-5.1 table rows for one granule (e.g. T16TCK)")
    ap.add_argument("--table", default=None, metavar="CSV",
                    help="write continuous (before integer storage) and LaSRC-used values for V0-S3")
    ap.add_argument("--products", default="s30 l30")
    a = ap.parse_args()

    if a.check:
        n = ok = 0
        for f in sorted(glob.glob(os.path.join(a.root, "records", "*", "S[123]", "*.json"))):
            if os.path.basename(f).startswith("._"):
                continue
            r = json.load(open(f))
            if r.get("exit_status") != 0 or r.get("uwv") is None:
                continue
            files = sorted(r["lads"], key=lambda x: x["name"])
            paths = [local_path(x["source"], a.root) for x in files]
            if not all(os.path.exists(p) for p in paths):
                continue
            v = lasrc_scene_centre(paths[0], paths[1], r["first_weight"], r["second_weight"],
                                   r["scene_lat"], r["scene_lon"])
            match = abs(v["uwv"] - r["uwv"]) < 1e-9 and abs(v["uoz"] - r["uoz"]) < 1e-9
            n += 1; ok += match
            print("%s %s %-62s uwv %.3f/%.3f uoz %.4f/%.4f %s" % (r["product"], r["scenario"], r["granule"],
                  v["uwv"], r["uwv"], v["uoz"], r["uoz"], "OK" if match else "DIFF"))
        print("check: %d of %d runs reproduce LaSRC's logged uwv and uoz exactly" % (ok, n))
        return

    if a.granule:
        # one granule, all scenarios: value at the scene centre (time interpolation only) and in
        # LaSRC's grid cell (time + spatial interpolation, what gets stored), then stored and used
        f = glob.glob(os.path.join(a.root, "records", "*", "S1", "*%s*.json" % a.granule))
        if not f:
            raise SystemExit("no S1 record matching %s" % a.granule)
        r = json.load(open(f[0]))
        acq = datetime.fromisoformat(r["acquisition"])
        print("| Granule | Run | PWAT at scene centre | PWAT LaSRC's cell | wv | uwv | TOZNE at scene centre | "
              "TOZNE LaSRC's cell | oz | uoz |")
        for scen in ("V0", "S1", "S2", "S3"):
            (n1, f1), (n2, f2) = pair(acq, scen, os.path.join(a.root, "ancillary_pool"), a.gdas_mirror)
            c = lasrc_scene_centre(f1, f2, r["first_weight"], r["second_weight"], r["scene_lat"], r["scene_lon"], 1)
            t = lasrc_scene_centre(f1, f2, r["first_weight"], r["second_weight"], r["scene_lat"], r["scene_lon"], 0)
            print("| %s | %s | %.3f | %.3f | %d | %.3f | %.2f | %.2f | %d | %.4f |"
                  % (a.granule, scen, t["pwat_interp"], c["pwat_interp"], int(c["pwat_interp"] * 10), c["uwv"],
                     t["tozne_interp"], c["tozne_interp"], int(c["tozne_interp"] * 0.4), c["uoz"]))
        return

    if a.table:
        # one row per granule x scenario: value before LaSRC's integer storage (continuous) and
        # after (what LaSRC uses); the logged value is added where LaSRC printed one
        import csv
        rows = []
        for p in a.products.split():
            for f in sorted(glob.glob(os.path.join(a.root, "records", p, "S1", "*.json"))):
                if os.path.basename(f).startswith("._"):
                    continue
                r = json.load(open(f))
                if r.get("scene_lat") is None or r.get("first_weight") is None:
                    continue
                acq = datetime.fromisoformat(r["acquisition"])
                for scen in ("V0", "S1", "S2", "S3"):
                    (n1, f1), (n2, f2) = pair(acq, scen, os.path.join(a.root, "ancillary_pool"), a.gdas_mirror)
                    if not (os.path.exists(f1) and os.path.exists(f2)):
                        print("missing GRIB for %s %s" % (scen, r["granule"])); continue
                    v = lasrc_scene_centre(f1, f2, r["first_weight"], r["second_weight"], r["scene_lat"], r["scene_lon"])
                    logf = os.path.join(a.root, "records", p, scen, r["granule"] + ".json")
                    lg = json.load(open(logf)) if os.path.exists(logf) and scen != "V0" else {}
                    rows.append(dict(product=r["product"], granule=r["granule"], scenario=scen,
                                     scene_lat=r["scene_lat"], scene_lon=r["scene_lon"],
                                     pwat_continuous_kgm2=round(v["pwat_interp"], 5),
                                     uwv_continuous=round(v["pwat_interp"] / 10, 6),
                                     uwv_lasrc=v["uwv"], uwv_logged=lg.get("uwv"),
                                     tozne_continuous_du=round(v["tozne_interp"], 4),
                                     uoz_lasrc_du=round(v["uoz"] * 1000, 2),
                                     uoz_logged_du=round(lg["uoz"] * 1000, 2) if lg.get("uoz") is not None else None,
                                     file1=f1, file2=f2))
                print("table: %s %s" % (r["product"], r["granule"]))
        out = a.table
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
        print("wrote %s (%d rows)" % (out, len(rows)))
        return

    for p in a.products.split():
        for f in sorted(glob.glob(os.path.join(a.root, "records", p, "S1", "*.json"))):
            if os.path.basename(f).startswith("._"):
                continue
            r = json.load(open(f))
            if r.get("scene_lat") is None or r.get("first_weight") is None:
                print("skip %s (no scene centre/weights in S1 record)" % r["granule"]); continue
            acq = datetime.fromisoformat(r["acquisition"])
            (n1, f1), (n2, f2) = pair(acq, a.scenario, os.path.join(a.root, "ancillary_pool"), a.gdas_mirror)
            v = lasrc_scene_centre(f1, f2, r["first_weight"], r["second_weight"], r["scene_lat"], r["scene_lon"])
            print("%s %s %-62s uwv %.3f uoz %.4f  (%s, %s)" % (r["product"], a.scenario, r["granule"], v["uwv"], v["uoz"],
                  os.path.basename(f1), os.path.basename(f2)))
            if a.write:
                rec = dict(product=r["product"], scenario=a.scenario, granule=r["granule"], acquisition=r["acquisition"],
                           exit_status=None, source="reproduced: scene_centre_values.py from the GRIB pair, no container run",
                           lads=[dict(name=n1, source=f1), dict(name=n2, source=f2)],
                           first_weight=r["first_weight"], second_weight=r["second_weight"],
                           scene_lat=r["scene_lat"], scene_lon=r["scene_lon"], uwv=v["uwv"], uoz=v["uoz"],
                           pres=r.get("pres"), pwat_interp=v["pwat_interp"], tozne_interp=v["tozne_interp"])
                out = os.path.join(a.root, "records", p, a.scenario, r["granule"] + ".json")
                os.makedirs(os.path.dirname(out), exist_ok=True)
                json.dump(rec, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
