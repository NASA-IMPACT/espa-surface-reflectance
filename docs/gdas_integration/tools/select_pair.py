#!/usr/bin/env python3
"""Pick the two GRIB files for one granule and scenario (large-scale ancillary test).

C1 = the 6-hourly cycle at/before the acquisition, C2 = C1 + 6 h, CP = C1 - 6 h.
  V0  GDAS C1 f000 (mirror)          + GDAS C2 f000 (mirror)
  S1  GDAS C1 f000 (mirror)          + GDAS C1 f006 relabelled to C2 (pool)
  S2  GFS  C1 f000 (pool)            + GFS  C1 f006 relabelled to C2 (pool)
  S3  GFS  CP f006 relabelled to C1  + GFS  CP f012 relabelled to C2 (pool)
Prints one line per file:  <LADS name>  <source path>
LADS names are gdas1.PGrbF00.YYMMDD.HHz by valid time, as the wrappers expect.

Usage: select_pair.py --acq 2025-01-17T14:37:41 --scenario S2 --pool POOL --gdas-mirror MIRROR
"""
import argparse
import os
import sys
from datetime import datetime, timedelta


def lads(t):
    return "gdas1.PGrbF00.%s.%sz" % (t.strftime("%y%m%d"), t.strftime("%H"))


def mirror(m, t):
    return os.path.join(m, "gdas.%s" % t.strftime("%Y%m%d"), t.strftime("%H"), "atmos",
                        "gdas.t%sz.pgrb2.1p00.f000" % t.strftime("%H"))


def pool(p, model, cyc, fh, valid=None):
    name = "%s.t%sz.1p00.f%03d" % (model, cyc.strftime("%H"), fh)
    if valid is not None:
        name += ".valid%s" % valid.strftime("%Y%m%d%H")
    return os.path.join(p, model, cyc.strftime("%Y%m%d"), cyc.strftime("%H"), name + ".grb2")


def pair(acq, scen, p, m):
    c1 = acq.replace(minute=0, second=0, microsecond=0, hour=acq.hour // 6 * 6)
    c2, cp = c1 + timedelta(hours=6), c1 - timedelta(hours=6)
    if scen == "V0":
        f1, f2 = mirror(m, c1), mirror(m, c2)
    elif scen == "S1":
        f1, f2 = mirror(m, c1), pool(p, "gdas", c1, 6, c2)
    elif scen == "S2":
        f1, f2 = pool(p, "gfs", c1, 0), pool(p, "gfs", c1, 6, c2)
    elif scen == "S3":
        f1, f2 = pool(p, "gfs", cp, 6, c1), pool(p, "gfs", cp, 12, c2)
    else:
        raise SystemExit("unknown scenario %s" % scen)
    return [(lads(c1), f1), (lads(c2), f2)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--acq", required=True)
    ap.add_argument("--scenario", required=True, choices=["V0", "S1", "S2", "S3"])
    ap.add_argument("--pool", required=True)
    ap.add_argument("--gdas-mirror", required=True)
    a = ap.parse_args()
    acq = datetime.fromisoformat(a.acq.replace("Z", ""))
    rc = 0
    for name, src in pair(acq, a.scenario, a.pool, a.gdas_mirror):
        if not os.path.isfile(src):
            sys.stderr.write("MISSING %s\n" % src)
            rc = 3
        print(name, src)
    sys.exit(rc)


if __name__ == "__main__":
    main()
