#!/bin/bash
set -o errexit
set -o nounset
set -o pipefail

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 GRANULE_ID" >&2
  exit 1
fi

granule="$1"
gfs_source_dir="${GFS_SOURCE_DIR:?set GFS_SOURCE_DIR to a local mirror of s3://noaa-gfs-bdp-pds}"
aux_lads_dir="${AUX_LADS_DIR:?set AUX_LADS_DIR to $LASRC_AUX_DIR/LADS}"

mkdir -p "$aux_lads_dir"

python3 - "$granule" "$gfs_source_dir" "$aux_lads_dir" <<'PY'
import datetime as dt
import shutil
import pathlib
import re
import sys

granule = sys.argv[1]
gfs_source_dir = pathlib.Path(sys.argv[2])
aux_lads_dir = pathlib.Path(sys.argv[3])

match = re.search(r"_(\d{8}T\d{6})_", granule)
if not match:
    raise SystemExit(f"Could not parse scene timestamp from granule: {granule}")

scene_dt = dt.datetime.strptime(match.group(1), "%Y%m%dT%H%M%S")
candidates = []
candidate_days = {
    (scene_dt - dt.timedelta(days=1)).strftime("%Y%m%d"),
    scene_dt.strftime("%Y%m%d"),
    (scene_dt + dt.timedelta(days=1)).strftime("%Y%m%d"),
}

for day in sorted(candidate_days):
    for path in sorted(
        gfs_source_dir.glob(f"gdas.{day}/*/atmos/gdas.t??z.pgrb2.1p00.f000")
    ):
        date_match = re.search(
            r"gdas\.(\d{8})/(\d{2})/atmos/gdas\.t\d{2}z\.pgrb2\.1p00\.f000$",
            str(path),
        )
        if not date_match:
            continue
        file_dt = dt.datetime.strptime("".join(date_match.groups()), "%Y%m%d%H")
        candidates.append((file_dt, path))

if len(candidates) < 2:
    raise SystemExit(
        "Found fewer than 2 GDAS files within +/- 1 day of "
        f"{scene_dt:%Y-%m-%d %H:%M:%S} under {gfs_source_dir}"
    )

before = [item for item in candidates if item[0] <= scene_dt]
after = [item for item in candidates if item[0] >= scene_dt]

if before and after:
    first = before[-1]
    second = after[0]
    if first[1] == second[1]:
        idx = candidates.index(first)
        if idx > 0:
            first = candidates[idx - 1]
        elif idx + 1 < len(candidates):
            second = candidates[idx + 1]
else:
    first, second = candidates[-2], candidates[-1]

selected = [first, second]
for file_dt, source_path in selected:
    link_name = f"gdas1.PGrbF00.{file_dt.strftime('%y%m%d')}.{file_dt.strftime('%H')}z"
    staged_path = aux_lads_dir / link_name
    if staged_path.exists():
        print(f"Keeping existing staged file {staged_path}")
        continue
    shutil.copy2(source_path, staged_path)
    print(f"Copied {source_path} -> {staged_path}")
PY
