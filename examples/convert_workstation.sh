#!/usr/bin/env bash
set -euo pipefail
root="${1:-/mnt/f/DATA/WCMC/ETSIhackathon}"
destination="${2:?Supply a fresh output directory as argument 2}"
execution="${3:-auto}"
workers="${4:-4}"
repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
converter="${ROOT_TO_PETSIRD:-$repo/cpp/build/root_to_petsird}"
mkdir -p "$destination"
for name in IEC_6min Hoffman_3min Hoffman_6min; do
  if [[ "$name" == IEC_6min ]]; then
    prefix="$root/ROOT_files/IEC/ETSIPETscanner_mIEC_"
    geometry=ETSIPETscanner_geometry.json
    norm=Sino_NCFs_ETSIPETscanner_span1_mrd39_df.Cdf
  else
    prefix="$root/ROOT_files/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_"
    geometry=ETSIPETscanner256mmAFOV_geometry.json
    norm=Sino_NCFs_ETSIPETscanner256mmAFOV_span1_mrd79_df.Cdf
  fi
  stop=360
  [[ "$name" != Hoffman_3min ]] || stop=180
  "$converter" --root-prefix "$prefix" --number-of-root-files 36 \
    -s "$root/GEOM_DATA/$geometry" -c "$root/NORM_files/$norm" \
    --petsird-file "$destination/${name}_prompts_delayeds_norm.petsird" \
    --execution "$execution" --workers "$workers" --time-mode auto --file-duration-ms 10000 \
    --selection all --delayeds required --delayed-offset-ps 500000 \
    --start-time-s 0 --end-time-s "$stop"
done
printf 'Completed three C++ conversions in %s\n' "$destination"
