#!/usr/bin/env bash
# Run inside the activated GATEROOT conda environment on WSL/Linux.
set -euo pipefail
: "${CONDA_PREFIX:?Activate the GATEROOT conda environment first}"
work="${1:-$HOME/gateroot-sdk-0.13.3}"
mkdir -p "$work"
if [[ ! -d "$work/PETSIRD/.git" ]]; then
  git clone --branch v0.13.3 --depth 1 https://github.com/ETSInitiative/PETSIRD.git "$work/PETSIRD"
fi
[[ "$(git -C "$work/PETSIRD" describe --tags --exact-match)" == v0.13.3 ]]
mkdir -p "$work/yardl"
curl -fsSL https://github.com/microsoft/yardl/releases/download/v0.6.7/yardl_0.6.7_linux_x86_64.tar.gz | tar -xz -C "$work/yardl"
(cd "$work/PETSIRD/model" && "$work/yardl/yardl" generate)
cmake -S "$work/PETSIRD/cpp" -B "$work/build" -G Ninja -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX"
cmake --build "$work/build" --parallel 2
cmake --install "$work/build"
python -m pip install --no-deps "$work/PETSIRD/python"
