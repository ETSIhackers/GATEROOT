#! /bin/bash
# shellcheck source=/dev/null

source /opt/conda/etc/profile.d/conda.sh
conda activate prd
source <(yardl completion bash)
source <(just --completions bash)

if [[ "${BASH_ENV:-}" == "$(readlink -f "${BASH_SOURCE[0]:-}")" ]]; then
    # We don't want subshells to unnecessarily source this again.
    unset BASH_ENV
fi
