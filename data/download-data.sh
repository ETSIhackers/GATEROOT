#!/bin/bash

set -euo pipefail

# Check if azcopy is installed
cd "$(dirname "$0")"

if ! command -v azcopy &> /dev/null
then
    wget -O azcopy.tar.gz https://aka.ms/downloadazcopy-v10-linux
    tar -xzvf azcopy.tar.gz
    AZCOPY_EXE=$(find ./ -name azcopy)
    cp $AZCOPY_EXE "${CONDA_PREFIX}/bin/"
fi

DATA_SAS_URI="https://etsinitiative.blob.core.windows.net/root?sp=rl&st=2023-11-13T23:17:25Z&se=2023-11-16T07:17:25Z&spr=https&sv=2022-11-02&sr=c&sig=nPzgxNOpaparisn8RG3rLllh%2F5PHaJmg1HyWltgjpgc%3D"
azcopy cp "${DATA_SAS_URI}" ./ --recursive
