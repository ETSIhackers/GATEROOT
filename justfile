set shell := ['bash', '-ceuo', 'pipefail']

@download-test-data:
    data/download-data.sh

@ensure-build-dir:
    mkdir -p cpp/build

@generate:
    cd PETSIRD/model; \
    yardl generate

@configure: generate ensure-build-dir
    cd cpp/build; \
    cmake -GNinja ..

@build: configure
    cd cpp/build; \
    ninja
