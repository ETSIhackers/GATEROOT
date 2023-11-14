set shell := ['bash', '-ceuo', 'pipefail']

default: run

@ensure-build-dir:
    mkdir -p cpp/build

@generate:
    cd PETSIRD/model; \
    yardl generate

@configure: generate ensure-build-dir
    cd cpp/build; \
    cmake -GNinja ..

@build: generate configure
    cd cpp/build; \
    ninja

@download-test-data:
    export AZURE_STORAGE_SAS_TOKEN="sp=rl&st=2023-11-14T18:39:07Z&se=2024-11-15T02:39:07Z&spr=https&sv=2022-11-02&sr=c&sig=0KVD7ORBM7Mx1%2BhVrVbqYcQycshhvT2XvdmrWVetiQM%3D"; \
    dvc pull

@run: build download-test-data
    cpp/build/root_to_petsird --root-prefix data/root/ETSIPETscanner_mIEC_ --petsird-file /dev/null
