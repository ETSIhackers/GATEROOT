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

@build: configure
    cd cpp/build; \
    ninja


#@run: build download-test-data
#    cpp/build/root_to_petsird --root-prefix data/root/ETSIPETscanner_mIEC_ --number-of-root-files 1 --petsird-file petsird.bin -s data/root/scanner_geometry.json

#@run-full: build download-test-data-full
#    cpp/build/root_to_petsird --root-prefix data/root/ETSIPETscanner_mIEC_ --number-of-root-files 36 --petsird-file petsird-full.bin -s data/root/scanner_geometry.json



@run_0min_mIEC: build
    cpp/build/root_to_petsird --root-prefix ROOT_DATA/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 0 --petsird-file petsird_ETSIPETscanner_mIEC_0min_v0.2.bin -s data/root/scanner_geometry.json

@run_10sec_mIEC: build
    cpp/build/root_to_petsird --root-prefix ROOT_DATA/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 1 --petsird-file petsird_ETSIPETscanner_mIEC_10sec_v0.2.bin -s data/root/scanner_geometry.json

@run-6min_mIEC: build
    cpp/build/root_to_petsird --root-prefix ROOT_DATA/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 36 --petsird-file petsird_ETSIPETscanner_mIEC_6min_v0.2.bin -s data/root/scanner_geometry.json


@run_0min_voxBrain: build
    cpp/build/root_to_petsird --root-prefix ROOT_DATA/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 0 --petsird-file petsird_ETSIPETscanner256mmAFOV_voxBrain_0min_v0.2.bin -s data/root/ETSIPETscanner256mmAFOV_geometry.json

@run_3min_voxBrain: build
    cpp/build/root_to_petsird --root-prefix ROOT_DATA/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 18 --petsird-file petsird_ETSIPETscanner256mmAFOV_voxBrain_3min_v0.2.bin -s data/root/ETSIPETscanner256mmAFOV_geometry.json

@run_6min_voxBrain: build
    cpp/build/root_to_petsird --root-prefix ROOT_DATA/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 36 --petsird-file petsird_ETSIPETscanner256mmAFOV_voxBrain_6min_v0.2.bin -s data/root/ETSIPETscanner256mmAFOV_geometry.json

@run_ONLY_0min_mIEC:
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 0 --petsird-file petsird_ETSIPETscanner_Norm_mIEC_0min_v0.7.2.bin -s data/root/ETSIPETscanner_geometry.json

@run_ONLY_10sec_mIEC:
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 1 --petsird-file petsird_ETSIPETscanner_Norm_mIEC_10sec_v0.7.2.bin -s data/root/ETSIPETscanner_geometry.json

@run_ONLY_6min_mIEC:
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 36 --petsird-file petsird_ETSIPETscanner_Norm_mIEC_6min_v0.7.2.bin -s data/root/ETSIPETscanner_geometry.json



@run_ONLY_0min_voxBrain:
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 0 --petsird-file petsird_ETSIPETscanner256mmAFOV_Norm_voxBrain_0min_v0.7.2.bin -s data/root/ETS$

@run_ONLY_3min_voxBrain:
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 18 --petsird-file petsird_ETSIPETscanner256mmAFOV_Norm_voxBrain_3min_v0.7.2.bin -s data/root/ET$

@run_ONLY_6min_voxBrain:
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 36 --petsird-file petsird_ETSIPETscanner256mmAFOV_Norm_voxBrain_6min_v0.7.2.bin -s data/root/ET$



#Run only (built assumed validated)
#@run: run_ONLY_0min_mIEC
#@run: run_ONLY_10sec_mIEC
@run: run_ONLY_6min_mIEC

#@run: run_ONLY_0min_voxBrain
#@run: run_ONLY_3min_voxBrain
#@run: run_ONLY_6min_voxBrain



# Build and run
#@run: run_0min_mIEC
#@run: run_10sec_mIEC
#@run: run_6min_mIEC

#@run: run_0min_voxBrain
#@run: run_3min_voxBrain
#@run: run_6min_voxBrain
