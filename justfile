set shell := ['bash', '-ceuo', 'pipefail']
version := "0.8.0"

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


#@run: build download-test-data
#    cpp/build/root_to_petsird --root-prefix data/root/ETSIPETscanner_mIEC_ --number-of-root-files 1 --petsird-file test.petsird -s data/root/scanner_geometry.json

#@run-full: build download-test-data-full
#    cpp/build/root_to_petsird --root-prefix data/root/ETSIPETscanner_mIEC_ --number-of-root-files 36 --petsird-file test-full.petsird -s data/root/scanner_geometry.json



@run_0min_mIEC: build
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 0 --petsird-file petsird_ETSIPETscanner_mIEC_0min_v{{version}}.petsird -s data/root/ETSIPETscanner_geometry.json

@run_10sec_mIEC: build
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 1 --petsird-file petsird_ETSIPETscanner_mIEC_10sec_v{{version}}.petsird -s data/root/ETSIPETscanner_geometry.json

@run_6min_mIEC: build
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 36 --petsird-file petsird_ETSIPETscanner_mIEC_6min_v{{version}}.petsird -s data/root/ETSIPETscanner_geometry.json


@run_0min_voxBrain: build
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 0 --petsird-file petsird_ETSIPETscanner256mmAFOV_voxBrain_0min_v{{version}}.petsird -s data/root/ETSIPETscanner256mmAFOV_geometry.json

@run_3min_voxBrain: build
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 18 --petsird-file petsird_ETSIPETscanner256mmAFOV_voxBrain_3min_v{{version}}.petsird -s data/root/ETSIPETscanner256mmAFOV_geometry.json

@run_6min_voxBrain: build
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 36 --petsird-file petsird_ETSIPETscanner256mmAFOV_voxBrain_6min_v{{version}}.petsird -s data/root/ETSIPETscanner256mmAFOV_geometry.json


@run_ONLY_0min_mIEC: 
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 0 --petsird-file petsird_ETSIPETscanner_realNorm_mIEC_0min_v{{version}}.petsird -s data/root/ETSIPETscanner_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner_span1_mrd39_df.Cdf

@run_ONLY_10sec_mIEC:
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 1 --petsird-file petsird_ETSIPETscanner_realNorm_mIEC_10sec_v{{version}}.petsird -s data/root/ETSIPETscanner_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner_span1_mrd39_df.Cdf

@run_ONLY_6min_mIEC:
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 36 --petsird-file petsird_ETSIPETscanner_realNorm_mIEC_6min_v{{version}}.petsird -s data/root/ETSIPETscanner_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner_span1_mrd39_df.Cdf



@run_ONLY_0min_voxBrain:
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 0 --petsird-file petsird_ETSIPETscanner256mmAFOV_realNorm_voxBrain_0min_v{{version}}.petsird -s data/root/ETSIPETscanner256mmAFOV_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner256mmAFOV_span1_mrd79_df.Cdf

@run_ONLY_3min_voxBrain:
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 18 --petsird-file petsird_ETSIPETscanner256mmAFOV_realNorm_voxBrain_3min_v{{version}}.petsird -s data/root/ETSIPETscanner256mmAFOV_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner256mmAFOV_span1_mrd79_df.Cdf

@run_ONLY_6min_voxBrain:
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 36 --petsird-file petsird_ETSIPETscanner256mmAFOV_realNorm_voxBrain_6min_v{{version}}.petsird -s data/root/ETSIPETscanner256mmAFOV_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner256mmAFOV_span1_mrd79_df.Cdf



@run_ONLY_ALLmin_mIEC:
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 0 --petsird-file petsird_ETSIPETscanner_realNorm_mIEC_0min_v{{version}}.petsird -s data/root/ETSIPETscanner_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner_span1_mrd39_df.Cdf
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 1 --petsird-file petsird_ETSIPETscanner_realNorm_mIEC_10sec_v{{version}}.petsird -s data/root/ETSIPETscanner_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner_span1_mrd39_df.Cdf
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 36 --petsird-file petsird_ETSIPETscanner_realNorm_mIEC_6min_v{{version}}.petsird -s data/root/ETSIPETscanner_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner_span1_mrd39_df.Cdf


@run_ONLY_ALLmin_voxBrain:
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 0 --petsird-file petsird_ETSIPETscanner256mmAFOV_realNorm_voxBrain_0min_v{{version}}.petsird -s data/root/ETSIPETscanner256mmAFOV_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner256mmAFOV_span1_mrd79_df.Cdf
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 18 --petsird-file petsird_ETSIPETscanner256mmAFOV_realNorm_voxBrain_3min_v{{version}}.petsird -s data/root/ETSIPETscanner256mmAFOV_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner256mmAFOV_span1_mrd79_df.Cdf
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 36 --petsird-file petsird_ETSIPETscanner256mmAFOV_realNorm_voxBrain_6min_v{{version}}.petsird -s data/root/ETSIPETscanner256mmAFOV_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner256mmAFOV_span1_mrd79_df.Cdf


@run_ONLY_ALLmin_ALLphantoms:
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 0 --petsird-file petsird_ETSIPETscanner_realNorm_mIEC_0min_v{{version}}.petsird -s data/root/ETSIPETscanner_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner_span1_mrd39_df.Cdf
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 1 --petsird-file petsird_ETSIPETscanner_realNorm_mIEC_10sec_v{{version}}.petsird -s data/root/ETSIPETscanner_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner_span1_mrd39_df.Cdf
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/IEC/ETSIPETscanner_mIEC_ --number-of-root-files 36 --petsird-file petsird_ETSIPETscanner_realNorm_mIEC_6min_v{{version}}.petsird -s data/root/ETSIPETscanner_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner_span1_mrd39_df.Cdf
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 0 --petsird-file petsird_ETSIPETscanner256mmAFOV_realNorm_voxBrain_0min_v{{version}}.petsird -s data/root/ETSIPETscanner256mmAFOV_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner256mmAFOV_span1_mrd79_df.Cdf
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 18 --petsird-file petsird_ETSIPETscanner256mmAFOV_realNorm_voxBrain_3min_v{{version}}.petsird -s data/root/ETSIPETscanner256mmAFOV_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner256mmAFOV_span1_mrd79_df.Cdf
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/voxBrain/ETSIPETscanner256mmAFOV_positronSource/ETSIPETscanner256mmAFOV_voxBrain_ --number-of-root-files 36 --petsird-file petsird_ETSIPETscanner256mmAFOV_realNorm_voxBrain_6min_v{{version}}.petsird -s data/root/ETSIPETscanner256mmAFOV_geometry.json -c ../NORM_files/Sino_NCFs_ETSIPETscanner256mmAFOV_span1_mrd79_df.Cdf

@run_ONLY_Norm_ETSIPETscanner1:
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/Norm/ETSIPETscanner_Norm_ --number-of-root-files 72 --petsird-file petsird_ETSIPETscanner_NormScan_18360sec_v{{version}}.petsird -s data/root/ETSIPETscanner_geometry.json

@run_ONLY_Norm_ETSIPETscanner2:
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/Norm/ETSIPETscanner256mmAFOV_Norm_ --number-of-root-files 36 --petsird-file petsird_ETSIPETscanner256mmAFOV_NormScan_18000sec_v{{version}}.petsird -s data/root/ETSIPETscanner256mmAFOV_geometry.json

@run_ONLY_Norm_ALLscanners:
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/Norm/ETSIPETscanner_Norm_ --number-of-root-files 72 --petsird-file petsird_ETSIPETscanner_NormScan_18360sec_v{{version}}.petsird -s data/root/ETSIPETscanner_geometry.json
    cpp/build/root_to_petsird --root-prefix ../ROOT_files/Norm/ETSIPETscanner256mmAFOV_Norm_ --number-of-root-files 36 --petsird-file petsird_ETSIPETscanner256mmAFOV_NormScan_18000sec_v{{version}}.petsird -s data/root/ETSIPETscanner256mmAFOV_geometry.json




#Run only (built assumed validated)
#@run: run_ONLY_0min_mIEC
#@run: run_ONLY_10sec_mIEC
#@run: run_ONLY_6min_mIEC

#@run: run_ONLY_0min_voxBrain
#@run: run_ONLY_3min_voxBrain
#@run: run_ONLY_6min_voxBrain

#@run: run_ONLY_ALLmin_mIEC
#@run: run_ONLY_ALLmin_voxBrain
#@run: run_ONLY_ALLmin_ALLphantoms

#@run: run_ONLY_Norm_ETSIPETscanner1
@run: run_ONLY_Norm_ETSIPETscanner2
#run: run_ONLY_Norm_ALLscanners

# Build and run
#@run: run_0min_mIEC
#@run: run_10sec_mIEC
#@run: run_6min_mIEC

#@run: run_0min_voxBrain
#@run: run_3min_voxBrain
#@run: run_6min_voxBrain
