# Standalone Windows Python converter

**GATEROOT 1.0.0 / PETSIRD 0.13.3.** For the current execution modes and clock/window
defaults see [the CLI guide](../docs/PARALLEL_CONVERSION.md). For the complete
Windows/WSL prompts+delayeds conversion and EMITRON scatter-corrected reconstruction
workflow see [the workstation guide](../docs/WORKSTATION.md).

This entry point uses the official PETSIRD 0.13 Python SDK (0.13.3 or later)
and uproot. It accepts the cylinder JSON and numbered ROOT prefixes used by the
C++ justfile.
The supported scanner profile is the one used by ETSIPETscanner1 and
ETSIPETscanner2, also with more than one ring of rsectors along z (`n_rsec_z`,
see the main README for the GATE repeaters): one DOI layer, one coincidence
energy bin, and no gaps or virtual crystals. It does not require CERN ROOT.

```powershell
conda create -n gateroot-100 -y python=3.12 pip
conda activate gateroot-100
python -m pip install -r .\python\requirements.txt
python .\python\root_to_petsird.py --help
```

`--selection all` retains trues, scatter and prompt randoms that pass the
geometry's energy and TOF windows. Set `--selection true` for the C++
example's same-event, no-phantom-Compton selection. Supply
`--delayed-offset-ps 500000` for the matching `delay` tree in these datasets.
The offset is added to `time1-time2` before endpoint canonicalization.

The default `--time-mode auto --file-duration-ms 10000` repairs the IEC's reset
10-second clocks and retains continuous Hoffman clocks. Explicit `per-file`
and `absolute` modes remain available. Time-window filtering occurs after repair.
The default execution mode is GPU when CuPy/CUDA is usable, otherwise CPU readers.
PETSIRD events use ascending time blocks and decreasing
detector IDs within each event, with the TOF direction changed on a swap.

`-c` takes the CASToR normalization Cdf, with its matching Cdh alongside
it. Files must be fully copied. The converter reads all NCF records even in
a limited ROOT test. It inverts positive NCFs, remaps ring-major CASToR IDs
to the GATE hierarchy, and averages inverse factors into rotation-symmetry
matrices. Zero values are retained in the pooling, and cells without any
input records have zero efficiency.

```powershell
$Base = 'F:\DATA\WCMC\ETSIhackathon'
$Out = Join-Path $Base 'GATEROOT_python_0132'
$Converter = '.\python\root_to_petsird.py'

python $Converter `
    --root-prefix (Join-Path $Base 'ROOT_files\IEC\ETSIPETscanner_mIEC_') `
    --number-of-root-files 2 --time-mode per-file --file-duration-ms 10000 `
    --petsird-file (Join-Path $Out 'IEC_smoke.petsird') `
    -s (Join-Path $Base 'GEOM_DATA\ETSIPETscanner_geometry.json') `
    -c (Join-Path $Base 'NORM_files\Sino_NCFs_ETSIPETscanner_span1_mrd39_df.Cdf') `
    --delayed-offset-ps 500000 --max-events-per-tree 100000 --verify

python $Converter `
    --root-prefix (Join-Path $Base 'ROOT_files\voxBrain\ETSIPETscanner256mmAFOV_positronSource\ETSIPETscanner256mmAFOV_voxBrain_') `
    --number-of-root-files 2 --time-mode absolute `
    --petsird-file (Join-Path $Out 'Hoffman_smoke.petsird') `
    -s (Join-Path $Base 'GEOM_DATA\ETSIPETscanner256mmAFOV_geometry.json') `
    -c (Join-Path $Base 'NORM_files\Sino_NCFs_ETSIPETscanner256mmAFOV_span1_mrd79_df.Cdf') `
    --delayed-offset-ps 500000 --max-events-per-tree 100000 --verify
```

Each output has a `.petsird.json` audit with event totals per file, normalization
record and SGID counts, and the SDK read-back result. Remove
`--max-events-per-tree 100000` and choose new filenames for full conversions.
Use 36 files for either six-minute acquisition. For the three-minute Hoffman
window test use all 36 files with `--end-time-s 180`. A file count of zero writes a header-only PETSIRD for
geometry and calibration testing.
