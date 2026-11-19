# IEC and Hoffman validation on Windows + WSL

These scripts create **all prompts** (trues + scatter + prompt randoms), store the
ROOT `delay` coincidences separately in PETSIRD, and embed the supplied
normalization. They use PETSIRD **0.13.3**. They do not modify previous conversions.

| Case | ROOT inputs | Corrected time window | Geometry |
|---|---|---|---|
| IEC_6min | All 36 IEC files | `[0,360)` seconds | ETSIPETscanner1 |
| Hoffman_3min | All 36 voxBrain files | `[0,180)` seconds | ETSIPETscanner2, 256 mm AFOV |
| Hoffman_6min | All 36 voxBrain files | `[0,360)` seconds | ETSIPETscanner2, 256 mm AFOV |

The 3-minute test deliberately specifies all 36 files, not just 18. The timing
prepass sees them all, repairs reset clocks if needed, and then applies the window.
Both formats use the same geometry JSON / expanded normalization Cdf+Cdh as your
existing justfile. The converter pools that expanded normalization into the same
rotation-symmetry representation as before; it is not a new normalization model.

## 1. Windows PowerShell: install the Python converter

Use a fresh checkout/environment to preserve ongoing EMITRON runs. Paste only the
commands, not terminal prompts. `--execution auto` selects CuPy/CUDA when usable;
explicit `gpu` below makes a missing GPU fail instead of quietly changing modes.
Your RTX 5090 needs a current NVIDIA Windows driver and a Blackwell-compatible
CuPy/CUDA runtime. ROOT decoding itself runs on CPU reader threads.

```powershell
conda create -n gateroot-100 -y python=3.12 pip
conda activate gateroot-100
$Base = 'F:\DATA\WCMC\ETSIhackathon'
$Repo = Join-Path $Base 'GATEROOT-v1.0.0'
git clone --branch v1.0.0 https://github.com/ETSIhackers/GATEROOT.git "$Repo"
if ($LASTEXITCODE -ne 0) { throw 'Clone failed' }
Set-Location -LiteralPath $Repo
python -m pip install -r .\python\requirements.txt
if ($LASTEXITCODE -ne 0) { throw 'Python dependencies failed' }
python -m pip install 'cupy-cuda12x>=13.4' 'nvidia-cuda-runtime-cu12>=12.8' 'nvidia-cuda-nvrtc-cu12>=12.8' pytest
if ($LASTEXITCODE -ne 0) { throw 'CUDA Python dependencies failed' }
python -c "import cupy as cp; import petsird; from importlib.metadata import version; print('PETSIRD',version('petsird')); print('GPUs',cp.cuda.runtime.getDeviceCount()); print(cp.arange(10).sum().get())"
if ($LASTEXITCODE -ne 0) { throw 'GPU preflight failed: check NVIDIA driver/CuPy runtime' }
python -m pytest -v .\test\test_parallel_conversion.py
if ($LASTEXITCODE -ne 0) { throw 'Python converter tests failed' }
```

The C++ tests skip in this Windows-only environment. A C++/CUDA runtime test runs
in step 3 under WSL. CuPy's CUDA_PATH advisory alone is not a failure if the device
operation succeeds. If your existing driver only supports a different CUDA major,
install the corresponding CuPy package instead; never install both CUDA wheel variants.

## 2. Convert the three datasets with Python

```powershell
$PythonOut = Join-Path $Base 'GATEROOT_100_python_prompts_delayeds'
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
& .\examples\convert_workstation.ps1 -Root $Base -Destination $PythonOut -Execution gpu -Workers 4
if ($LASTEXITCODE -ne 0) { throw 'Python conversion failed' }
Get-ChildItem -LiteralPath $PythonOut -Filter '*.petsird*'
```

Use a new destination if you have already run this exact example. To benchmark CPU
or serial mode use another destination and `-Execution cpu` / `-Execution none`.
Do not compare wall times while other conversions/reconstructions compete for disk/GPU.

The supplied IEC delay tree has approximately 500 ns shifted `time2`; `500000` ps
recenters its TOF difference. The scripts set this explicitly for both scanners.
Verify it if your GATE sorter configuration changes. `delayed_randoms_scale=1` in
reconstruction assumes equal prompt/delayed coincidence-window widths, not equal
prompt and delayed event counts. No randoms percentage is fabricated.

## 3. WSL Ubuntu terminal: build and validate C++ CUDA

These are **Ubuntu commands**, not PowerShell commands. From PowerShell use
`wsl -d Ubuntu` once to enter Ubuntu; do not run `wsl` inside Ubuntu.
The NVIDIA Windows driver supplies WSL GPU access: `nvidia-smi` must work in WSL.
The following assumes conda/miniforge is already initialized in Ubuntu, as in
your previous C++ installation. Installing the CUDA toolkit here does not install
a Linux display driver.

```bash
cd "$HOME"
git clone --branch v1.0.0 https://github.com/ETSIhackers/GATEROOT.git GATEROOT-v1.0.0
cd GATEROOT-v1.0.0
conda env create -f environment.yml
conda activate GATEROOT
bash examples/install_cpp_sdk.sh
conda install -y -c nvidia cuda-nvcc=12.8 cuda-cudart-dev=12.8 cuda-cccl=12.8
sudo apt-get update
sudo apt-get install -y g++-13
nvidia-smi
cmake -S cpp -B cpp/build -G Ninja -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" -DGATEROOT_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=120 -DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++-13
cmake --build cpp/build --parallel 2
GATEROOT_TEST_CPP_GPU=1 python -m pytest -v test
```

CUDA 12.8 and architecture `120` target the RTX 5090. On an older Ubuntu release
without g++-13, use a compiler supported by CUDA 12.8 (for example g++-12), changing
both package name and `CMAKE_CUDA_HOST_COMPILER`. Do not use conda GCC 15 with
CUDA 12.8. A CPU-only build uses `-DGATEROOT_CUDA=OFF` and needs neither CUDA nor
the host-compiler/architecture arguments. Do not run the C++ GPU parity test with
a CPU-only build.

## 4. Convert with C++ in the Ubuntu terminal

```bash
bash examples/convert_workstation.sh /mnt/f/DATA/WCMC/ETSIhackathon /mnt/f/DATA/WCMC/ETSIhackathon/GATEROOT_100_cpp_prompts_delayeds gpu 4
```

You can use separate terminals for Python/C++, with distinct destinations. On the
external hard disk, run them sequentially for useful timing comparisons and less
I/O contention. The final single writer preserves chronology in every mode.
Reports beside the outputs record backend, counts, clock shifts and wall times.

## 5. EMITRON preliminary reconstruction, scatter and final reconstruction

Use your existing **EMITRON development environment**, not `gateroot-100` or the
WSL converter environment. The source checkout must include sampled SSS,
`activity_estimation`, optimized exact sensitivity (`collapse_tof_bins`),
precorrected-additive PETSIRD sidecars and persistent decoded caches. This is the
`EMITRON-gate-root-petsird013` checkout used in your recent tests.

The example uses that environment's supported PETSIRD reader. The 0.13.2 reader
in that EMITRON revision can decode the schema-compatible 0.13.3 converter output;
the planning command checks this before scanning events. Do not upgrade the SDK
inside an EMITRON revision that explicitly pins its reader to 0.13.2. The converter
itself is installed/tested with 0.13.3 in its separate environment.

The two supplied original STIR headers are read from:

```
F:\DATA\WCMC\ETSIhackathon\GATEROOT_python_0132_trues\umap_mIEC_ETSIPETscanner.hv
F:\DATA\WCMC\ETSIhackathon\GATEROOT_python_0132_trues\umap_ETSIPETscanner256mmAFOV_brainSimplePhantom.hv
```

Their paired `.v` payloads must be present. The example centers each attenuation
grid at `(0,0,0)`, matching the centered maps used in your earlier GATE reconstructions;
it preserves voxel values/spacing without interpolation or axis flips. This is a
registration assumption for these particular phantoms: inspect the saved centered
umaps before using it for other data. Explicitly select their actual units. The
commands below assume **cm^-1**; for a map in mm^-1 use `--mu-units 'mm^-1'`, which
multiplies values by 10 before declaring cm^-1. Do not guess units from filename.

For each case the runner executes:

1. Extract PETSIRD delayed events into a per-detector-pair, non-TOF random-count histogram.
2. Run **2 OSEM iterations / 10 subsets**, attenuation and embedded normalization,
   with delayed randoms but no scatter correction. Generate optimized exact TOF
   sensitivity once and save it. No true activity distribution is provided.
3. Resample that preliminary image to the attenuation grid and compute a sampled
   single-scatter distribution. Fit its amplitude in attenuation-defined tails
   after subtracting estimated randoms.
4. Compile `q=(scatter_rate+random_rate)/multiplicative_factor` into the ordered
   additive-only PETSIRD correction sidecar. Keep normalization Ã— attenuation in
   the independent multiplicative model and sensitivity.
5. Start a fresh **2-iteration / 10-subset** final reconstruction from the configured
   initial image, reusing the preliminary sensitivity, and save only iteration 2.

The final reconstruction is two complete iterations, not continuation from
iteration 2 of the preliminary reconstruction. Scatter/random additive changes
do not change sensitivity. Each implementation/case has separate output folders.
The script currently generates one preliminary sensitivity per case, then reuses
it for that case's final reconstruction; it does not silently reuse another
case's differently fingerprinted normalization/attenuation configuration.

Reconstruction grid: `360Ã—360Ã—64` (IEC) / `360Ã—360Ã—128` (Hoffman), 2-mm isotropic,
centered on the scanner. TOF resolution and timing edges come from the PETSIRD
header. No arbitrary PSF or post-filter is added. Physical frame durations are
360/180/360 seconds; no unprovided isotope or absolute Bq/ml calibration is assumed.

Randoms are pooled across TOF bins **per detector pair**, then distributed uniformly
across the coincidence timing window. This is explicit variance reduction, not
an assertion that randoms are uniform across detector pairs. The SSS example uses
16 angular detectors Ã— 8 rings, `[5,5,1]` umap coarsening, and tail fitting. This is
EMITRON's approximate single-Compton model, not a validated STIR/openSSS implementation
or a model of multiple scattering. Validate sampling convergence and tail-fit
residuals before drawing quantitative conclusions. It retains the full accepted
pair domain for tail exposure and sensitivity, rather than only observed LORs.

In a **Windows PowerShell** terminal:

```powershell
conda activate emitron-optimized-tof
$Base = 'F:\DATA\WCMC\ETSIhackathon'
$ConverterRepo = Join-Path $Base 'GATEROOT-v1.0.0'
$EmitronRepo = 'F:\DATA\WCMC\UHBPET\EMITRON\EMITRON-gate-root-petsird013'
Set-Location -LiteralPath $EmitronRepo
$Runner = Join-Path $ConverterRepo 'examples\reconstruct_emitron.py'
$env:PYTHONPATH = "$EmitronRepo;$env:PYTHONPATH"
$MuDir = Join-Path $Base 'GATEROOT_python_0132_trues'
$Cache = Join-Path $env:LOCALAPPDATA 'EMITRON\gate100-cache'
$CppOut = Join-Path $Base 'GATEROOT_100_cpp_prompts_delayeds'
$PythonOut = Join-Path $Base 'GATEROOT_100_python_prompts_delayeds'
$CppRecon = Join-Path $CppOut 'EMITRON_SSS_2i10s'
$PythonRecon = Join-Path $PythonOut 'EMITRON_SSS_2i10s'
$Common = @('--mu-dir',$MuDir,'--mu-units','cm^-1','--cache-dir',$Cache)
python $Runner --input-dir $CppOut --output $CppRecon @Common --plan
if ($LASTEXITCODE -ne 0) { throw 'C++ data/header/umap planning failed' }
python $Runner --input-dir $PythonOut --output $PythonRecon @Common --plan
if ($LASTEXITCODE -ne 0) { throw 'Python data/header/umap planning failed' }
```

Inspect `plan.json` and `mu511_centered.hdr` in each case folder. Then run the six
pipelines **sequentially** to avoid the memory exhaustion seen with concurrent jobs:

```powershell
python $Runner --input-dir $CppOut --output $CppRecon @Common
if ($LASTEXITCODE -ne 0) { throw 'C++ PETSIRD reconstruction batch failed' }
python $Runner --input-dir $PythonOut --output $PythonRecon @Common
if ($LASTEXITCODE -ne 0) { throw 'Python PETSIRD reconstruction batch failed' }
```

To run only IEC initially append `--cases IEC_6min`. Later use
`--cases Hoffman_3min Hoffman_6min`. A completed `image_iter002.npy` is retained;
completed scatter estimates and sidecars are reused. Inputs are fingerprinted;
use a new output directory when changing inputs or settings. Incomplete sidecar
creation is not automatically deleted: preserve logs and choose a new folder
if compilation was interrupted. Keep substantial free space in the cache drive
and output drive (raw event caches, delayed histograms and sidecars can be many GB).

Final images are `EMITRON_SSS_2i10s/<CASE>/image_iter002.npy` plus Interfile outputs.
Preliminary activity, sensitivity, sampled scatter/tail-fit audit, randoms
descriptor, sidecar, final JSON and input registration audit remain alongside.
The large six-case SSS/reconstruction workflow is an executable local validation
example; it is not claimed to have been run to completion in CI.
