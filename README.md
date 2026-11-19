# GATE ROOT to PETSIRD converter

## Background
The [Emission Tomography Standardization Initiative (ETSI)](https://etsinitiative.org/)
is working towards establishing a standard for PET Raw Data, called PETSIRD ("PET ETSI Raw Data").

This repository implements a convertor from output of [GATE](https://opengate.readthedocs.io/) (a Monte Carlo simulator) to PETSIRD

## Version information

Converter **1.0.0**, tested with **PETSIRD 0.13.3**. Both implementations support
ordered parallel conversion, automatic clock-reset repair, acquisition windows,
prompts and delayed coincidences, explicit unit detector efficiencies and embedded geometric normalization.

- [Execution modes, timestamp rules and CLI](docs/PARALLEL_CONVERSION.md)
- [Windows Python / WSL C++ installation and three complete EMITRON examples](docs/WORKSTATION.md)

## Installation instructions
1. Get the software
   ```sh
   git clone --branch v1.0.0 https://github.com/ETSIhackers/GATEROOT.git
   ```
2. Install dependencies, easiest via conda
   ```sh
   conda env create -f environment.yml
   conda activate GATEROOT
   ```
3. Build
   ```sh
   bash examples/install_cpp_sdk.sh
   just build
   ```

### Updating from a previous version

PETSIRD is not a submodule. The environment supplies native dependencies; the SDK
installation script generates and builds the pinned PETSIRD 0.13.3 release so an
older conda package cannot silently determine the binary interface.

- update the code as follows
  ```sh
  git pull
  ```
- update your conda environment

- remove existing build (as it might point to wrong environment) and build from scratch
  ```sh
   rm -rf cpp/build
   just build
   ```

## Example JSON file specifying geometry
```
{
  "model_name": "ETSIPETscanner1",
  "ax_phys_crystal_num": 5,
  "ax_virtual_crystal_num": 0,
  "n_cry_layers": 1,
  "n_cry_xy": 5,
  "n_cry_z": 5,
  "n_det": 600,
  "n_rings": 40,
  "n_rsec_xy": 60,
  "n_rsec_z": 1,
  "n_mod_xy": 2,
  "n_mod_z": 8,
  "n_smod_xy": 1,
  "n_smod_z": 1,
  "cry_ax_gap": 0,
  "cry_tx_gap": 0,
  "smod_ax_gap": 0,
  "smod_tx_gap": 0,
  "mod_ax_gap": 0,
  "mod_tx_gap": 0,
  "rsec_ax_gap": 0,
  "rsec_tx_gap": 0,
  "number_of_energy_bins": 1,
  "number_of_TOF_bins": 71,
  "TOF_bin_width_mm": 12.75,
  "radius": 308.0,
  "s_width": 380,
  "tx_phys_crystal_num": 10,
  "tx_virtual_crystal_num": 0,
  "detector_x_dim": 20.0,
  "detector_y_dim": 3.2,
  "detector_z_dim": 3.2,
  "energy_LLD": 430.0,
  "energy_ULD": 650.0,
  "EnergyResolutionAt511": 0.11,
  "TOF_resolution_mm": 30,
  "LM_time_block_duration_ms": 1
}
```
Units:
- dimensions and radius are in mm
- TOF resolution and bin width are in mm (to convert from ps, multiply with 0.29979245 / 2)
- LLD and ULD are in keV
- energy resolution is a fraction of 511 keV
- LM_time_block_duration_ms is in ms

### GATE geometry with more than one axial rsector (`n_rsec_z > 1`)

GATE writes `rsectorID` as `axial * n_rsec_xy + angular` (the angular index varies fastest), and the
converter relies on this. GATE only follows this convention if the axial copies of the rsector are made
with a `linear` repeater, or with an axial repeater placed before the `ring` repeater.
Do not use `ring` followed by `cubicArray`: GATE then writes `rsectorID` in a different order, and its own
coincidence sorter computes wrong sector differences, so the ROOT output itself is not correct.
