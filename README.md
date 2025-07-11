# GATE ROOT to PETSIRD converter

## Background
The [Emission Tomography Standardization Initiative (ETSI)](https://etsinitiative.org/)
is working towards establishing a standard for PET Raw Data, called PETSIRD ("PET ETSI Raw Data").

This repository implements a convertor from output of [GATE](https://opengate.readthedocs.io/) (a Monte Carlo simulator) to PETSIRD

## Version information

This compiles with PETSIRD v0.7.2.

## Installation instructions
1. Get the software
   ```sh
   git clone --recurse-submodules https://github.com/ETSIhackers/GATEROOT.git
   ```
2. Install dependencies, easiest via conda
   ```sh
   conda env create -f environment.yml
   conda activate GATEROOT
   ```
3. Install yardl 0.6.3 and add it to your path

4. Build
   ```sh
   just build
   ```

### Updating from a previous version

As we use a PETSIRD submodule, you have to take this into account

- update the code as follows
  ```sh
  git pull
  git submodule update --init
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
  "n_cry_layers": 1,
  "n_cry_xy": 5,
  "n_cry_z": 5,
  "n_mod_xy": 2,
  "n_mod_z": 8,
  "n_rsec_xy": 60,
  "n_rsec_z": 1,
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
  "radius": 308.0,
  "detector_x_dim": 20.0,
  "detector_y_dim": 3.2,
  "detector_z_dim": 3.2,
  "number_of_TOF_bins": 62,
  "TOF_bin_width_mm": 1.1,
  "TOF_resolution_mm": 3.3,
  "number_of_energy_bins": 1,
  "energy_LLD": 430.0,
  "energy_ULD": 650.0,
  "EnergyResolutionAt511": 0.11,
  "LM_time_block_duration_ms": 1
}
```
Units:
- dimensions and radius are in mm
- TOF resolution and bin width are in mm (to convert from ps, multiply with 0.29979245 / 2)
- LLD and ULD are in keV
- energy resolution is a fraction of 511 keV
- LM_time_block_duration_ms is in ms

