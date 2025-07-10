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



