# Test data: GATE simulation of the PTSRC point sources

Small GATE simulation used by `test/test_gate_point_sources.py`, with the 4 point-source positions of the
PTSRC protocol in [PETSIRDTestData](https://github.com/ETSInitiative/PETSIRDTestData/blob/main/doc/AcquisitionProtocols.md)
(L = 256 mm, rho = 120 mm). Every coincidence contains the true emission position (`sourcePosX/Y/Z`) and
the source number (`sourceID`), so the converted PETSIRD file can be checked against them.

Unlike the protocol, the 4 sources are in the same acquisition; they can be separated with `sourceID`.

- `ptsrc_1.root`: GATE output (`Coincidences` tree only), 23,533 coincidences of which 23,461 trues
- `main.mac`: GATE macro used to generate it
- `points.txt`: true positions of the 4 point sources (mm, GATE world coordinates), in `sourceID` order
- `geometry.json`: scanner description for `root_to_petsird`

## Simulation

- GATE 9.2 (Docker image `opengatecollaboration/gate:9.2-docker`), physics list `emstandard_opt4`, 1 mm cuts.
- Scanner: ETSI PET scanner with 256 mm axial FOV (as in Zenodo record 15930738): 60 rsectors,
  modules 1 x 2 x 16, 1 submodule, crystals 1 x 5 x 5 of 20 x 3.2 x 3.2 mm LYSO, inner radius 308 mm.
- Sources: 4 back-to-back 511 keV gamma point sources (no positron range, no acollinearity), 90 kBq each, in air.
- Digitizer: adder (every crystal read out on its own), energy blurring 11 % FWHM at 511 keV,
  energy window 430-650 keV, 141 ps FWHM per single (200 ps coincidence timing resolution),
  coincidence window 4.3 ns, multiples policy `takeWinnerOfGoods`.
- Acquisition: 2 s, random seed 20261116.

To regenerate: put `main.mac` in `mac/`, create an empty `output/` directory, download `GateMaterials.db` from
https://raw.githubusercontent.com/OpenGATE/Gate/v9.2/GateMaterials.db and run

```sh
docker run --rm -v "$PWD":/APP opengatecollaboration/gate:9.2-docker mac/main.mac
```
