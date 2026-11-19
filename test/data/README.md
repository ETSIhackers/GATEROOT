# Test data: GATE simulation of the PTSRC point sources

Small GATE simulation used by `test/test_gate_point_sources.py`, with the 4 point-source positions of the
PTSRC protocol in [PETSIRDTestData](https://github.com/ETSInitiative/PETSIRDTestData/blob/main/doc/AcquisitionProtocols.md)
(L = 256 mm, rho = 120 mm). Every coincidence contains the true emission position (`sourcePosX/Y/Z`) and
the source number (`sourceID`), so the converted PETSIRD file can be checked against them.

The protocol asks for 1 source per acquisition, so there is one ROOT file per source (`Coincidences` tree only).
They come from one simulation of the 4 sources, split by `sourceID` (see below); the 53 coincidences between
two different sources (randoms) are not in any of them.

| file | source | position (mm) | coincidences | trues |
|---|---|---|---|---|
| `ptsrc_P1_1.root` | P1 | (0, -10, 115.2) | 1,211 | 1,208 |
| `ptsrc_P2_1.root` | P2 | (120, 0, 0) | 10,081 | 10,077 |
| `ptsrc_P3_1.root` | P3 | (10, 0, 0) | 10,786 | 10,782 |
| `ptsrc_P4_1.root` | P4 | (120, 0, 115.2) | 1,407 | 1,399 |

The scanner is made of 2 rings of rsectors along z (`n_rsec_z = 2`), so that the tests also check the order
of GATE's `rsectorID` (axial and angular index, see #42). The crystals are at the same places as with 1 ring of
rsectors of 16 modules along z.

- `main.mac`: GATE macro used for the simulation
- `points.txt`: true positions of the 4 point sources (mm, GATE world coordinates), in `sourceID` order (P1 to P4)
- `geometry.json`: scanner description for `root_to_petsird`

## Simulation

- GATE 9.2 (Docker image `opengatecollaboration/gate:9.2-docker`), physics list `emstandard_opt4`, 1 mm cuts.
- Scanner: ETSI PET scanner with 256 mm axial FOV (as in Zenodo record 15930738), but made of 2 rings of
  60 rsectors along z (`ring` then `linear` repeater, 128 mm apart), modules 1 x 2 x 8 per rsector,
  1 submodule, crystals 1 x 5 x 5 of 20 x 3.2 x 3.2 mm LYSO, inner radius 308 mm.
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

and split `output/ptsrc_1.root` by source with ROOT:

```python
import ROOT
simulation = ROOT.TFile.Open("output/ptsrc_1.root")
tree = simulation.Get("Coincidences")
for k, name in enumerate(["P1", "P2", "P3", "P4"]):
    out = ROOT.TFile(f"ptsrc_{name}_1.root", "RECREATE")
    tree.CopyTree(f"sourceID1 == {k} && sourceID2 == {k}").Write("", ROOT.TObject.kOverwrite)
    out.Close()
```
