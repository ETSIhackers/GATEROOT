"""Check root_to_petsird against a small GATE simulation of the 4 PTSRC point sources (see data/README.md).

The converted file is compared with the GATE ground truth (true source position of every coincidence)
and with the scanner description in the JSON file. Positions of the crystals come only from the
geometry in the PETSIRD header (petsird.helpers.geometry), not from the converter code.

The converter writes the true coincidences in the order of the ROOT file, so the k-th PETSIRD event
corresponds to the k-th true coincidence in the ROOT file.

Run from the repository root after building:  python -m pytest test
(set ROOT_TO_PETSIRD to use another executable than cpp/build/root_to_petsird)
"""
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import petsird
import pytest
import ROOT
from petsird.helpers import expand_detection_bin
from petsird.helpers.geometry import get_ADP

DATA = Path(__file__).parent / "data"
CONVERTER = Path(os.environ.get("ROOT_TO_PETSIRD", Path(__file__).parents[1] / "cpp" / "build" / "root_to_petsird"))
TYPE_OF_MODULE = petsird.TypeOfModule(0)


@pytest.fixture(scope="module")
def converted(tmp_path_factory):
    """header and prompt events of the converted test file"""
    out = tmp_path_factory.mktemp("petsird") / "ptsrc.petsird"
    subprocess.run([str(CONVERTER), "-r", str(DATA / "ptsrc_"), "-n", "1",
                    "-s", str(DATA / "geometry.json"), "-p", str(out)], check=True, capture_output=True)
    with petsird.BinaryPETSIRDReader(str(out)) as reader:
        header = reader.read_header()
        events = [event for time_block in reader.read_time_blocks()
                  if isinstance(time_block, petsird.TimeBlock.EventTimeBlock)
                  for event in time_block.value.prompt_events[TYPE_OF_MODULE][TYPE_OF_MODULE]]
    return header, events


@pytest.fixture(scope="module")
def true_sources():
    """true emission position and source number (GATE sourceID) of every true coincidence in the ROOT file"""
    df = ROOT.RDataFrame("Coincidences", str(DATA / "ptsrc_1.root"))
    df = df.Filter("eventID1 == eventID2 && comptonPhantom1 == 0 && comptonPhantom2 == 0")
    cols = df.AsNumpy(["sourcePosX1", "sourcePosY1", "sourcePosZ1", "sourceID1"])
    positions = np.stack([cols["sourcePosX1"], cols["sourcePosY1"], cols["sourcePosZ1"]], axis=1).astype(float)
    return positions, cols["sourceID1"].astype(int)


@pytest.fixture(scope="module")
def lors(converted):
    """average detection points of detection_bins[0] (A) and detection_bins[1] (B) of every event"""
    header, events = converted
    cache = {}

    def adp(detection_bin):
        expanded = expand_detection_bin(header.scanner, TYPE_OF_MODULE, detection_bin)
        key = (expanded.module_index, expanded.element_index)
        if key not in cache:
            cache[key] = np.array(get_ADP(header.scanner, TYPE_OF_MODULE, expanded).c, dtype=float)
        return cache[key]

    A = np.array([adp(event.detection_bins[0]) for event in events])
    B = np.array([adp(event.detection_bins[1]) for event in events])
    return A, B


def test_number_of_events(converted, true_sources):
    assert len(converted[1]) == len(true_sources[0])


def test_lors_pass_through_sources(lors, true_sources):
    """the LOR of every event must pass close to its source (crystals are 3.2 mm wide)"""
    A, B = lors
    S, source_id = true_sources
    distance = np.linalg.norm(np.cross(S - A, B - A), axis=1) / np.linalg.norm(B - A, axis=1)
    points = np.loadtxt(DATA / "points.txt")  # in sourceID order
    for k, point in enumerate(points):
        median = np.median(distance[source_id == k])
        assert median < 2.0, f"source at {point} mm: median distance to the LORs is {median:.2f} mm"


def test_tof(converted, lors, true_sources):
    """PETSIRD stores (t1 - t2) * c / 2, with t1 the arrival time on detection_bins[0]"""
    header, events = converted
    A, B = lors
    S, _ = true_sources
    edges = np.asarray(header.scanner.tof_bin_edges[TYPE_OF_MODULE][TYPE_OF_MODULE].edges)
    centres = 0.5 * (edges[:-1] + edges[1:])
    measured = centres[np.array([event.tof_idx for event in events])]
    expected = 0.5 * (np.linalg.norm(S - A, axis=1) - np.linalg.norm(S - B, axis=1))
    difference = measured - expected
    # expected spread: TOF resolution (30 mm FWHM) and bin width (30 mm), about 15 mm
    assert abs(np.median(difference)) < 2.0
    assert np.std(difference) < 20.0


def test_crystal_size(converted):
    header, _ = converted
    geometry = json.loads((DATA / "geometry.json").read_text())
    box = header.scanner.scanner_geometry.replicated_modules[TYPE_OF_MODULE].object.detecting_elements.object.shape
    corners = np.array([corner.c for corner in box.corners])
    sizes = corners.max(axis=0) - corners.min(axis=0)
    expected = [geometry["detector_x_dim"], geometry["detector_y_dim"], geometry["detector_z_dim"]]
    np.testing.assert_allclose(sizes, expected, atol=1e-3)


def test_crystals_do_not_overlap(converted):
    header, _ = converted
    elements = header.scanner.scanner_geometry.replicated_modules[TYPE_OF_MODULE].object.detecting_elements
    corners = np.array([corner.c for corner in elements.object.shape.corners])
    sizes = corners.max(axis=0) - corners.min(axis=0)
    # inside a module the crystals are only translated
    shifts = np.array([[transform.matrix[i][3] for i in range(3)] for transform in elements.transforms])
    overlap = np.all(np.abs(shifts[:, None, :] - shifts[None, :, :]) < sizes - 1e-3, axis=2)
    np.fill_diagonal(overlap, False)
    assert not overlap.any(), f"{overlap.sum() // 2} pairs of crystals overlap"
