"""Check the converters against a small GATE simulation of the 4 PTSRC point sources (see data/README.md).

As the PTSRC protocol asks for 1 source per acquisition, each source (P1 to P4) is in its own ROOT file.
The backprojection test runs on each of them. The other tests do not depend on the position of the source
and only use P4, which is away from the centre in x and in z (the TOF check needs a source away from the centre).

1. The C++ converter (cpp/build/root_to_petsird) and the Python converter (python/root_to_petsird.py)
   must write the same PETSIRD content: same events, time blocks and header (geometry within 1e-3 mm).
2. PTSRC test of the acquisition protocols in PETSIRDTestData (doc/AcquisitionProtocols.md):
   geometric backprojection of the whole file (all efficiencies 1), without and with TOF, compared with
   the true position of the source. The protocol asks for the centre of mass of the backprojection, but
   for a source away from the centre of the FOV the centre of mass is pulled by the long, asymmetric
   tails of the backprojection (without TOF, about half way to the centre of the FOV), so the position
   of the maximum is used instead.
3. Checks that the backprojection does not cover well: the TOF sign, event by event, and the crystal
   boxes (the backprojection only uses the average detection points).

The converted files are the ones of the C++ converter. Positions of the crystals come only from the
geometry in the PETSIRD header (petsird.helpers.geometry), not from the converter code.
The converter writes the true coincidences in the order of the ROOT file, so the k-th PETSIRD event
corresponds to the k-th true coincidence in the ROOT file (used by the TOF check).

Run from the repository root after building:  python -m pytest test
(set ROOT_TO_PETSIRD to use another executable than cpp/build/root_to_petsird)
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import petsird
import pytest
import ROOT
from parallelproj.projectors import ListmodePETProjector, TOFParameters
from petsird.helpers import expand_detection_bin
from petsird.helpers.geometry import get_ADP

DATA = Path(__file__).parent / "data"
REPOSITORY = Path(__file__).parents[1]
CONVERTER = Path(os.environ.get("ROOT_TO_PETSIRD", REPOSITORY / "cpp" / "build" / "root_to_petsird"))
PYTHON_CONVERTER = REPOSITORY / "python" / "root_to_petsird.py"
TYPE_OF_MODULE = petsird.TypeOfModule(0)

# PTSRC sources of the protocol, each in data/ptsrc_<name>_1.root; true positions in data/points.txt (same order)
SOURCES = ["P1", "P2", "P3", "P4"]
POINTS = np.loadtxt(DATA / "points.txt")

# image for the backprojection: 1 mm voxels, 340 x 340 x 280 mm, centred on the scanner
VOXEL_MM = 1.0
IMAGE_SHAPE = (340, 340, 280)
IMAGE_ORIGIN = np.array([-(n - 1) / 2 * VOXEL_MM for n in IMAGE_SHAPE])  # centre of the first voxel


def read_petsird(path):
    """header, events (detection bins and tof_idx) and time blocks (start, stop, number of events)"""
    with petsird.BinaryPETSIRDReader(str(path)) as reader:
        header = reader.read_header()
        events, blocks = [], []
        for time_block in reader.read_time_blocks():
            if isinstance(time_block, petsird.TimeBlock.EventTimeBlock):
                block = time_block.value
                block_events = block.prompt_events[TYPE_OF_MODULE][TYPE_OF_MODULE]
                events += block_events
                blocks.append((block.time_interval.start, block.time_interval.stop, len(block_events)))
    bins = np.array([event.detection_bins for event in events], dtype=np.int64).reshape(-1, 2)
    tof_idx = np.array([event.tof_idx for event in events], dtype=np.int64)
    return header, bins, tof_idx, np.array(blocks, dtype=np.int64).reshape(-1, 3)


def root_prefix(source):
    return DATA / f"ptsrc_{SOURCES[source]}_"


@pytest.fixture(scope="module")
def source(request):
    """index of the source (in SOURCES and POINTS): P4, unless the test is parametrized over the sources"""
    return getattr(request, "param", SOURCES.index("P4"))


@pytest.fixture(scope="module")
def cpp_file(source, tmp_path_factory):
    out = tmp_path_factory.mktemp(f"cpp_{SOURCES[source]}") / "ptsrc.petsird"
    subprocess.run([str(CONVERTER), "-r", str(root_prefix(source)), "-n", "1",
                    "-s", str(DATA / "geometry.json"), "-p", str(out)], check=True, capture_output=True)
    return out


@pytest.fixture(scope="module")
def python_file(source, tmp_path_factory):
    out = tmp_path_factory.mktemp(f"python_{SOURCES[source]}") / "ptsrc.petsird"
    # same selection as the C++ converter: true coincidences only
    subprocess.run([sys.executable, str(PYTHON_CONVERTER), "--root-prefix", str(root_prefix(source)),
                    "--number-of-root-files", "1", "-s", str(DATA / "geometry.json"),
                    "--time-mode", "absolute", "--selection", "true", "--petsird-file", str(out)],
                   check=True, capture_output=True)
    return out


@pytest.fixture(scope="module")
def converted(cpp_file):
    return read_petsird(cpp_file)


@pytest.fixture(scope="module")
def true_positions(source):
    """true emission position of every true coincidence in the ROOT file"""
    df = ROOT.RDataFrame("Coincidences", f"{root_prefix(source)}1.root")
    df = df.Filter("eventID1 == eventID2 && comptonPhantom1 == 0 && comptonPhantom2 == 0")
    cols = df.AsNumpy(["sourcePosX1", "sourcePosY1", "sourcePosZ1"])
    return np.stack([cols["sourcePosX1"], cols["sourcePosY1"], cols["sourcePosZ1"]], axis=1).astype(float)


@pytest.fixture(scope="module")
def lors(converted):
    """average detection points of detection_bins[0] (A) and detection_bins[1] (B) of every event"""
    header, bins, _, _ = converted
    cache = {}

    def adp(detection_bin):
        if detection_bin not in cache:
            expanded = expand_detection_bin(header.scanner, TYPE_OF_MODULE, detection_bin)
            cache[detection_bin] = np.array(get_ADP(header.scanner, TYPE_OF_MODULE, expanded).c, dtype=float)
        return cache[detection_bin]

    A = np.array([adp(int(b)) for b in bins[:, 0]]).reshape(-1, 3)
    B = np.array([adp(int(b)) for b in bins[:, 1]]).reshape(-1, 3)
    return A, B


def test_python_converter_gives_the_same_output(converted, python_file):
    header_cpp, bins_cpp, tof_cpp, blocks_cpp = converted
    header_py, bins_py, tof_py, blocks_py = read_petsird(python_file)

    np.testing.assert_array_equal(bins_py, bins_cpp, "detection bins")
    np.testing.assert_array_equal(tof_py, tof_cpp, "tof_idx")
    np.testing.assert_array_equal(blocks_py, blocks_cpp, "time blocks (start, stop, number of events)")

    scanner_cpp, scanner_py = header_cpp.scanner, header_py.scanner
    assert scanner_py.model_name == scanner_cpp.model_name
    module_cpp = scanner_cpp.scanner_geometry.replicated_modules[TYPE_OF_MODULE]
    module_py = scanner_py.scanner_geometry.replicated_modules[TYPE_OF_MODULE]

    def matrices(transforms):
        return np.array([np.asarray(transform.matrix) for transform in transforms])

    def corners(module):
        return np.array([corner.c for corner in module.object.detecting_elements.object.shape.corners])

    np.testing.assert_allclose(matrices(module_py.transforms), matrices(module_cpp.transforms),
                               atol=1e-3, err_msg="module transforms")
    np.testing.assert_allclose(matrices(module_py.object.detecting_elements.transforms),
                               matrices(module_cpp.object.detecting_elements.transforms),
                               atol=1e-3, err_msg="crystal transforms")
    np.testing.assert_allclose(corners(module_py), corners(module_cpp), atol=1e-3, err_msg="crystal box")
    np.testing.assert_allclose(scanner_py.tof_bin_edges[0][0].edges, scanner_cpp.tof_bin_edges[0][0].edges)
    assert np.isclose(scanner_py.tof_resolution[0][0], scanner_cpp.tof_resolution[0][0])
    np.testing.assert_allclose(scanner_py.event_energy_bin_edges[0].edges, scanner_cpp.event_energy_bin_edges[0].edges)
    assert np.isclose(scanner_py.energy_resolution_at_511[0], scanner_cpp.energy_resolution_at_511[0])
    assert scanner_py.prompt_event_policy == scanner_cpp.prompt_event_policy
    assert scanner_py.delayed_event_policy == scanner_cpp.delayed_event_policy
    assert [(m.id, m.name) for m in scanner_py.bulk_materials] == [(m.id, m.name) for m in scanner_cpp.bulk_materials]


def test_number_of_events(converted, true_positions):
    assert len(converted[1]) == len(true_positions)


@pytest.mark.parametrize("source", range(len(SOURCES)), ids=SOURCES, indirect=True)
@pytest.mark.parametrize("tof", [False, True], ids=["non-TOF", "TOF"])
def test_backprojection_maximum(source, converted, lors, tof):
    """PTSRC: backproject all the events with all efficiencies 1; the maximum of the image must be
    within 1 voxel (1 mm) of the true position on each axis.

    For TOF, PETSIRD stores (t1 - t2) * c / 2, with t1 the arrival time on detection_bins[0], and
    parallelproj measures the TOF position from the start to the end point of the LOR. With the
    symmetric PETSIRD TOF bins, bin j has the same centre in both, so tof_idx is used as is.
    """
    header, _, tof_idx, _ = converted
    A, B = lors
    scanner = header.scanner
    projector = ListmodePETProjector(A.astype(np.float32), B.astype(np.float32), IMAGE_SHAPE,
                                     (VOXEL_MM,) * 3, img_origin=IMAGE_ORIGIN.astype(np.float32))
    if tof:
        edges = np.asarray(scanner.tof_bin_edges[TYPE_OF_MODULE][TYPE_OF_MODULE].edges)
        projector.tof_parameters = TOFParameters(
            num_tofbins=len(edges) - 1, tofbin_width=float(edges[1] - edges[0]),
            sigma_tof=float(scanner.tof_resolution[TYPE_OF_MODULE][TYPE_OF_MODULE]) / 2.355)
        projector.event_tofbins = tof_idx.astype(np.int16)
        projector.tof = True
    image = np.asarray(projector.adjoint(np.ones(len(A), dtype=np.float32)))
    found = IMAGE_ORIGIN + np.array(np.unravel_index(np.argmax(image), image.shape)) * VOXEL_MM
    error = found - POINTS[source]
    assert np.all(np.abs(error) <= VOXEL_MM), \
        f"source at {POINTS[source]} mm: maximum of the backprojection at {found} mm (error {error.round(2)} mm)"


def test_tof(converted, lors, true_positions):
    """PETSIRD stores (t1 - t2) * c / 2, with t1 the arrival time on detection_bins[0].

    Checked event by event: the TOF backprojection hardly sees a wrong TOF sign, because for sources
    near the centre the mirrored TOF position is close to the source itself.
    """
    header, _, tof_idx, _ = converted
    A, B = lors
    S = true_positions
    edges = np.asarray(header.scanner.tof_bin_edges[TYPE_OF_MODULE][TYPE_OF_MODULE].edges)
    centres = 0.5 * (edges[:-1] + edges[1:])
    difference = centres[tof_idx] - 0.5 * (np.linalg.norm(S - A, axis=1) - np.linalg.norm(S - B, axis=1))
    # mean, not median: there is no TOF bin centred on 0 (the central ones are at -15 and +15 mm), and for a
    # source near the centre almost all events are in these 2 bins, so the median jumps to about +-5 mm
    assert abs(np.mean(difference)) < 2.0
    # expected spread: TOF resolution (30 mm FWHM) and bin width (30 mm), about 15 mm
    assert np.std(difference) < 20.0


def test_crystal_size(converted):
    header = converted[0]
    geometry = json.loads((DATA / "geometry.json").read_text())
    box = header.scanner.scanner_geometry.replicated_modules[TYPE_OF_MODULE].object.detecting_elements.object.shape
    corners = np.array([corner.c for corner in box.corners])
    sizes = corners.max(axis=0) - corners.min(axis=0)
    expected = [geometry["detector_x_dim"], geometry["detector_y_dim"], geometry["detector_z_dim"]]
    np.testing.assert_allclose(sizes, expected, atol=1e-3)


def test_crystals_do_not_overlap(converted):
    header = converted[0]
    elements = header.scanner.scanner_geometry.replicated_modules[TYPE_OF_MODULE].object.detecting_elements
    corners = np.array([corner.c for corner in elements.object.shape.corners])
    sizes = corners.max(axis=0) - corners.min(axis=0)
    # inside a module the crystals are only translated
    shifts = np.array([[transform.matrix[i][3] for i in range(3)] for transform in elements.transforms])
    overlap = np.all(np.abs(shifts[:, None, :] - shifts[None, :, :]) < sizes - 1e-3, axis=2)
    np.fill_diagonal(overlap, False)
    assert not overlap.any(), f"{overlap.sum() // 2} pairs of crystals overlap"
