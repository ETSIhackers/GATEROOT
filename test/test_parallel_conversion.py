"""Parity and ordering checks derived from the repository's GATE point-source data."""
import importlib.util
import json
import io
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import petsird
import pytest
import uproot

REPO = Path(__file__).parents[1]
sys.path.insert(0, str(REPO / 'python'))
from parallel_conversion import offsets, backend


@pytest.fixture
def input_pair(tmp_path):
    source = REPO / 'test/data/ptsrc_P1_1.root'
    geometry = REPO / 'test/data/geometry.json'
    if not source.exists():
        source = Path(os.environ.get('GATEROOT_REFERENCE_ROOT', str(source)))
        geometry = Path(os.environ.get('GATEROOT_REFERENCE_GEOMETRY', str(geometry)))
    if not source.exists():
        pytest.skip('Repository GATE fixture unavailable')
    with uproot.open(io.BytesIO(source.read_bytes())) as root:
        fields = [key+side for side in ('1','2') for key in
                  ('rsectorID','moduleID','submoduleID','crystalID','layerID','gantryID',
                   'energy','time','eventID','comptonPhantom','globalPosX','globalPosY','globalPosZ')]
        arrays = root['Coincidences'].arrays(fields, entry_stop=200, library='np')
    data = {name:arrays[name].copy() for name in fields}
    count = len(data['time1'])
    original_dt = data['time1']-data['time2']
    # Deliberately scramble timestamps within batches, including exact ties.
    order = np.arange(count)[::-1]
    clock = (np.arange(count)//2)*.005
    data['time1'] = clock[order]
    data['time2'] = data['time1']-original_dt
    delayed = {name:value.copy() for name,value in data.items()}
    delayed['time2'] += 500e-9
    for index in (1,2):
        with uproot.recreate(tmp_path / f'events_{index}.root') as root:
            for name,rows in [('Coincidences',data),('delay',delayed)]:
                root.mktree(name, {key:value.dtype for key,value in rows.items()})
                root[name].extend(rows)
    return tmp_path/'events_', geometry


def decoded(path):
    events, blocks = [], []
    with petsird.BinaryPETSIRDReader(str(path)) as reader:
        reader.read_header()
        for union in reader.read_time_blocks():
            block = union.value
            blocks.append((block.time_interval.start,block.time_interval.stop))
            for kind,rows in enumerate((block.prompt_events,block.delayed_events)):
                for row in rows:
                    for group in row:
                        events += [(block.time_interval.start,kind,*event.detection_bins,event.tof_idx)
                                   for event in group]
    return events,blocks


def convert(tmp_path,pair,mode,implementation='python',window=()):
    prefix, geometry = pair
    output = tmp_path / f'{implementation}_{mode}_{len(window)}.petsird'
    command = ([sys.executable,str(REPO/'python/root_to_petsird.py')] if implementation=='python'
               else [os.environ.get('ROOT_TO_PETSIRD',str(REPO/'cpp/build/root_to_petsird'))])
    command += ['--root-prefix',str(prefix),'--number-of-root-files','2','-s',str(geometry),
                '--petsird-file',str(output),'--execution',mode,'--workers','2','--batch-size','37',
                '--file-duration-ms','2000','--delayeds','required','--delayed-offset-ps','500000',*window]
    completed = subprocess.run(command,capture_output=True,text=True,timeout=180)
    if completed.returncode:
        raise RuntimeError(completed.stdout+completed.stderr)
    return decoded(output)


def test_clock_offsets():
    assert offsets([(0,9999,True),(0,9999,True),(0,9999,True)],'auto',10000)==[0,10000,20000]
    assert offsets([(0,9999,True),(10000,19999,True)],'auto',10000)==[0,0]
    assert offsets([(0,9999,True),(0,9999,True)],'absolute',10000)==[0,0]


def test_python_cpu_serial_parity_and_window(tmp_path,input_pair):
    serial = convert(tmp_path,input_pair,'none')
    assert convert(tmp_path,input_pair,'cpu')==serial
    events,blocks=serial
    assert events and all(b>a for a,b in blocks)
    assert all(blocks[i][0]>=blocks[i-1][1] for i in range(1,len(blocks)))
    assert any(row[1]==1 for row in events)
    window=['--start-time-s','2','--end-time-s','2.25']
    chosen, _=convert(tmp_path,input_pair,'cpu',window=window)
    assert chosen == [row for row in events if 2000<=row[0]<2250]


def test_python_gpu_parity(tmp_path,input_pair):
    if backend('auto')[0]!='gpu':
        pytest.skip('CUDA GPU and CuPy are required')
    assert convert(tmp_path,input_pair,'gpu')==convert(tmp_path,input_pair,'none')
    window=['--start-time-s','2','--end-time-s','2.25']
    assert convert(tmp_path,input_pair,'gpu',window=window)==convert(tmp_path,input_pair,'none',window=window)


def test_cpp_cpu_python_parity(tmp_path,input_pair):
    executable=Path(os.environ.get('ROOT_TO_PETSIRD',REPO/'cpp/build/root_to_petsird'))
    if not executable.exists():
        pytest.skip('Build the C++ converter first')
    reference=convert(tmp_path,input_pair,'none')
    assert convert(tmp_path,input_pair,'none','cpp')==reference
    assert convert(tmp_path,input_pair,'cpu','cpp')==reference
    assert convert(tmp_path,input_pair,'cpu','cpp',window=['--start-time-s','2','--end-time-s','2.25']) == convert(
        tmp_path,input_pair,'cpu',window=['--start-time-s','2','--end-time-s','2.25'])


def test_cpp_gpu_python_parity(tmp_path,input_pair):
    if os.environ.get('GATEROOT_TEST_CPP_GPU')!='1':
        pytest.skip('Set GATEROOT_TEST_CPP_GPU=1 on a CUDA-enabled C++ build')
    assert convert(tmp_path,input_pair,'gpu','cpp')==convert(tmp_path,input_pair,'none')
    window=['--start-time-s','2','--end-time-s','2.25']
    assert convert(tmp_path,input_pair,'gpu','cpp',window=window)==convert(tmp_path,input_pair,'none',window=window)
