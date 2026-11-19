"""Bounded parallel ROOT readers and a deterministic, single PETSIRD writer."""
from concurrent.futures import ThreadPoolExecutor, wait
from contextlib import nullcontext
import heapq
import math
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import time

import numpy as np
import petsird as p
import uproot

RECORD = np.dtype([('time', '<f8'), ('file', '<u4'), ('kind', 'u1'),
                   ('entry', '<u8'), ('a', '<u4'), ('b', '<u4'), ('tof', '<u4')])
FIELDS = [key + side for side in ('1', '2') for key in
          ('rsectorID', 'moduleID', 'submoduleID', 'crystalID', 'layerID', 'gantryID', 'energy', 'time')]


def open_root(path):
    # The outer reader pool owns parallelism, including in serial mode.
    from uproot.source.futures import TrivialExecutor
    return uproot.open(path, handler=uproot.source.file.MemmapSource,
                       num_workers=1, decompression_executor=TrivialExecutor(),
                       interpretation_executor=TrivialExecutor())


def backend(requested):
    if requested in ('auto', 'gpu'):
        try:
            import cupy as cp
            if cp.cuda.runtime.getDeviceCount() > 0:
                return 'gpu', cp
        except Exception as error:
            if requested == 'gpu':
                raise RuntimeError('GPU mode needs a working CuPy CUDA installation') from error
        if requested == 'gpu':
            raise RuntimeError('GPU mode requested but no CUDA device is available')
    return ('cpu' if requested == 'auto' else requested), np


def bounds(path, include_delayeds, limit):
    lo, hi, has_delays = math.inf, -math.inf, False
    with open_root(path) as root:
        if 'Coincidences' not in root:
            raise ValueError(f'Missing Coincidences tree in {path}')
        has_delays = 'delay' in root
        for name in ('Coincidences', 'delay') if include_delayeds and has_delays else ('Coincidences',):
            for rows in root[name].iterate(['time1'], entry_stop=limit, step_size=500000, library='np'):
                values = rows['time1'] * 1000.
                if not np.isfinite(values).all() or np.any(values < 0):
                    raise ValueError(f'Invalid timestamps in {path}, tree {name}')
                if len(values):
                    lo, hi = min(lo, float(values.min())), max(hi, float(values.max()))
    return lo, hi, has_delays


def offsets(spans, mode, duration_ms):
    result, offset, previous = [], 0., -math.inf
    for index, (lo, hi, _) in enumerate(spans):
        if mode == 'per-file':
            offset = index * duration_ms
            if hi >= duration_ms:
                raise ValueError('Per-file clock exceeds --file-duration-ms')
        elif mode == 'auto' and lo + offset < previous:
            # Anchor reset segments on the specified acquisition-file period.
            offset = math.ceil(previous / duration_ms) * duration_ms
        result.append(offset)
        if math.isfinite(hi):
            previous = max(previous, hi + offset)
    return result


def transform(rows, geometry, categories, offset, start, stop, delayed_ps, xp):
    values = {name: xp.asarray(rows[name]) for name in rows.dtype.names} if isinstance(rows, np.ndarray) and rows.dtype.names else {
        name: xp.asarray(value) for name, value in rows.items()}
    if any(not bool(xp.all(xp.isfinite(values['time'+side]))) for side in ('1','2')):
        raise ValueError('Nonfinite ROOT coincidence timestamps')
    ids = []
    for side in ('1', '2'):
        fields = [values[key + side] for key in ('rsectorID', 'moduleID', 'submoduleID', 'crystalID', 'layerID')]
        sizes = [geometry.modules, *(x*z for x,z in zip(geometry.xy, geometry.z)), 1]
        if bool(xp.any(values['gantryID' + side] != 0)) or any(
                bool(xp.any((v < 0) | (v >= n))) for v,n in zip(fields, sizes)):
            raise ValueError('ROOT hierarchy index outside scanner JSON')
        fields[0] = fields[0] % geometry.rsec_xy * geometry.rsec_z + fields[0] // geometry.rsec_xy
        detector = xp.zeros(len(fields[0]), dtype=xp.int64)
        for value, size in zip(fields, sizes):
            detector = detector * size + value
        ids.append(detector)
    a, b = ids
    clock = values['time1'].astype(xp.float64) * 1000. + offset
    dt = (values['time1'].astype(xp.float64) - values['time2'].astype(xp.float64))*1e12 + delayed_ps
    displacement = xp.where(a < b, -dt, dt) * .299792458 / 2
    tof = xp.searchsorted(xp.asarray(geometry.tof), displacement, side='right') - 1
    keep = (clock >= start) & (clock < stop) & (a != b) & (tof >= 0) & (tof < len(geometry.tof)-1)
    for side in ('1', '2'):
        energy = values['energy'+side].astype(xp.float64) * 1000.
        keep &= (energy >= geometry.energy[0]) & (energy < geometry.energy[-1])
    if categories is not None:
        category = xp.where(values['eventID1'] != values['eventID2'], 3,
                            xp.minimum(values['comptonPhantom1'] + values['comptonPhantom2'], 2))
        keep &= xp.isin(category, xp.asarray(categories))
    take = xp.flatnonzero(keep)
    order = xp.lexsort(xp.stack((take, clock[take])))
    take = take[order]
    host = (lambda v: xp.asnumpy(v)) if xp is not np else np.asarray
    return tuple(host(v) for v in (take, clock[take], xp.maximum(a,b)[take], xp.minimum(a,b)[take], tof[take]))


def spool_file(job):
    index, path, geometry, args, offset, has_delays, directory, mode = job
    xp = np
    if mode == 'gpu':
        import cupy as xp
    streams = [('Coincidences', 0)] + ([('delay', 1)] if has_delays else [])
    runs, totals = [], [0, 0]
    stream = xp.cuda.Stream(non_blocking=True) if mode == 'gpu' else nullcontext()
    with stream, open_root(path) as root:
        for name, kind in streams:
            fields = FIELDS.copy()
            categories = None
            if kind == 0 and len(args.selection) != 4:
                fields += [key+side for side in ('1','2') for key in ('eventID','comptonPhantom')]
                categories = [('true','single-scatter','multiple-scatter','random').index(c) for c in args.selection]
            entry = 0
            for rows in root[name].iterate(fields, entry_stop=args.max_events_per_tree,
                                           step_size=args.batch_size, library='np'):
                take, clock, a, b, tof = transform(rows, geometry, categories, offset,
                    args.start_time_s*1000., args.end_time_s*1000.,
                    (args.delayed_offset_ps or 0.) if kind else 0., xp)
                data = np.empty(len(take), dtype=RECORD)
                data['time'], data['file'], data['kind'] = clock, index, kind
                data['entry'], data['a'], data['b'], data['tof'] = take+entry, a, b, tof
                if len(data):
                    target = Path(directory) / f'{index:06d}_{kind}_{entry:012d}.npy'
                    np.save(target, data)
                    runs.append(str(target))
                    totals[kind] += len(data)
                entry += len(rows['time1'])
    return runs, {'file': str(path), 'offset_ms': offset, 'prompts': totals[0], 'delayeds': totals[1]}


def ordered_blocks(runs, block_ms):
    """Merge bounded runs at block granularity, sorting ties by original file and entry."""
    arrays = [np.load(path, mmap_mode='r') for path in runs]
    heap = [(float(a['time'][0]), i, 0) for i,a in enumerate(arrays) if len(a)]
    heapq.heapify(heap)
    try:
        while heap:
            tag = math.floor(heap[0][0]/block_ms)*block_ms
            parts = []
            while heap and heap[0][0] < tag+block_ms:
                _, index, begin = heapq.heappop(heap)
                data = arrays[index]
                end = int(np.searchsorted(data['time'], tag+block_ms, side='left'))
                parts.append(data[begin:end])
                if end < len(data):
                    heapq.heappush(heap, (float(data['time'][end]), index, end))
            data = np.concatenate(parts)
            order = np.lexsort((data['entry'], data['kind'], data['file'], data['time']))
            yield tag, data[order]
    finally:
        for a in arrays:
            a._mmap.close()


def run(args, geometry, header, report):
    mode, _ = backend(args.execution)
    workers = 1 if mode == 'none' else (args.workers or min(8, os.cpu_count() or 1))
    paths = [Path(f'{args.root_prefix}{i}.root') for i in range(1,args.number_of_root_files+1)]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    output = args.petsird_file
    partial, manifest = output.with_suffix('.partial'), output.with_suffix('.petsird.json')
    if any(path.exists() for path in (output,partial,manifest)):
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    with (ThreadPoolExecutor(max_workers=workers) if mode != 'none' else nullcontext(None)) as pool:
        probe = lambda path: bounds(path, args.delayeds != 'off', args.max_events_per_tree)
        spans = list(pool.map(probe, paths)) if pool else [probe(path) for path in paths]
        if args.delayeds == 'required' and any(not s[2] for s in spans):
            raise ValueError('A required delay tree is missing')
        shifts = offsets(spans, args.time_mode, args.file_duration_ms)
        has_delays = [span[2] and args.delayeds != 'off' for span in spans]
        header.scanner.delayed_event_policy = p.CoincidencePolicy.OTHER if any(has_delays) else p.CoincidencePolicy.NONE
        report.update(execution_requested=args.execution, execution_used=mode, workers=workers,
                      selected_interval_s=[args.start_time_s, args.end_time_s if math.isfinite(args.end_time_s) else None],
                      clock_offsets_ms=shifts)
        print(f'Execution: {mode}, {workers} ROOT readers; timestamp ordering before PETSIRD serialization', flush=True)
        with TemporaryDirectory(prefix='gateroot-', dir=args.temp_dir) as temporary:
            jobs = [(i,path,geometry,args,shifts[i],has_delays[i],temporary,mode) for i,path in enumerate(paths)]
            futures = [pool.submit(spool_file, job) for job in jobs] if pool else []
            runs = []
            try:
                results = (future.result() for future in futures) if pool else map(spool_file,jobs)
                for paths_written, item in results:
                    runs += paths_written
                    report['files'].append(item)
                    print(item, flush=True)
            finally:
                for future in futures:
                    future.cancel()
                wait(futures)
            report['read_transform_seconds'] = time.perf_counter()-started
            try:
                with p.BinaryPETSIRDWriter(str(partial)) as writer:
                    writer.write_header(header)
                    for tag, data in ordered_blocks(runs, geometry.block_ms):
                        if tag < 0 or tag+geometry.block_ms > np.iinfo(np.uint32).max:
                            raise ValueError('Acquisition timestamps exceed PETSIRD uint32 milliseconds')
                        events = []
                        for kind in (0,1):
                            rows = data[data['kind'] == kind]
                            events.append([p.CoincidenceEvent(detection_bins=[int(row['a']),int(row['b'])],
                                                              tof_idx=int(row['tof'])) for row in rows])
                        block = p.EventTimeBlock(time_interval=p.TimeInterval(start=tag,stop=tag+geometry.block_ms),
                                                 prompt_events=[[events[0]]],
                                                 delayed_events=[[events[1]]] if any(has_delays) else [])
                        writer.write_time_blocks([p.TimeBlock.EventTimeBlock(block)])
                    writer.write_time_blocks([])
                if args.verify:
                    from root_to_petsird import verify
                    report.update(verify(partial, geometry, report))
                partial.rename(output)
            except BaseException:
                partial.unlink(missing_ok=True)
                raise
    report['total_seconds'] = time.perf_counter()-started
    import json
    manifest.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(f'Saved {output} and {manifest}', flush=True)
