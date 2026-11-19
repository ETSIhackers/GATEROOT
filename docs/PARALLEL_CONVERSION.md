# Ordered parallel conversion â€” GATEROOT 1.0.0

Both Python and C++ use the official PETSIRD 0.13.3 SDK. Existing geometry,
selection, ROOT prefix, file count and normalization arguments remain available.
The defaults for execution, clock handling and delayed inclusion below are new.

| Argument | Default | Meaning |
|---|---|---|
| `--execution auto\|gpu\|cpu\|none` | `auto` | Use CUDA if usable, otherwise parallel CPU readers; `none` uses one reader and CPU transforms |
| `--workers N` | `0` | 0 chooses up to 8 CPU cores; explicit positive N controls independent ROOT readers |
| `--batch-size N` | 500000 | Maximum input events per reader batch |
| `--temp-dir PATH` | system temporary directory | Sorted event spool storage, cleaned after success or failure |
| `--time-mode auto\|absolute\|per-file` | `auto` | Repair detected inter-file clock resets, preserve raw clocks, or explicitly place each file at its numbered segment |
| `--file-duration-ms N` | 10000 | Segment spacing used by reset repair / per-file clocks; must match the simulation segments |
| `--start-time-s S` | 0 | Inclusive acquisition-window start, after clock correction |
| `--end-time-s E` | unbounded | Exclusive acquisition-window end, after clock correction |
| `--selection ...` | `all` | Prompt classes: all, true, scatter, single-scatter, multiple-scatter, random; comma-separated unions supported |
| `--delayeds auto\|required\|off` | `auto` | Store `delay` ROOT trees if present, require them in every file, or explicitly omit them |
| `--delayed-offset-ps D` | 0 | Add D to delayed `time1-time2` before detector-order canonicalization; providing D requires a delay tree |
| `--version` | â€” | Print converter / SDK version information |

GPU mode is a hybrid pipeline: ROOT decompression and file I/O run in independent
CPU readers. Python uses CuPy for event classification, detector mapping, time/energy/TOF
filtering and batch ordering; C++ uses a CUDA kernel for event filtering, TOF
assignment and endpoint ordering, and Thrust for sorting. C++ detector hierarchy
mapping and truth classification remain on CPU. Disk I/O and official SDK
serialization do not execute on the GPU. C++ bounds simultaneous device work to
one batch; Python uses one CUDA stream per reader. Explicit `gpu` fails clearly
if unavailable; only `auto` falls back. A CPU-only C++ build cannot select GPU.
Python GPU mode needs the separately installed matching CuPy package.

This does not guarantee GPU is faster: disk bandwidth, decompression, small
batches, normalization loading, merging and PETSIRD writing can dominate. Start
with 2â€“4 readers on an external hard disk. More readers can increase disk seeks.
Serial mode disables the converter's reader pool and CUDA; it does not promise
that every linked runtime uses exactly one operating-system thread.

## Chronology and acquisition windows

Input files are numbered `PREFIX1.root` through `PREFIXN.root`, in numeric order.
A timestamp-only prepass examines `time1` of prompts and enabled delayeds before
filtering. In auto mode, a file whose minimum timestamp falls below the previous
corrected maximum is shifted to the next segment boundary (ceiling of the
previous corrected maximum divided by `file-duration-ms`). Subsequent continuous
files retain that segment's offset. Empty files do not advance the clock.
This assumes successive acquisition segments; genuinely overlapping independent
streams should use `absolute` to retain their original clock values.

For example two reset 10-second files have corrected spans approximately
`[0,10)` and `[10,20)`. Selecting `[12,15)` retains only the requested three seconds
of the second file. Continuous-clock files receive no unnecessary shift. The
prepass still checks every specified input file; files wholly outside the
selected interval then skip full event decoding. Omit both time arguments to
convert the complete corrected acquisition.

Selected events are ordered by corrected `time1`, then input file number, stream
(prompts before delayeds), and ROOT entry number for ties. Sorted batch spools
are merged by one writer. Worker completion order never determines output order.
PETSIRD stores a time interval per event block, not each event's original floating
timestamp: selection uses the full original corrected timestamp **before** block
quantization. Boundary blocks may straddle a fractional requested boundary, but
contain only events that passed the exact half-open selection. Within each
prompt/delayed vector the original corrected chronology is retained. Interleaving
between the two vector types is represented by their shared block, not a separate
per-event timestamp field.

The clock shift changes acquisition times, not the within-coincidence arrival
difference. Delayed sorter recentering is a separate operation. The ETSIPET data
in the workstation example use a 500 ns delayed offset (`500000` ps). A different
sorter may require a different value. Delayeds are not filtered by the requested
prompt truth class; both streams receive the same energy/TOF/time-window cuts.
There is no synthetic delayed stream if the ROOT tree is absent.

Output detector IDs retain the converter's existing descending endpoint order;
swapping endpoints reverses the TOF sign. Prompt truth labels are selection
inputs, not separate PETSIRD true/scatter/random fields. Prompts remain prompts;
the delayed stream is stored in PETSIRD `delayed_events`.

## Resources, reports and validation

Conversion uses bounded reader batches plus disk-backed sorted runs. Memory grows
with workers Ã— batch size, normalization matrices, and events in one output time
block; a huge time-block duration can still create a huge block. Temporary disk
holds roughly 40 bytes per selected coincidence and, during merge passes, may
need an additional copy. Allow generous free space for large acquisitions, plus
the output. Use distinct output and temporary directories for concurrent runs.
The single C++ writer merges externally with at most 32 inputs per merge pass.
Neither converter overwrites an existing output. An interrupted `.partial` file
must be inspected/removed or avoided by choosing a fresh output filename.

Each output has an adjacent `.petsird.json` conversion report with execution mode,
workers, per-file clock offsets, selected interval, prompt/delayed counts and
elapsed times. C++ and Python are tested for decoded header/event/block equality;
identical file bytes or equal file sizes are not the correctness criterion.
Both explicitly store unit detector efficiencies when normalization is supplied.

Run `python -m pytest -v test` after building C++. Existing four-point-source tests
verify geometry, TOF sign and image localization. New tests derive shuffled,
reset-clock prompt/delayed inputs from the included GATE ROOT data, compare
serial/CPU implementations and post-shift windows, and check GPU parity when
hardware is present. `GATEROOT_TEST_CPP_GPU=1` enables C++ CUDA runtime parity.
CI builds CPU and CUDA variants with SDK 0.13.3; hosted CI has no GPU and therefore
does not run device kernels. The local Python CUDA test was run on an RTX 5090.
