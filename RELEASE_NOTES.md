# GATEROOT 1.0.0

Official PETSIRD 0.13.3 support for both standalone implementations.

- Select automatic GPU/CPU, explicit GPU, CPU reader multithreading or serial execution.
- Independent CPU ROOT readers, optional CUDA/CuPy transforms and deterministic single-writer output.
- Bounded sorted spools and external merging preserve corrected chronological order across files and batches.
- Repair inter-file acquisition clock resets by default; preserve continuous clocks. Explicit absolute/per-file modes remain available.
- Select a half-open acquisition time window after clock correction; omit it to retain the full acquisition.
- Include all prompts and available delayed coincidences by default; explicit delayed offset, required/off options and truth selections remain available.
- Preserve latest multi-axial-rsector geometry mapping and embedded geometric normalization; explicit unit detector efficiencies in both converters.
- Conversion audit reports include execution mode, workers, clock shifts, per-file event totals and timing.
- CPU/CUDA build CI, ROOT-derived serial/parallel parity, window/delayed tests, and existing four-point-source localization/TOF tests.
- Windows Python and WSL C++ examples for IEC 6 min and Hoffman 3/6 min, followed by EMITRON preliminary OSEM, SSS, additive-only sidecar and final OSEM with sensitivity reuse.

GPU acceleration does not move ROOT decompression or official PETSIRD serialization onto the GPU.
Explicit GPU mode needs working GPU dependencies; automatic mode falls back to CPU.
Python CUDA runtime parity was tested on RTX 5090; hosted CI compiles C++ CUDA but has no device.
Run the documented C++ GPU parity test on your workstation before large device conversions.
The six large SSS/reconstruction examples are provided for local validation and are not CI numerical validation of scatter correction.

See [the execution guide](docs/PARALLEL_CONVERSION.md) and [workstation instructions](docs/WORKSTATION.md).
