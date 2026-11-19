GATEROOT 1.0.0 supports PETSIRD 0.13.3 through the Python and C++ converters.

- `--execution auto|gpu|cpu|none` selects processing. Auto prefers CUDA, then CPU readers. ROOT decoding uses CPUs. CUDA accelerates selection, TOF processing and sorting; Python also maps detector IDs on the GPU.
- Independent readers use bounded batches and temporary sorted runs. One writer merges timestamps deterministically. Equal timestamps retain file/entry order. Prompt and delayed streams retain their own chronological order in PETSIRD.
- `--time-mode auto` repairs inter-file clock resets. `per-file` places file i at `(i-1)*file_duration_ms`. `absolute` preserves timestamps and merges overlapping files.
- `--start-time-s` and `--end-time-s` select a half-open interval after clock repair. The default covers all supplied files.
- Delay trees are detected automatically. `--delayeds required` requires their presence. `--delayed-offset-ps` recenters the delayed window before TOF encoding. Prompts default to true + single scatter + multiple scatter + randoms.
- Explicit unit detection-bin efficiencies accompany expanded normalization in both implementations.
- Tests cover serial/CPU parity, reset-clock windows, stable ordering, delayeds, and optional GPU/cross-language comparisons. CI builds CUDA and CPU variants. GPU runtime tests require hardware and are explicitly skipped on CPU-only runners.

Examples cover IEC 6 minutes, voxBrain 6 minutes, and the first 180 seconds selected from all 36 voxBrain files. The optional EMITRON workflow estimates scatter from a preliminary attenuation-corrected two-iteration image and compiles an additive-only correction sidecar for the final reconstruction.
