## Context

`PALDecoder.cs` (~1500 lines) processes PAL TV signals captured by a HackRF SDR. The decoding pipeline is logically sequential:

> IQ bytes → AM demod → frame sync → H-align → field split → luma/chroma FIR → chroma PLL → crop → YUV→RGB → interleave → display

Currently all these stages live inside a single `DecodePALSignal` loop body, sharing mutable field-level buffer references (`_field1Buffer`, `_luminanceBufferF1`, …) and profiling accumulators. There are no unit tests. The `ColorPll` class is already extracted; everything else is private methods on `PALDecoder`.

The main entry point `PALDecoder.DecodePALSignal(int sampleRate, Stream stream)` must remain unchanged so existing callers are unaffected.

## Goals / Non-Goals

**Goals:**
- Extract each logical processing stage into a dedicated, named class/struct with a clear input/output contract.
- Keep the same buffer-reuse and SIMD optimisation strategy so performance is not regressed.
- Provide an xUnit test project that exercises each stage with synthetic inputs (no hardware required).
- Mark any suspicious algorithm constants or invariants with `// TODO: check the correctness: <reason>`.

**Non-Goals:**
- Changing any numerical algorithm (filter taps, PLL gains, sync thresholds, BT.601 coefficients).
- Introducing async/dataflow (`System.Threading.Channels`) – the pipeline stays synchronous/parallel in the same thread model as today.
- Adding a formal `IPipelineStage<TIn, TOut>` generic interface – a simple set of well-named classes is sufficient and avoids over-engineering.
- Porting the Python scripts.

## Decisions

### D1: Stage granularity – classes, not lambdas

**Decision:** Each stage becomes a `public sealed class` (or `internal sealed class` for stages that need no external visibility) inside a new `PALDecoder.Pipeline` namespace within the same file (or a `Pipeline/` subfolder).

**Rationale:** Classes allow the test project to `new StageX(params).Process(input)` without reflection or mocking frameworks. Lambdas or local functions cannot be tested in isolation.

**Alternative considered:** A single `IPipelineStage<TIn, TOut>` generic interface. Rejected because the type variety of stage I/O (byte[], double[], Complex[], byte[,,]) makes a uniform generic interface either too abstract or forced.

### D2: Buffer ownership – stages accept pre-allocated destination buffers

**Decision:** Stages that benefit from buffer reuse accept optional `ref double[]?` destination parameters (same `EnsureBuffer` pattern already present). Stages that always produce new output (e.g., `YuvToRgbConverter`) allocate internally and the orchestrator holds the buffer refs.

**Rationale:** Mirrors the existing performance strategy; avoids GC pressure on hot path. Tests can pass `null` and receive freshly allocated results.

### D3: `PALDecoder` becomes the orchestrator

**Decision:** `PALDecoder` retains the public API (`DecodePALSignal`) and holds all reusable buffer fields and profiling accumulators. Internally it replaces the large inline code blocks with stage instances instantiated in the constructor (filter-cache and PLL are stateful).

**Rationale:** Zero-change surface for existing callers; all complexity moves into the individually-testable stage classes.

### D4: Test project location and framework

**Decision:** New project `PALDecoder.Tests/PALDecoder.Tests.csproj` using **xUnit** (already commonly used in the .NET ecosystem, no new framework). Added to `uno_palyground.sln`. References `uno_palyground` project directly (or a shared `PALDecoder.Pipeline` project if extraction warrants it).

**Rationale:** xUnit requires no external test runner configuration beyond the SDK; `dotnet test` works out-of-the-box.

### D5: Suspicious constants annotated, not changed

**Decision:** During extraction, any constant, threshold, or formula that cannot be trivially verified from the PAL spec receives a `// TODO: check the correctness: <reason>` comment. Examples:
- `skipUntil + delta` offset arithmetic in `DecodePALSignal` (empirical delta calculation).
- `2.5 * mad` threshold in `EstimateHorizontalOffset` (tuned heuristic).
- BFF/TFF assumption default (`BottomFieldFirst`) vs. actual broadcast norm.
- `field2StartLine = linesPerFieldAll` under Strategy A (assumes no VBI in frameData).

## Risks / Trade-offs

- **[Risk] Parallel buffer access** – `Parallel.Invoke` on F1/F2 uses separate `_*F1` / `_*F2` buffer pairs. If a stage class holds instance-level scratch buffers and is shared between parallel branches, data races would occur. → **Mitigation:** Each field branch gets its own stage instance (or stage methods are static / accept caller-supplied buffers).

- **[Risk] Filter cache thread safety** – `ConcurrentDictionary`-based `_filterCache` is already thread-safe. Extraction must not break this. → **Mitigation:** Move `_filterCache` to a `static` shared cache or pass the same cache instance to stages that need it.

- **[Trade-off] File size** – Keeping all stages in one file keeps diffs reviewable but makes the file long. Splitting into `Pipeline/*.cs` improves discoverability but adds files. Decision: split into separate files under `PySDRGuide/Pipeline/`.

## Migration Plan

1. Create `PALDecoder.Tests` project and ensure it builds (empty) before any stage extraction.
2. Extract stages one at a time, running `dotnet build` and all existing tests after each.
3. After all stages are extracted, run the app manually on a recorded IQ file to verify no visual regression.
4. Rollback: revert the extracted stage class and restore the inline code; the public API is unchanged throughout.

## Open Questions

- Should `ColorPll` (already a separate class) move into the `Pipeline/` subfolder? Preference: yes, for consistency.
- Should the test project reference `uno_palyground` directly (easy) or require extracting `PALDecoder` into its own class library (cleaner but larger scope)? Decision deferred to tasks — start with direct project reference.
