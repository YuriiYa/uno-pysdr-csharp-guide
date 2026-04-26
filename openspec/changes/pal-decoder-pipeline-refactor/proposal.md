## Why

`PALDecoder` is a single monolithic class (~1500 lines) whose `DecodePALSignal` method weaves together reading, demodulation, sync detection, field separation, luma/chroma separation, chroma decode, cropping, YUV→RGB conversion, field interleaving, and display into one large loop body. This makes individual stages impossible to test in isolation, hard to reason about, and difficult to extend. Refactoring to an explicit **streaming pipeline** gives each stage a clear contract, makes the data-flow self-documenting, and unlocks unit testing of every processing step without hardware.

## What Changes

- Introduce a `IPalPipelineStage<TIn, TOut>` interface (or equivalent functional boundary) so each processing step is an independently named and testable unit.
- Extract the following named pipeline stages from `PALDecoder`:
  - **IqDemodulator** – reads raw IQ bytes → `double[]` AM-magnitude video samples, DC-removed.
  - **FrameSynchronizer** – detects VBI / V-sync in a video buffer → frame-start sample index.
  - **HorizontalAligner** – detects H-sync per line → per-line active-pixel column offset.
  - **FieldSplitter** – splits a video frame buffer into Field-1 and Field-2 contiguous sample arrays.
  - **LumaChromaSeparator** – FIR LPF (luma) + FIR BPF (chroma) → `(double[] Y, double[] C)`.
  - **ChromaDecoder** – ColorPLL demodulation + LPF → `(double[] U, double[] V)`.
  - **ActiveAreaCropper** – removes blanking from Y/U/V planes → active-region arrays.
  - **YuvToRgbConverter** – BT.601 matrix → `byte[,,]` RGB field.
  - **FieldInterleaver** – weaves two RGB fields into a full-frame `byte[,,]`.
- `PALDecoder` becomes a **pipeline orchestrator** that wires these stages together; its core loop shrinks to a sequence of clearly-named stage calls.
- Add an **xUnit test project** (`PALDecoder.Tests`) covering each stage with synthetic inputs:
  - Pure-tone IQ inputs for `IqDemodulator`.
  - Synthetic V-sync waveforms for `FrameSynchronizer`.
  - Synthetic H-sync lines for `HorizontalAligner`.
  - Known Y/U/V values for `YuvToRgbConverter` (spot-check BT.601 matrix).
  - Identity / passthrough cases for `ActiveAreaCropper` and `FieldInterleaver`.
- Logic is **preserved exactly** – no algorithm changes, only structural extraction. Suspicious numerical constants or invariants receive `// TODO: check the correctness:` comments inline.

## Capabilities

### New Capabilities

- `pal-pipeline-stages`: Named, independently testable PAL decoding pipeline stages extracted from `PALDecoder`.
- `pal-decoder-tests`: xUnit test project covering each pipeline stage with synthetic signal inputs.

### Modified Capabilities

*(none – no public API or spec-level behavior changes)*

## Impact

- **`uno_palyground/PySDRGuide/PALDecoder.cs`** – refactored; existing callers (`SecondPage.xaml.cs` or similar) call the same `PALDecoder.DecodePALSignal(sampleRate, stream)` entry point unchanged.
- **`NetHackRF/`** – unchanged.
- **New project**: `PALDecoder.Tests/PALDecoder.Tests.csproj` (xUnit, net9.0).
- **`uno_palyground.sln`** – gains the new test project reference.
- No new NuGet packages beyond xUnit and Microsoft.NET.Test.Sdk.
