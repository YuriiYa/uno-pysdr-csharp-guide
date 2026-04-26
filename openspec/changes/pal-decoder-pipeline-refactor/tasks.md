## 1. Test Project Scaffold

- [ ] 1.1 Create `PALDecoder.Tests/PALDecoder.Tests.csproj` (xUnit, net9.0) with references to `Microsoft.NET.Test.Sdk`, `xunit`, `xunit.runner.visualstudio`, and a `<ProjectReference>` to `uno_palyground`
- [ ] 1.2 Add `PALDecoder.Tests` to `uno_palyground.sln`
- [ ] 1.3 Verify `dotnet build PALDecoder.Tests/PALDecoder.Tests.csproj` succeeds and `dotnet test` finds zero tests (passing)

## 2. Create Pipeline Stage Classes

- [ ] 2.1 Create `uno_palyground/PySDRGuide/Pipeline/` folder and add `IqDemodulator.cs` — extract `ReadAndDemodFrame` logic; expose `bool TryRead(Stream source, int samplesPerFrame, double[] frameBuffer)`
- [ ] 2.2 Add `FrameSynchronizer.cs` — extract `FindFrameStart` and `RefineToFirstActiveLine` and `GoertzelPower`; expose `int FindFrameStart(double[] videoSignal, int sampleRate, int samplesPerLine)`
- [ ] 2.3 Add `HorizontalAligner.cs` — extract `EstimateHorizontalOffset`, `LightLowPass`, `MinMax`, `Median`; expose `int EstimateOffset(double[] videoSignal, int frameStart, int samplesPerLine, int sampleRate, int linesToUse = 24)`
- [ ] 2.4 Add `FieldSplitter.cs` — extract field-copy block; expose `bool TrySplit(double[] frameData, int samplesPerLine, double[] field1, double[] field2)` using Strategy A constants
- [ ] 2.5 Add `LumaChromaSeparator.cs` — extract `SeparateLumaChromaPooled`, `ApplyTwoFiltersStaticToDest`, `ApplyFilterStaticToDest`, `CreateLowPassFilter`, `CreateBandPassFilter`, and `FilterReverseCache`; expose `void Separate(double[] field, double[] lumaOut, double[] chromaOut)`
- [ ] 2.6 Move `ColorPll` into `Pipeline/ColorPll.cs` (already a separate class — move file, keep namespace `uno_palyground.PySDRGuide` or add `Pipeline` sub-namespace)
- [ ] 2.7 Add `ChromaDecoder.cs` — extract `DecodeChroma` and `FusedChromaDemodAndFilter`, `ApplySameFilterTwoSignalsStatic`; expose `void Decode(double[] chroma, int sampleRate, int samplesPerLine, int startLineOffset, double[] uOut, double[] vOut)`
- [ ] 2.8 Add `ActiveAreaCropper.cs` — extract `CropToActiveInto`; expose `void CropInto(double[] signal, int samplesPerLine, double[] dest, int copyWidth)`
- [ ] 2.9 Add `YuvToRgbConverter.cs` — extract `ConvertYUVToRGB_BT601_OptimizedInto` and `Clamp`; expose `void Convert(double[] y, double[] u, double[] v, int samplesPerLine, byte[,,] dest)`
- [ ] 2.10 Add `FieldInterleaver.cs` — extract `InterleaveFieldsInto`; expose `void Interleave(byte[,,] field1, byte[,,] field2, byte[,,] dest, FieldOrder order)`

## 3. Annotate Suspicious Logic

- [ ] 3.1 In `DecodePALSignal` (or `FieldSplitter`): add `// TODO: check the correctness: delta calculation mixes blanking intervals in a non-obvious way; verify against actual broadcast timing` next to the `delta` / `skipUntil + delta` lines
- [ ] 3.2 In `HorizontalAligner.EstimateOffset`: add `// TODO: check the correctness: 2.5*MAD threshold is empirically tuned, not derived from spec` next to `double thr = med - 2.5 * mad`
- [ ] 3.3 In `FrameSynchronizer.FindFrameStart`: add `// TODO: check the correctness: 3.0*MAD threshold is empirically tuned` next to `double syncThreshold = median - 3.0 * mad`
- [ ] 3.4 In `FieldSplitter`: add `// TODO: check the correctness: Strategy A assumes frameData starts at Field-1 active line 0 with no VBI; verify this holds for all stream sources` near the field-start-line constants
- [ ] 3.5 In `PALDecoder` constructor default: add `// TODO: check the correctness: BottomFieldFirst is the default but PAL broadcast may use TopFieldFirst depending on the source` next to the `FieldOrder.BottomFieldFirst` default

## 4. Refactor PALDecoder as Orchestrator

- [ ] 4.1 Replace inline `ReadAndDemodFrame` call in `DecodePALSignal` loop with `_iqDemodulator.TryRead(...)`
- [ ] 4.2 Replace inline `FindFrameStart` / `RefineToFirstActiveLine` calls with `_frameSynchronizer.FindFrameStart(...)`
- [ ] 4.3 Replace inline `EstimateHorizontalOffset` call with `_horizontalAligner.EstimateOffset(...)`
- [ ] 4.4 Replace field-copy block with `_fieldSplitter.TrySplit(...)`
- [ ] 4.5 Replace `Parallel.Invoke` inline luma/chroma blocks with calls to `_lumachromaSeparatorF1.Separate(...)` and `_lumachromaSeparatorF2.Separate(...)` (two instances for parallel field processing)
- [ ] 4.6 Replace inline chroma decode block with `_chromaDecoderF1.Decode(...)` and `_chromaDecoderF2.Decode(...)`
- [ ] 4.7 Replace `CropToActiveInto` calls with `_activeAreaCropper.CropInto(...)`
- [ ] 4.8 Replace `ConvertYUVToRGB_BT601_OptimizedInto` calls with `_yuvToRgbConverter.Convert(...)`
- [ ] 4.9 Replace `InterleaveFieldsInto` call with `_fieldInterleaver.Interleave(...)`
- [ ] 4.10 Verify `dotnet build uno_palyground.sln` succeeds with zero errors after all replacements

## 5. Unit Tests

- [ ] 5.1 Write `IqDemodulatorTests.cs`: test generic-stream path with a synthetic `MemoryStream` of known int8 IQ pairs; assert output magnitudes and DC removal
- [ ] 5.2 Write `IqDemodulatorTests.cs`: test short-read returns `false`
- [ ] 5.3 Write `FrameSynchronizerTests.cs`: build a synthetic video buffer with two broad pulses at a known position; assert `FindFrameStart` returns sample index within ±1 line
- [ ] 5.4 Write `FrameSynchronizerTests.cs`: test signal-too-short returns 0
- [ ] 5.5 Write `HorizontalAlignerTests.cs`: build a synthetic line with a 4.7 µs negative-polarity H-sync at a known column; assert returned offset is within ±2 samples
- [ ] 5.6 Write `FieldSplitterTests.cs`: fill a frame buffer with line-index markers; assert field1 and field2 boundaries are correct
- [ ] 5.7 Write `FieldSplitterTests.cs`: test insufficient-lines returns false
- [ ] 5.8 Write `LumaChromaSeparatorTests.cs`: DC-input → chroma output near zero test
- [ ] 5.9 Write `LumaChromaSeparatorTests.cs`: 4.43 MHz tone → chroma amplitude larger than luma test (at 10 MHz sample rate)
- [ ] 5.10 Write `ChromaDecoderTests.cs`: verify V negation on odd absolute lines using a controlled two-line chroma buffer
- [ ] 5.11 Write `ActiveAreaCropperTests.cs`: verify per-line copy copies exactly `copyWidth` samples from each line
- [ ] 5.12 Write `YuvToRgbConverterTests.cs`: white (Y=1,U=0,V=0) → (255,255,255); black (Y=0,U=0,V=0) → (0,0,0); clamped Y=1.5 → 255
- [ ] 5.13 Write `FieldInterleaverTests.cs`: BottomFieldFirst → row 0 = field2, row 1 = field1
- [ ] 5.14 Write `FieldInterleaverTests.cs`: TopFieldFirst → row 0 = field1, row 1 = field2
- [ ] 5.15 Run `dotnet test PALDecoder.Tests/PALDecoder.Tests.csproj` and confirm all tests pass

## 6. Verification

- [ ] 6.1 Run `dotnet build uno_palyground.sln` — zero errors, zero warnings introduced
- [ ] 6.2 Confirm `PALDecoder.DecodePALSignal(int, Stream)` and `PALDecoder.DecodePALSignal(int, FileStream)` remain the only public entry points with identical signatures
