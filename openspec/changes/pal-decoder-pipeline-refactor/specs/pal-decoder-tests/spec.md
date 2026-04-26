## ADDED Requirements

### Requirement: Test project exists and builds
A `PALDecoder.Tests` xUnit test project SHALL exist, target `net9.0`, and build cleanly with `dotnet build` without modifying any production project.

#### Scenario: Empty test suite compiles
- **WHEN** `dotnet build PALDecoder.Tests/PALDecoder.Tests.csproj` is executed
- **THEN** the build exits with code 0

#### Scenario: dotnet test runs without errors
- **WHEN** `dotnet test PALDecoder.Tests/PALDecoder.Tests.csproj` is executed
- **THEN** all tests pass (zero failures, zero errors)

---

### Requirement: Each pipeline stage has at least one passing unit test
Every stage class extracted from `PALDecoder` (`IqDemodulator`, `FrameSynchronizer`, `HorizontalAligner`, `FieldSplitter`, `LumaChromaSeparator`, `ChromaDecoder`, `ActiveAreaCropper`, `YuvToRgbConverter`, `FieldInterleaver`) SHALL have at least one xUnit `[Fact]` or `[Theory]` test exercising its primary functionality with a synthetic (no-hardware) input.

#### Scenario: IqDemodulator test
- **WHEN** the IqDemodulator test is executed with a synthetic in-memory stream of interleaved int8 IQ bytes
- **THEN** the test passes without accessing any hardware or file system

#### Scenario: YuvToRgbConverter BT.601 white test
- **WHEN** Y=1, U=0, V=0 values fill the input arrays
- **THEN** all output RGB bytes equal 255 and the test passes

#### Scenario: FieldInterleaver order test
- **WHEN** two synthetic field arrays with distinct pixel values are interleaved
- **THEN** the test asserts the correct BottomFieldFirst or TopFieldFirst row ordering and passes

---

### Requirement: Tests are independent of hardware and UI
All tests in `PALDecoder.Tests` SHALL run without a connected HackRF device, without `hackrf.dll` present, and without a `DispatcherQueue` / WinUI context.

#### Scenario: Tests run in CI headless environment
- **WHEN** `dotnet test` is executed on a machine without HackRF hardware
- **THEN** all tests pass (no `PlatformNotSupportedException` or missing-DLL errors)

---

### Requirement: TODO comments are present for suspicious logic
During stage extraction, any numerical constant or algorithmic invariant that cannot be directly derived from the PAL specification SHALL be annotated with a `// TODO: check the correctness: <reason>` comment in the production source.

#### Scenario: delta offset comment exists
- **WHEN** the `DecodePALSignal` code computing `skipUntil + delta` is reviewed
- **THEN** a `// TODO: check the correctness:` comment is present explaining the empirical delta calculation

#### Scenario: sync threshold MAD multiplier comment exists
- **WHEN** the `EstimateHorizontalOffset` or `FindFrameStart` source is reviewed
- **THEN** `// TODO: check the correctness:` comments are present next to the `2.5 * mad` and `3.0 * mad` heuristic thresholds explaining they are tuned values
