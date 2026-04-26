## ADDED Requirements

### Requirement: IqDemodulator stage
The `IqDemodulator` class SHALL read raw interleaved signed-int8 IQ bytes from a `Stream` (or accept a `FileStream` via `IQWavReader`) and produce a DC-removed `double[]` array of AM-magnitude video samples. It SHALL return `false` when the stream ends before the requested sample count is filled.

#### Scenario: File stream produces correct magnitude and removes DC
- **WHEN** a `FileStream` pointing to a WAV IQ file is passed with a known single-frequency tone
- **THEN** the output array length equals `samplesRequested`, all values are non-negative magnitudes, and the mean of the output array is approximately zero (DC removed)

#### Scenario: Generic stream (live HackRF) reads interleaved int8 IQ pairs
- **WHEN** a non-FileStream containing 2N bytes of interleaved signed-int8 (I,Q) pairs is passed with `samplesRequested = N`
- **THEN** the output array contains N magnitude values equal to `sqrt((I/128)²+(Q/128)²)` for each pair, and the array mean is approximately zero

#### Scenario: Short read returns false
- **WHEN** the stream contains fewer bytes than required to fill `samplesRequested` samples
- **THEN** `IqDemodulator.TryRead` returns `false`

---

### Requirement: FrameSynchronizer stage
The `FrameSynchronizer` class SHALL detect the start of a PAL video frame within a `double[]` video signal by locating the Vertical Blanking Interval (VBI) broad-pulse V-sync pattern, returning the sample index of the first active line.

#### Scenario: VBI detected in synthetic waveform
- **WHEN** a synthetic video signal is constructed with two consecutive broad pulses (≥ 20 µs at negative polarity) at a known position
- **THEN** `FrameSynchronizer.FindFrameStart` returns a sample index within ±1 line of the known VBI start

#### Scenario: Signal too short falls back to zero
- **WHEN** a video signal shorter than 50 lines is provided
- **THEN** `FrameSynchronizer.FindFrameStart` returns 0

#### Scenario: Inverted polarity is handled
- **WHEN** the video signal has the sync pulses at positive polarity instead of negative
- **THEN** auto-polarity detection inverts the signal and the result is the same as the non-inverted case

---

### Requirement: HorizontalAligner stage
The `HorizontalAligner` class SHALL estimate the horizontal active-pixel start offset by detecting H-sync pulses across a configurable number of lines and returning the median active-column index relative to the expected `(4.7 + 5.8) µs` blanking duration.

#### Scenario: H-sync detected across multiple lines
- **WHEN** a synthetic video signal is constructed with a 4.7 µs negative-polarity H-sync pulse per line at a known column offset
- **THEN** `HorizontalAligner.EstimateOffset` returns an integer within ±2 samples of the expected active-video start column

#### Scenario: No H-sync found returns zero
- **WHEN** a flat (constant) video buffer is provided
- **THEN** `HorizontalAligner.EstimateOffset` returns 0

---

### Requirement: FieldSplitter stage
The `FieldSplitter` class SHALL split a contiguous `double[]` video frame buffer into two contiguous field arrays using Strategy A (field 1 starts at line 0, field 2 starts at line `PAL_LINES_PER_FRAME / 2`).

#### Scenario: Correct field boundaries extracted
- **WHEN** a frame buffer filled with a known per-line marker value is split for a given `samplesPerLine`
- **THEN** `field1[0..samplesPerLine-1]` equals the first line of the source and `field2[0..samplesPerLine-1]` equals line 312 of the source

#### Scenario: Insufficient lines returns false
- **WHEN** the frame buffer is shorter than `field2StartLine + fieldLinesVisible` lines
- **THEN** `FieldSplitter.TrySplit` returns `false`

---

### Requirement: LumaChromaSeparator stage
The `LumaChromaSeparator` class SHALL apply a FIR low-pass filter for luminance and a FIR band-pass filter for chrominance to a video field `double[]` buffer in a single dual-filter pass, writing results into caller-supplied destination arrays.

#### Scenario: DC-only input passes through luma, suppressed in chroma
- **WHEN** a constant (DC) signal is passed through the separator
- **THEN** the luma output equals the DC value (after FIR settling) and the chroma output is near zero

#### Scenario: 4.43 MHz tone is passed to chroma, attenuated in luma
- **WHEN** a pure 4.43 MHz sine wave signal at the sample rate is passed through the separator
- **THEN** the chroma output amplitude is significantly larger than the luma output amplitude at that frequency

---

### Requirement: ChromaDecoder stage
The `ChromaDecoder` class SHALL demodulate a chrominance `double[]` buffer using a `ColorPll` instance and FIR low-pass filtering to produce U and V component arrays. On lines where V is inverted (odd absolute line number), it SHALL negate the demodulated V sample before filtering.

#### Scenario: V component inverted on odd lines
- **WHEN** a chrominance buffer representing two lines is decoded with `startLineOffset = 0`
- **THEN** the V output for samples on line 1 (index `samplesPerLine` to `2*samplesPerLine-1`) is the negation of the unmodified demod imaginary part

#### Scenario: U and V outputs are filtered
- **WHEN** a high-frequency noise burst is appended to a clean chroma signal and decoded
- **THEN** the LPF attenuates the noise burst in both U and V outputs

---

### Requirement: ActiveAreaCropper stage
The `ActiveAreaCropper` class SHALL copy the active-pixel columns of each line from a source signal array into a pre-allocated destination array, starting from column 0 for the given `copyWidth`.

#### Scenario: Correct columns extracted per line
- **WHEN** a signal buffer has each line filled with its line-index value and `copyWidth < samplesPerLine`
- **THEN** the destination buffer contains the first `copyWidth` samples of each line, with no data from the blanking interval columns

#### Scenario: Safe guard on short last line
- **WHEN** the last partial line would overflow the source buffer
- **THEN** `ActiveAreaCropper.CropInto` stops before copying the partial line without throwing

---

### Requirement: YuvToRgbConverter stage
The `YuvToRgbConverter` class SHALL convert Y, U, V double arrays into a `byte[,,]` RGB array using the BT.601 matrix. Output channel values SHALL be clamped to [0, 255].

#### Scenario: BT.601 white point (Y=1, U=0, V=0) maps to RGB (255, 255, 255)
- **WHEN** Y=1.0, U=0.0, V=0.0 are provided for all pixels
- **THEN** all RGB output bytes are 255

#### Scenario: BT.601 black point (Y=0, U=0, V=0) maps to RGB (0, 0, 0)
- **WHEN** Y=0.0, U=0.0, V=0.0 are provided
- **THEN** all RGB output bytes are 0

#### Scenario: Out-of-range values are clamped
- **WHEN** Y=1.5 is provided (above [0,1] range)
- **THEN** the R, G, B output bytes are 255 (clamped, not overflowing)

---

### Requirement: FieldInterleaver stage
The `FieldInterleaver` class SHALL weave two RGB field `byte[,,]` arrays into a full-frame `byte[,,]` according to the configured `FieldOrder`. In `BottomFieldFirst` mode the even output rows come from field 2 and odd rows from field 1; in `TopFieldFirst` mode this is reversed.

#### Scenario: BottomFieldFirst interleave order is correct
- **WHEN** field1 rows are filled with value 1 and field2 rows are filled with value 2, with `BottomFieldFirst`
- **THEN** output row 0 contains field2 pixels and output row 1 contains field1 pixels

#### Scenario: TopFieldFirst interleave order is correct
- **WHEN** field1 rows are filled with value 1 and field2 rows are filled with value 2, with `TopFieldFirst`
- **THEN** output row 0 contains field1 pixels and output row 1 contains field2 pixels

---

### Requirement: PALDecoder orchestrator preserves public API
The `PALDecoder` class SHALL continue to expose `DecodePALSignal(int sampleRate, Stream stream)` and `DecodePALSignal(int sampleRate, FileStream fs)` with identical behavior to the pre-refactor implementation.

#### Scenario: Entry point signature unchanged
- **WHEN** existing call sites compile against the refactored `PALDecoder`
- **THEN** no compilation errors occur (verified by `dotnet build`)
