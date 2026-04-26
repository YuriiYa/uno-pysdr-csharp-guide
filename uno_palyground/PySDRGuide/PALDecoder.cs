using ScottPlot;
using System.Diagnostics;
using Microsoft.UI.Dispatching;
using uno_palyground.PySDRGuide;
using System.Threading.Tasks;
using uno_palyground.PySDRGuide.Pipeline;

// ## Abbreviations:
// ADC: Analog-to-Digital Converter – digitizes incoming analog IF/baseband samples.
// DAC: Digital-to-Analog Converter – converts processed digital samples back to analog (for TX).
// DX: Long-Distance Reception/Communication – receiving signals from far-away stations (often weak).
// LNA: Low-Noise Amplifier – first RF stage, boosts weak signals with minimal added noise.  LNA: 0–40 dB in 8 dB steps.
// VGA: Variable Gain Amplifier – later (IF/baseband) gain stage for amplitude leveling. VGA: 0–62 dB in 2 dB steps
// AGC: Automatic Gain Control – dynamic adjustment loop managing overall receiver gain.
// RF: Radio Frequency – original high-frequency spectrum before downconversion.
// IF: Intermediate Frequency – shifted frequency after first mixing stage; easier to filter.
// LO: Local Oscillator – generates mixing frequency for up/down conversion.
// PLL: Phase-Locked Loop – locks oscillator phase/frequency to a reference (e.g., color burst, carrier).
// SNR: Signal-to-Noise Ratio – signal power relative to noise floor.
// NF: Noise Figure – measure of noise degradation added by a receiver chain.
// IMD: Intermodulation Distortion – spurious products from nonlinear mixing of strong signals.
// IP3: Third-Order Intercept Point – metric of linearity; higher implies better strong-signal handling.
// MDS: Minimum Discernible Signal – weakest signal detectable above noise.
// BW: Bandwidth – frequency span of interest or filter passband width.
// IQ: In-phase / Quadrature – orthogonal baseband components representing complex signals.
// DSP: Digital Signal Processing – algorithmic manipulation of sampled data.
// FIR: Finite Impulse Response (filter) – convolution-based filter with finite coefficients.
// FFT: Fast Fourier Transform – efficient computation of discrete spectrum.
// LO Leakage: Residual LO feedthrough at baseband (DC spur).
// DC Offset: Unwanted constant bias in I/Q channels (often corrected).
// AGC Loop: Control mechanism adjusting LNA/VGA for target amplitude.
// AFE: Analog Front End – collective RF input circuitry (filters, LNA, mixers).
// BB: Baseband – zero-IF representation after full downconversion.
// FOM: Figure of Merit – generalized performance composite metric (context-dependent).
// FEC: Forward Error Correction – coding to recover from bit errors (more for digital comms).
// BER: Bit Error Rate – fraction of received bits in error.
// QAM: Quadrature Amplitude Modulation – constellation using amplitude + phase states.
// AM: Amplitude Modulation – information encoded in amplitude variations.
// FM: Frequency Modulation – information encoded in instantaneous frequency shifts.
// SSB: Single Sideband – modulation using one spectral sideband.
// CW: Continuous Wave – unmodulated carrier (e.g., Morse transmissions).
// VCO: Voltage-Controlled Oscillator – tunable oscillator often inside PLL.
// PA: Power Amplifier – final high-power stage for transmission.
// DNR: Dynamic Noise Reduction (contextual DSP technique).
// TDD: Time Division Duplex – alternating TX/RX time slots (less Pluto-specific but SDR-relevant).

// PAL D type parameters:

// Video bandwidth: ~6.0 MHz (System D/K)
// FM sound offset: +6.5 MHz (for awareness; don’t band-pass near it)
// Luma LPF cutoff: ~5.0–5.5 MHz (I suggest 5.0 MHz for margin)
// Chroma BPF: center at 4.43361875 MHz, tighten to ~3.9–4.95 MHz
// Channel bandwidth: ~8.0 MHz (for awareness; don’t band-pass near it)


// 0 Hz    1 MHz    3 MHz   3.5 MHz 4.43 MHz  5.5 MHz    8 MHz
// |-------|--------|--------|--------|-------|----------|-------|
// |  DC   | Luma   | Guard  |     Chroma     |  Sound   | Rest  |
// |       | (Y)    | Band   |    (U & V)     |          |       |

// 0 Hz (DC Component)
// - What it is: The average value of the entire signal, representing the average brightness of the whole video frame.
// - Purpose: It sets the overall brightness level. In processing, this is almost always removed (as your AMDemodulation function does) to focus on the signal's variations, which contain the actual picture information.

// 0 Hz – ~3.0 MHz (Luminance Signal - Y)
// - What it is: This is the primary video signal, containing all the brightness and detail information. It is essentially the high-resolution black and white television picture.
// - Structure: A baseband analog signal where voltage level corresponds directly to pixel brightness. The lowest voltage represents the "blacker-than-black" sync pulses, and the highest voltage represents peak white.

// ~3.0 MHz – ~3.5 MHz (Guard Band)
// - What it is: An intentionally empty frequency range.
// - Purpose: To act as a buffer zone. It prevents the high-frequency components of the luminance signal (fine details) from interfering with the low-frequency components of the chrominance signal (color information). This is crucial for preventing "dot crawl" and other artifacts.

// ~3.5 MHz – ~5.5 MHz (Chrominance Signal - U & V)
// - What it is: This contains all the color information for the picture. It is cleverly "superimposed" on top of the luminance signal.
// - Structure: This is not a simple baseband signal. It is created using Quadrature Amplitude Modulation (QAM).
//   - Subcarrier: A precise sine wave at 4.43361875 MHz acts as the "carrier" for the color.
//   - U Component (B-Y): The "blue difference" signal. It is modulated onto the subcarrier.
//   - V Component (R-Y): The "red difference" signal. It is also modulated onto the subcarrier, but 90 degrees out of phase with the U component.
//   - PAL Specifics (Phase Alternating Line): This is the key feature of PAL. The phase of the V component is inverted (+90° to -90°) on every other line. This clever trick allows the receiver to average out and cancel any phase errors that occur during transmission, resulting in more stable color than the NTSC system.

// ~5.5 MHz – ~6.5 MHz (Sound Carrier & Guard Band)
// - What it is: The audio signal for the broadcast.
// - Structure: The audio is frequency modulated (FM) onto its own carrier, typically located at 5.5 MHz or 6.0 MHz above the main video carrier. This is why you can filter out the video without affecting the sound.

// ~6.5 MHz – 8 MHz (Rest / Unused)
// - What it is: The remainder of the allocated 8 MHz channel.
// - Purpose: This is another guard band to prevent interference with the next adjacent television channel. Your SDR captures this, but it contains no useful information for decoding this specific PAL signal.

// Detailed Time Domain Structure (One Video Line - 64µs)
// This is the most critical part for understanding your frame alignment problem. This describes the structure of the signal over time for a single horizontal line.

// <--1.5µs-->|<----4.7µs---->|<----5.8µs---->|<------------------ 52µs ------------------>
//            |               |              |
//   Front    |   H-Sync      |  Back Porch  |                 Active Video
//   Porch    |   Pulse       |              |
// -----------|---------------|--------------|-----------------------------------------------
//    0V      |               |    0V        |  Luma (0V to +0.7V) + Chroma (Subcarrier)
//            |     -0.3V     |              |
//            |_______________|--Breezeway---|
//                                |--Burst---|

// Front Porch (1.5 µs)

// - Voltage: 0V (Black Level).
// - Purpose: A brief, stable period of black level that ensures the previous line's active video has finished before the critical sync pulse arrives. It allows the receiver's circuitry to prepare for the sync.

// Horizontal Sync Pulse (4.7 µs)
// - Voltage: -0.3V ("Blacker-than-black").
// - Purpose: This is the master timing pulse for the horizontal line. Its sharp leading edge is the primary trigger that tells the receiver to end the current line and begin the horizontal retrace (flyback). This is the pulse your line-synchronization code should lock onto.

// Back Porch (5.8 µs)
// - This is the most complex part of the blanking interval.
// - Breezeway (~0.6 µs): A very short gap at 0V immediately after the sync pulse. It prevents the sync pulse from interfering with the color burst.
// - Color Burst (2.25 µs): This is the color reference key.
//   - Content: It consists of 10±1 cycles of the pure color subcarrier (4.43361875 MHz).
//   - Purpose: It provides a reference for both the phase and frequency of the color signal. The receiver uses a Phase-Locked Loop (PLL) to lock onto this burst. All color information in the active video is then measured relative to the phase of this burst.
//   - PAL V-Axis Information: The phase of the burst itself is shifted to signal the state of the V-component for the upcoming line.
//     - On lines where V is positive, the burst's phase is 135°.
//     - On lines where V is negative, the burst's phase is 225°.

// Active Video (52 µs)
// - Purpose: This is the only part of the signal that contains the visible picture information for that line.
// - Structure: It's a composite signal made of two parts added together:
//   - Luminance (Y): The baseband brightness signal. Its voltage varies between 0V (black) and +0.7V (peak white), drawing the details of the image.
//   - Chrominance (C): The color information, which is a 4.43 MHz sine wave "riding" on top of the luminance signal.
//     - Saturation is determined by the amplitude (size) of this sine wave. No amplitude means no color (grayscale).
//     - Hue is determined by the phase of this sine wave relative to the reference color burst.
//     - V-Axis Switch: On one line, the color is represented by U + V. On the very next line, it's represented by U - V. This is the "Phase Alternating Line" from which PAL gets its name.

// Level 2: The PAL 8-Field Sequence Explained
// This sequence arises from three interacting cycles that only align every 8 fields (4 frames).

// 25 fps / 50 fields: Standard interlaced video.
// 625 lines: An odd number, meaning each field has 312.5 lines. This causes the second field of a frame to start halfway across the screen, creating the interlace.
// Quarter-line offset of f_sc: The color subcarrier frequency is chosen to be (283.75) * f_line + 25 Hz. The crucial part is the .75, which means the subcarrier phase shifts by -90° (or +270°) from the start of one line to the start of the next.
// V-Axis Switch: The + V / -V switch happens every line.
// The combination of these factors means it takes 4 full frames (8 fields) for the exact same line-start/subcarrier-phase/V-switch relationship to occur on line 1 of the frame.

// Here is the detailed 8-field structure:
// Field Frame	Lines in Field	V-Switch on Line 1	Subcarrier Phase on Line 1	Notes
// 1	    1	    1 to 312.5	        +V	                0°	                    Start of the sequence.
// 2	    1	    313 to 625	        -V	                180°	                Field starts on a half-line. Subcarrier is inverted.
// 3	    2	    1 to 312.5	        -V	                0°	                    V-switch state is now opposite to Field 1.
// 4	    2	    313 to 625	        +V	                180°	                V-switch is opposite to Field 2.
// 5	    3	    1 to 312.5	        +V	                180°	                This is the key change. Subcarrier is now inverted compared to Field 1.
// 6	    3	    313 to 625	        -V	                0°	                    Subcarrier is inverted compared to Field 2.
// 7	    4	    1 to 312.5	        -V	                180°	                V-switch is opposite, subcarrier is inverted compared to Field 1.
// 8	    4	    313 to 625	        +V	                0°	                    V-switch is opposite, subcarrier is inverted compared to Field 2.
// 9	    5	    1 to 312.5	        +V	                0°	                    Sequence Repeats. Identical to Field 1.

// This complex relationship is why professional PAL editing equipment had to be "8-field aware" to ensure edits were seamless and didn't disrupt the color sequence. For your decoder, while you don't need to implement all of this, understanding it explains why simple frame-by-frame processing can sometimes lead to color or timing artifacts.
// Add system profiles above the PALDecoder class or inside it (make static if inside)

// Each field has several lines that do not contain any picture information. These lines occur during the Vertical Blanking Interval (VBI) when the electron beam in the TV tube returns from the bottom to the top of the screen. The VBI also covers several lines which carry data such as Closed Captions and Teletext.
// ITU-R BT.601-5 and 656-4 describe a digital active area. This is used when the analogue video signal is converted to a digital format, and it does not exactly coincide with the analogue active area.
// detailed explanation of the format https://dvmp.co.uk/digital-video.htm
// TvSystem, SystemProfile, FieldOrder are defined in Pipeline/PalTypes.cs

public class PALDecoder
{
    private readonly Plot _plot;
    private readonly Plot _plotModulation;
    private readonly DispatcherQueue _dispatcherQueue;
    private readonly SystemProfile _profile;
    private readonly FieldOrder _fieldOrder;

    // PAL I specifications — see PalConstants for authoritative values

    // Tunable FIR tap lengths — see PalConstants.LumaLpfTaps / ChromaSeparationTaps / ChromaBasebandLpfTaps

    public PALDecoder(Plot plot, Plot plotModulation, DispatcherQueue dispatcherQueue, TvSystem system = TvSystem.PAL_DK, FieldOrder fieldOrder = FieldOrder.BottomFieldFirst)
    {
        _plot = plot;
        _plotModulation = plotModulation;
        _dispatcherQueue = dispatcherQueue;
        _profile = SystemProfile.For(system);
        _fieldOrder = fieldOrder;
    }

    // Reusable buffers to minimize per-frame allocations
    private double[]? _field1Buffer;
    private double[]? _field2Buffer;
    private byte[,,]? _rgbField1Buffer;
    private byte[,,]? _rgbField2Buffer;
    private byte[,,]? _rgbInterleavedFrameBuffer;
    // Active-area reusable buffers (post-crop) per component per field
    private double[]? _lum1Active, _lum2Active;
    private double[]? _u1Active, _u2Active;
    private double[]? _v1Active, _v2Active;
    // Per-field parallel buffers (avoid contention when processing the two fields concurrently)
    private double[]? _luminanceBufferF1;
    private double[]? _chrominanceBufferF1;
    private double[]? _luminanceBufferF2;
    private double[]? _chrominanceBufferF2;
    private double[]? _uFilteredF1;
    private double[]? _vFilteredF1;
    private double[]? _uFilteredF2;
    private double[]? _vFilteredF2;
    // Pipeline stages (lazy-initialized per sample rate)
    private IqDemodulator? _iqDemodulator;
    private FrameSynchronizer? _frameSynchronizer;
    private HorizontalAligner? _horizontalAligner;
    private FieldSplitter? _fieldSplitter;
    private LumaChromaSeparator? _lumaSeparatorF1;
    private LumaChromaSeparator? _lumaSeparatorF2;
    private ChromaDecoder? _chromaDecoderF1;
    private ChromaDecoder? _chromaDecoderF2;
    private ActiveAreaCropper? _activeAreaCropper;
    private YuvToRgbConverter? _yuvToRgbConverter;
    private FieldInterleaver? _fieldInterleaver;

    // Profiling accumulators (ticks)
    private long _ticksFieldCopy;
    private long _ticksLumaChroma;
    private long _ticksChromaDecode;
    private long _ticksCrop;
    private long _ticksRgbConvert;
    private long _ticksInterleave;
    private int _profiledFrames;
    private const bool EnableStageProfiling = true; // toggle to disable detailed stage timing

    private static void EnsureBuffer(ref double[]? buffer, int requiredLength)
    {
        if (buffer == null || buffer.Length < requiredLength)
            buffer = new double[requiredLength];
    }

    private static void EnsureRgbBuffer(ref byte[,,]? buffer, int height, int width)
    {
        if (buffer == null || buffer.GetLength(0) != height || buffer.GetLength(1) != width)
            buffer = new byte[height, width, 3];
    }

    // Legacy entry point retained for convenience; delegates to Stream overload.
    public void DecodePALSignal(int sampleRate, FileStream fs) => DecodePALSignal(sampleRate, (Stream)fs);

    public void DecodePALSignal(int sampleRate, Stream stream)
    {
        int samplesPerLine = (int)(PalConstants.PAL_LINE_DURATION * sampleRate);
        int samplesPerFrame = samplesPerLine * PalConstants.PAL_LINES_PER_FRAME;
        bool isFile = stream is FileStream;
        var file = stream as FileStream;

        if (isFile)
        {
            try { IQWavReader.ConfigureRingBuffer(samplesPerFrame * 2); } catch { /* ignore */ }
        }

        // Lazy-initialize pipeline stages
        _iqDemodulator ??= new IqDemodulator();
        _frameSynchronizer ??= new FrameSynchronizer();
        _horizontalAligner ??= new HorizontalAligner();
        _fieldSplitter ??= new FieldSplitter();
        _activeAreaCropper ??= new ActiveAreaCropper();
        _yuvToRgbConverter ??= new YuvToRgbConverter();
        _fieldInterleaver ??= new FieldInterleaver();

        double[] videoSignal = new double[samplesPerFrame];
        if (!_iqDemodulator.TryRead(stream, samplesPerFrame, videoSignal))
        {
            Console.WriteLine("No samples available for PAL decoding (initial frame).");
            return;
        }
        int frameStart = _frameSynchronizer.FindFrameStart(videoSignal, sampleRate, samplesPerLine);
        int autoHOffset = _horizontalAligner.EstimateOffset(videoSignal, frameStart, samplesPerLine, sampleRate);
        int skipUntil = frameStart + autoHOffset;

        var nonVideoData = (int)Math.Round((1.5 + 4.7 + 5.8) * 1e-6 * sampleRate);
        // TODO: check the correctness: delta calculation mixes blanking intervals in a non-obvious way; verify against actual broadcast timing
        var delta = (samplesPerLine - nonVideoData) / 2 + nonVideoData;
        if (skipUntil < 0) skipUntil = 0;
        if (skipUntil > samplesPerFrame) skipUntil = samplesPerFrame;

        SkipSamplesStreamingGeneric(stream, skipUntil + delta);

        long numberOfFrames;
        if (isFile)
        {
            long remainingComplexSamples = (file!.Length - file.Position) / 4;
            numberOfFrames = remainingComplexSamples / samplesPerFrame;
            if (numberOfFrames <= 0)
            {
                Console.WriteLine("No complete frames remain after initial sync frame.");
                return;
            }
            Console.WriteLine($"Decoding {numberOfFrames} frames (samplesPerLine={samplesPerLine}, samplesPerFrame={samplesPerFrame}).");
        }
        else
        {
            numberOfFrames = long.MaxValue;
            Console.WriteLine("Decoding continuous stream (unknown frame count)...");
        }

        double[] frameData = new double[samplesPerFrame];
        int preActiveWidth = (int)Math.Round(52e-6 * sampleRate);
        int preDesiredActiveCol = (int)Math.Round((4.7 + 5.8) * 1e-6 * sampleRate);
        int copyWidth = Math.Max(preActiveWidth, Math.Max(0, samplesPerLine - preDesiredActiveCol));
        int activeLinesPerField = PalConstants.PAL_VISIBLE_LINES / 2; // 288

        int activeBufferLenPerField = activeLinesPerField * copyWidth;
        EnsureBuffer(ref _lum1Active, activeBufferLenPerField);
        EnsureBuffer(ref _lum2Active, activeBufferLenPerField);
        EnsureBuffer(ref _u1Active, activeBufferLenPerField);
        EnsureBuffer(ref _u2Active, activeBufferLenPerField);
        EnsureBuffer(ref _v1Active, activeBufferLenPerField);
        EnsureBuffer(ref _v2Active, activeBufferLenPerField);

        for (long frameIndex = 0; frameIndex < numberOfFrames; frameIndex++)
        {
            byte[,,]? rgbFrame = null;
            var stateResult = TimeLapseHelper.PrintTime(() =>
            {
                Console.WriteLine($"Decoding frame {frameIndex + 1}/{numberOfFrames}...");
                long frameFieldCopyTicks = 0, frameLumaChromaTicks = 0, frameChromaDecodeTicks = 0, frameCropTicks = 0, frameRgbTicks = 0, frameInterleaveTicks = 0;

                if (!isFile && stream is HackRF.Namespace.TeeIqStream tee)
                {
                    tee.DrainToLatestFrame(samplesPerFrame * 2);
                }
                if (!_iqDemodulator.TryRead(stream, samplesPerFrame, frameData))
                {
                    Console.WriteLine($"Short read at frame {frameIndex}; expected {samplesPerFrame} samples. Stopping.");
                    return (breakResult: true, continueResult: false);
                }

                long t0 = Stopwatch.GetTimestamp();
                int fieldSamples = activeLinesPerField * samplesPerLine;
                EnsureBuffer(ref _field1Buffer, fieldSamples);
                EnsureBuffer(ref _field2Buffer, fieldSamples);
                if (!_fieldSplitter!.TrySplit(frameData, samplesPerLine, _field1Buffer!, _field2Buffer!))
                {
                    Console.WriteLine($"Frame {frameIndex}: insufficient lines for field split. Skipping frame.");
                    return (breakResult: false, continueResult: true);
                }
                long t1 = Stopwatch.GetTimestamp();
                frameFieldCopyTicks = t1 - t0;

                var field1 = _field1Buffer!;
                var field2 = _field2Buffer!;
                double[]? lum1 = null, u1 = null, v1 = null;
                double[]? lum2 = null, u2 = null, v2 = null;

                // Lazy-initialize per-field separator and decoder pairs
                _lumaSeparatorF1 ??= new LumaChromaSeparator(_profile, sampleRate, PalConstants.LumaLpfTaps, PalConstants.ChromaSeparationTaps);
                _lumaSeparatorF2 ??= new LumaChromaSeparator(_profile, sampleRate, PalConstants.LumaLpfTaps, PalConstants.ChromaSeparationTaps);
                _chromaDecoderF1 ??= new ChromaDecoder(sampleRate, PalConstants.ChromaBasebandLpfTaps);
                _chromaDecoderF2 ??= new ChromaDecoder(sampleRate, PalConstants.ChromaBasebandLpfTaps);

                t0 = Stopwatch.GetTimestamp();
                Parallel.Invoke(
                    () =>
                    {
                        int len1 = field1.Length;
                        EnsureBuffer(ref _luminanceBufferF1, len1);
                        EnsureBuffer(ref _chrominanceBufferF1, len1);
                        _lumaSeparatorF1!.Separate(field1, _luminanceBufferF1!, _chrominanceBufferF1!);
                        lum1 = _luminanceBufferF1!;
                        EnsureBuffer(ref _uFilteredF1, len1);
                        EnsureBuffer(ref _vFilteredF1, len1);
                        _chromaDecoderF1!.Decode(_chrominanceBufferF1!, sampleRate, samplesPerLine, 0, _uFilteredF1!, _vFilteredF1!);
                        u1 = _uFilteredF1!; v1 = _vFilteredF1!;
                    },
                    () =>
                    {
                        int len2 = field2.Length;
                        EnsureBuffer(ref _luminanceBufferF2, len2);
                        EnsureBuffer(ref _chrominanceBufferF2, len2);
                        _lumaSeparatorF2!.Separate(field2, _luminanceBufferF2!, _chrominanceBufferF2!);
                        lum2 = _luminanceBufferF2!;
                        EnsureBuffer(ref _uFilteredF2, len2);
                        EnsureBuffer(ref _vFilteredF2, len2);
                        _chromaDecoderF2!.Decode(_chrominanceBufferF2!, sampleRate, samplesPerLine, PalConstants.PAL_LINES_PER_FRAME / 2, _uFilteredF2!, _vFilteredF2!);
                        u2 = _uFilteredF2!; v2 = _vFilteredF2!;
                    }
                );
                t1 = Stopwatch.GetTimestamp();
                long combined = t1 - t0;
                double prevSep = _ticksLumaChroma > 0 && _ticksChromaDecode > 0 ? _ticksLumaChroma : 6;
                double prevChroma = _ticksChromaDecode > 0 ? _ticksChromaDecode : 4;
                double totalPrev = prevSep + prevChroma;
                frameLumaChromaTicks = (long)(combined * (prevSep / totalPrev));
                frameChromaDecodeTicks = combined - frameLumaChromaTicks;

                t0 = Stopwatch.GetTimestamp();
                _activeAreaCropper!.CropInto(lum1!, samplesPerLine, _lum1Active!, copyWidth);
                _activeAreaCropper!.CropInto(u1!, samplesPerLine, _u1Active!, copyWidth);
                _activeAreaCropper!.CropInto(v1!, samplesPerLine, _v1Active!, copyWidth);
                _activeAreaCropper!.CropInto(lum2!, samplesPerLine, _lum2Active!, copyWidth);
                _activeAreaCropper!.CropInto(u2!, samplesPerLine, _u2Active!, copyWidth);
                _activeAreaCropper!.CropInto(v2!, samplesPerLine, _v2Active!, copyWidth);
                t1 = Stopwatch.GetTimestamp();
                frameCropTicks = t1 - t0;
                int activeWidth = copyWidth;

                int fieldHeight = activeLinesPerField;
                EnsureRgbBuffer(ref _rgbField1Buffer, fieldHeight, activeWidth);
                EnsureRgbBuffer(ref _rgbField2Buffer, fieldHeight, activeWidth);
                var rgbField1 = _rgbField1Buffer!;
                var rgbField2 = _rgbField2Buffer!;
                t0 = Stopwatch.GetTimestamp();
                _yuvToRgbConverter!.Convert(_lum1Active!, _u1Active!, _v1Active!, activeWidth, rgbField1);
                _yuvToRgbConverter!.Convert(_lum2Active!, _u2Active!, _v2Active!, activeWidth, rgbField2);
                t1 = Stopwatch.GetTimestamp();
                frameRgbTicks = t1 - t0;

                int finalHeight = Math.Min(PalConstants.PAL_VISIBLE_LINES, fieldHeight * 2);
                EnsureRgbBuffer(ref _rgbInterleavedFrameBuffer, finalHeight, activeWidth);
                t0 = Stopwatch.GetTimestamp();
                _fieldInterleaver!.Interleave(rgbField1, rgbField2, _rgbInterleavedFrameBuffer!, _fieldOrder);
                t1 = Stopwatch.GetTimestamp();
                frameInterleaveTicks = t1 - t0;
                rgbFrame = _rgbInterleavedFrameBuffer;

                if (EnableStageProfiling)
                {
                    _ticksFieldCopy += frameFieldCopyTicks;
                    _ticksLumaChroma += frameLumaChromaTicks;
                    _ticksChromaDecode += frameChromaDecodeTicks;
                    _ticksCrop += frameCropTicks;
                    _ticksRgbConvert += frameRgbTicks;
                    _ticksInterleave += frameInterleaveTicks;
                    _profiledFrames++;
                    double ToMs(long ticks) => ticks * 1000.0 / Stopwatch.Frequency;
                    Console.WriteLine($"Stage ms (frame {frameIndex}): copy={ToMs(frameFieldCopyTicks):F3} luma/chroma={ToMs(frameLumaChromaTicks):F3} chromaDecode={ToMs(frameChromaDecodeTicks):F3} crop={ToMs(frameCropTicks):F3} rgb={ToMs(frameRgbTicks):F3} interleave={ToMs(frameInterleaveTicks):F3}");
                    if (_profiledFrames % 10 == 0)
                    {
                        Console.WriteLine($"Averages over {_profiledFrames} frames: copy={ToMs(_ticksFieldCopy / _profiledFrames):F3} luma/chroma={ToMs(_ticksLumaChroma / _profiledFrames):F3} chromaDecode={ToMs(_ticksChromaDecode / _profiledFrames):F3} crop={ToMs(_ticksCrop / _profiledFrames):F3} rgb={ToMs(_ticksRgbConvert / _profiledFrames):F3} interleave={ToMs(_ticksInterleave / _profiledFrames):F3}");
                    }
                }
                return (breakResult: false, continueResult: false);
            });
            if (stateResult.breakResult) break;
            if (stateResult.continueResult) continue;

            if (rgbFrame == null)
            {
                Console.WriteLine($"Frame {frameIndex}: failed to interleave fields.");
                continue;
            }

            DisplayVideoFrame(rgbFrame);

            if (!isFile && rgbFrame == null)
            {
                Console.WriteLine("Continuous stream ended.");
                break;
            }
        }
    }

    private static void SkipSamplesStreaming(FileStream fs, int samplesToSkip)
    {
        int remaining = samplesToSkip;
        while (remaining > 0)
        {
            int chunk = Math.Min(remaining, 8192);
            var seg = IQWavReader.ReadIQIntoRingOptimized(fs, chunk);
            if (seg.IsEmpty) break; // EOF
            remaining -= seg.SamplesRead;
        }
    }

    // Generic skip: discard raw bytes for non-FileStream sources (assumes 2 bytes per complex sample int8 I,Q).
    private static void SkipSamplesStreamingGeneric(Stream stream, int samplesToSkip)
    {
        if (samplesToSkip <= 0) return;
        if (stream is FileStream fs)
        {
            SkipSamplesStreaming(fs, samplesToSkip);
            return;
        }
        const int bytesPerSample = 2;
        long bytesToDiscard = (long)samplesToSkip * bytesPerSample;
        byte[] buffer = new byte[Math.Min(bytesToDiscard, 8192)];
        while (bytesToDiscard > 0)
        {
            int want = (int)Math.Min(bytesToDiscard, buffer.Length);
            int read = stream.Read(buffer, 0, want);
            if (read <= 0) break;
            bytesToDiscard -= read;
        }
    }


    private void DisplayVideoFrame(byte[,,] rgbFrame)
    {
        _dispatcherQueue.TryEnqueue(() =>
       {
           int height = rgbFrame.GetLength(0);
           int width = rgbFrame.GetLength(1);

           // Convert to format suitable for ScottPlot
           double[,] grayScale = new double[height, width];

           for (int row = 0; row < height; row++)
           {
               for (int col = 0; col < width; col++)
               {
                   // Use the original order (no horizontal flip)
                   double gray = 0.299 * rgbFrame[row, col, 0] +
                                 0.587 * rgbFrame[row, col, 1] +
                                 0.114 * rgbFrame[row, col, 2];
                   grayScale[row, col] = gray;
               }
           }

           _plot.Clear();
           var hm = _plot.Add.Heatmap(grayScale);
           hm.Colormap = new ScottPlot.Colormaps.Grayscale();

           _plot.XLabel("Pixel");
           _plot.YLabel("Line");
           _plot.Title("PAL Video Frame");
           _plot.Axes.SetLimitsX(0, width);
           _plot.Axes.SetLimitsY(0, height);
           

           _plot.Axes.Link(_plotModulation, x: true, y: false);
           _plotModulation.PlotControl?.Refresh();
           _plot.PlotControl?.Refresh();

           Console.WriteLine($"Decoded PAL frame: {width}x{height} pixels");
       });

    }
}
