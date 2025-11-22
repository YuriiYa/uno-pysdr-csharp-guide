using ScottPlot;
using System.Numerics;
using System.Diagnostics;
using System.Runtime.CompilerServices;
using Microsoft.UI.Dispatching;
using uno_palyground.PySDRGuide;
using System.Threading.Tasks;
using System.Collections.Concurrent;

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
public enum TvSystem { PAL_I, PAL_DK }

public readonly record struct SystemProfile(
    double LumaCutoffHz,
    double ChromaLowHz,
    double ChromaHighHz,
    double AudioOffsetHz,
    double VideoBandwidthHz
)
{
    public static SystemProfile For(TvSystem s) => s switch
    {
        TvSystem.PAL_DK => new(
            LumaCutoffHz: 5.0e6,        // keep luma safely below 6 MHz video BW
            ChromaLowHz: 3.9e6,        // tighter around 4.4336 MHz
            ChromaHighHz: 4.95e6,
            AudioOffsetHz: 6.5e6,       // FM sound offset (awareness)
            VideoBandwidthHz: 6.0e6
        ),
        TvSystem.PAL_I => new(
            LumaCutoffHz: 4.8e6,
            ChromaLowHz: 3.9e6,
            ChromaHighHz: 4.95e6,
            AudioOffsetHz: 6.0e6,
            VideoBandwidthHz: 5.5e6
        ),
        _ => throw new ArgumentOutOfRangeException(nameof(s))
    };
}

public enum FieldOrder
{
    TopFieldFirst,
    BottomFieldFirst
}

public class PALDecoder
{
    private readonly Plot _plot;
    private readonly Plot _plotModulation;
    private readonly DispatcherQueue _dispatcherQueue;
    private readonly SystemProfile _profile;
    private readonly FieldOrder _fieldOrder;

    // PAL I specifications
    public const double PAL_FRAME_RATE = 25.0; // 25 fps
    public const int PAL_LINES_PER_FRAME = 625;
    public const int PAL_VISIBLE_LINES = 576;
    public const double PAL_LINE_DURATION = 64e-6; // 64 microseconds per line
    public const double PAL_COLOR_CARRIER_FREQ = 4433618.75; // Hz
    public const double PAL_VIDEO_BANDWIDTH = 5.5e6; // 5.5 MHz

    // Tunable FIR tap lengths (performance vs quality)
    // Use the SAME tap length for luma & chroma separation so we can exploit the single-pass dual filter.
    // If you observe luma/chroma leakage or softness, raise to 91 or 101; if performance is still insufficient,
    // you can try lowering to 73. Keeping them equal is key for ApplyTwoFiltersStaticToDest to avoid two passes.
    private const int LUMA_LPF_TAPS = 81;           // unified length (was 101)
    private const int CHROMA_SEPARATION_TAPS = 81;  // unified length (was 81 already)
    private const int CHROMA_BASEBAND_LPF_TAPS = 63; // after demod low-pass; can be shorter

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
    // Filter / color separation reusable buffers
    private double[]? _luminanceBuffer;
    private double[]? _chrominanceBuffer;
    private double[]? _uScratch;
    private double[]? _vScratch;
    private double[]? _uFiltered;
    private double[]? _vFiltered;
    // Per-field parallel buffers (avoid contention when processing the two fields concurrently)
    private double[]? _luminanceBufferF1;
    private double[]? _chrominanceBufferF1;
    private double[]? _luminanceBufferF2;
    private double[]? _chrominanceBufferF2;
    private double[]? _uScratchF1;
    private double[]? _vScratchF1;
    private double[]? _uFilteredF1;
    private double[]? _vFilteredF1;
    private double[]? _uScratchF2;
    private double[]? _vScratchF2;
    private double[]? _uFilteredF2;
    private double[]? _vFilteredF2;
    // Reversed filter cache for SIMD FIR (stores reversed copies to enable forward vector dot products)
    private static class FilterReverseCache
    {
        // ConditionalWeakTable is thread-safe but Add can race if we do a TryGet+Add pattern.
        // Use GetValue factory which is atomic per key.
        private static readonly ConditionalWeakTable<double[], double[]> Reverse = new();
        public static double[] Get(double[] filter)
        {
            return Reverse.GetValue(filter, static f =>
            {
                var r = (double[])f.Clone();
                Array.Reverse(r);
                return r;
            });
        }
    }
    private static double[] GetReversedFilter(double[] filter) => FilterReverseCache.Get(filter);

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

    // New generic Stream variant. If the stream is not a FileStream (e.g. live HackRF),
    // frame count is unknown and decoding proceeds until reads fail.
    public void DecodePALSignal(int sampleRate, Stream stream)
    {
        int samplesPerLine = (int)(PAL_LINE_DURATION * sampleRate);
        int samplesPerFrame = samplesPerLine * PAL_LINES_PER_FRAME;
        bool isFile = stream is FileStream;
        var file = stream as FileStream;

        if (isFile)
        {
            // Configure ring buffer only for file-based replay.
            try { IQWavReader.ConfigureRingBuffer(samplesPerFrame * 2); } catch { /* ignore */ }
        }

        double[] videoSignal = new double[samplesPerFrame];
        if (!ReadAndDemodFrame(stream, samplesPerFrame, videoSignal))
        {
            Console.WriteLine("No samples available for PAL decoding (initial frame).");
            return;
        }
        int frameStart = FindFrameStart(videoSignal, sampleRate, samplesPerLine);
        int autoHOffset = EstimateHorizontalOffset(videoSignal, frameStart, samplesPerLine, sampleRate);
        int skipUntil = frameStart + autoHOffset;

        var nonVideoData = (int)Math.Round((1.5 + 4.7 + 5.8) * 1e-6 * sampleRate);
        var delta = (samplesPerLine - nonVideoData) / 2 + nonVideoData;
        if (skipUntil < 0) skipUntil = 0;
        if (skipUntil > samplesPerFrame) skipUntil = samplesPerFrame;

        // Advance / discard alignment samples.
        SkipSamplesStreamingGeneric(stream, skipUntil+delta);

        long numberOfFrames;
        if (isFile)
        {
            long remainingComplexSamples = (file!.Length - file.Position) / 4; // assumes 4 bytes per complex sample in file path
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
            numberOfFrames = long.MaxValue; // continuous
            Console.WriteLine("Decoding continuous stream (unknown frame count)...");
        }

        double[] frameData = new double[samplesPerFrame];
        // Precompute active-width parameters (constant across frames for given sampleRate)
        int preActiveWidth = (int)Math.Round(52e-6 * sampleRate);
        int preDesiredActiveCol = (int)Math.Round((4.7 + 5.8) * 1e-6 * sampleRate);
        // replicate original CropToActive copyWidth logic
        int copyWidth = Math.Max(preActiveWidth, Math.Max(0, samplesPerLine - preDesiredActiveCol));
        int activeLinesPerField = PAL_VISIBLE_LINES / 2; // 288

        // Ensure active buffers once (size: lines * copyWidth)
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
                // For continuous live sources, discard backlog to minimize latency.
                if (!isFile && stream is HackRF.Namespace.TeeIqStream tee)
                {
                    tee.DrainToLatestFrame(samplesPerFrame * 2); // raw bytes (I,Q int8) per frame
                }
                if (!ReadAndDemodFrame(stream, samplesPerFrame, frameData))
                {
                    Console.WriteLine($"Short read at frame {frameIndex}; expected {samplesPerFrame} samples. Stopping.");
                    return (breakResult: true, continueResult: false);
                }

                // Separate fields
                int fieldLines = PAL_VISIBLE_LINES / 2; // 288
                                                        // After you fill frameData with exactly one PAL frame (625 lines worth) OR one visible frame window,
                                                        // split fields as contiguous blocks, not every-other-line.
                int linesPerFieldAll = PAL_LINES_PER_FRAME / 2;   // 312
                int fieldLinesVis = PAL_VISIBLE_LINES / 2;        // 288

                // If frameStart already points to Field-1 active top, use 0 and +linesPerFieldAll.
                // Otherwise add per-field VBI offsets (~22–25 lines) to land on active video:
                // const int VBI_LINES_FIELD1 = 23;
                // const int VBI_LINES_FIELD2 = 23;

                // Choose one of the two strategies:

                // Strategy A: frameData starts at Field-1 active line 0
                int field1StartLine = 0;
                int field2StartLine = linesPerFieldAll;

                // Strategy B: frameData starts at the very beginning of the frame (includes VBI)
                //int field1StartLine = VBI_LINES_FIELD1;
                //int field2StartLine = linesPerFieldAll + VBI_LINES_FIELD2;

                int fieldSamples = fieldLines * samplesPerLine; // contiguous samples per field (visible portion)
                EnsureBuffer(ref _field1Buffer, fieldSamples);
                EnsureBuffer(ref _field2Buffer, fieldSamples);
                var field1 = _field1Buffer!;
                var field2 = _field2Buffer!;

                int availableLines = frameData.Length / samplesPerLine;
                if (availableLines < (field2StartLine + fieldLinesVis))
                {
                    Console.WriteLine($"Frame {frameIndex}: insufficient lines (available={availableLines}) for expected second field start {field2StartLine}. Skipping frame.");
                    return (breakResult: false, continueResult: true);
                }
                // Copy each field as a single contiguous block (eliminates per-line loop & allocations)
                // Visible lines are already contiguous for each field given strategy A
                int srcField1OffsetSamples = field1StartLine * samplesPerLine;
                int srcField2OffsetSamples = field2StartLine * samplesPerLine;
                int copySamples = fieldLinesVis * samplesPerLine;
                int bytesToCopy = copySamples * sizeof(double);
                long t0 = Stopwatch.GetTimestamp();
                if (srcField1OffsetSamples + copySamples <= frameData.Length)
                    Buffer.BlockCopy(frameData, srcField1OffsetSamples * sizeof(double), field1, 0, bytesToCopy);
                if (srcField2OffsetSamples + copySamples <= frameData.Length)
                    Buffer.BlockCopy(frameData, srcField2OffsetSamples * sizeof(double), field2, 0, bytesToCopy);
                long t1 = Stopwatch.GetTimestamp();
                frameFieldCopyTicks = t1 - t0;
                // Pre-warm filters once (outside parallel region) to avoid dictionary races
                double lumaCutoff = Math.Min(_profile.LumaCutoffHz, 0.45 * sampleRate);
                var lumaFilter = CreateLowPassFilter(lumaCutoff, sampleRate, LUMA_LPF_TAPS);
                double chromaLow = _profile.ChromaLowHz;
                double chromaHigh = _profile.ChromaHighHz;
                var chromaFilter = CreateBandPassFilter(chromaLow, chromaHigh, sampleRate, CHROMA_SEPARATION_TAPS);
                var chromaLPF = CreateLowPassFilter(1.3e6, sampleRate, CHROMA_BASEBAND_LPF_TAPS);

                // Process each field in parallel: luma/chroma separation + chroma decode
                double[]? lum1 = null, chr1 = null, u1 = null, v1 = null;
                double[]? lum2 = null, chr2 = null, u2 = null, v2 = null;

                t0 = Stopwatch.GetTimestamp();
                Parallel.Invoke(
                    () =>
                    {
                        // Field 1 luma/chroma
                        int len1 = field1.Length;
                        EnsureBuffer(ref _luminanceBufferF1, len1);
                        EnsureBuffer(ref _chrominanceBufferF1, len1);
                        ApplyTwoFiltersStaticToDest(field1, lumaFilter, chromaFilter, _luminanceBufferF1!, _chrominanceBufferF1!);
                        lum1 = _luminanceBufferF1!; chr1 = _chrominanceBufferF1!;
                        // Fused chroma demod+LPF field 1
                        EnsureBuffer(ref _uScratchF1, len1);
                        EnsureBuffer(ref _vScratchF1, len1);
                        EnsureBuffer(ref _uFilteredF1, len1);
                        EnsureBuffer(ref _vFilteredF1, len1);
                        var revLPF = GetReversedFilter(chromaLPF);
                        FusedChromaDemodAndFilter(chr1!, sampleRate, samplesPerLine, field1StartLine, chromaLPF, revLPF, _uScratchF1!, _vScratchF1!, _uFilteredF1!, _vFilteredF1!);
                        u1 = _uFilteredF1!; v1 = _vFilteredF1!;
                    },
                    () =>
                    {
                        // Field 2 luma/chroma
                        int len2 = field2.Length;
                        EnsureBuffer(ref _luminanceBufferF2, len2);
                        EnsureBuffer(ref _chrominanceBufferF2, len2);
                        ApplyTwoFiltersStaticToDest(field2, lumaFilter, chromaFilter, _luminanceBufferF2!, _chrominanceBufferF2!);
                        lum2 = _luminanceBufferF2!; chr2 = _chrominanceBufferF2!;
                        // Fused chroma demod+LPF field 2
                        EnsureBuffer(ref _uScratchF2, len2);
                        EnsureBuffer(ref _vScratchF2, len2);
                        EnsureBuffer(ref _uFilteredF2, len2);
                        EnsureBuffer(ref _vFilteredF2, len2);
                        var revLPF2 = GetReversedFilter(chromaLPF);
                        FusedChromaDemodAndFilter(chr2!, sampleRate, samplesPerLine, field2StartLine, chromaLPF, revLPF2, _uScratchF2!, _vScratchF2!, _uFilteredF2!, _vFilteredF2!);
                        u2 = _uFilteredF2!; v2 = _vFilteredF2!;
                    }
                );
                t1 = Stopwatch.GetTimestamp();
                // The combined time includes both luma/chroma and chroma decode for both fields due to parallel execution.
                // Attribute proportionally: keep prior profiling buckets (split roughly by previous ratio) to retain stage visibility.
                long combined = t1 - t0;
                // Estimate split using previous frame's proportion if available; fallback 60/40 (sep/decode)
                double prevSep = _ticksLumaChroma > 0 && _ticksChromaDecode > 0 ? _ticksLumaChroma : 6;
                double prevChroma = _ticksChromaDecode > 0 ? _ticksChromaDecode : 4;
                double totalPrev = prevSep + prevChroma;
                frameLumaChromaTicks = (long)(combined * (prevSep / totalPrev));
                frameChromaDecodeTicks = combined - frameLumaChromaTicks;

                // Assign refs for downstream steps
                var lum1Ref = lum1!; var lum2Ref = lum2!; var u1Ref = u1!; var v1Ref = v1!; var u2Ref = u2!; var v2Ref = v2!;


                // OPTIONAL: crop Y/U/V after decode (safe; burst already used)
                //    var activeWidth = samplesPerLine;
                t0 = Stopwatch.GetTimestamp();
                CropToActiveInto(lum1Ref, samplesPerLine, sampleRate, _lum1Active!, copyWidth);
                CropToActiveInto(u1Ref, samplesPerLine, sampleRate, _u1Active!, copyWidth);
                CropToActiveInto(v1Ref, samplesPerLine, sampleRate, _v1Active!, copyWidth);
                CropToActiveInto(lum2Ref, samplesPerLine, sampleRate, _lum2Active!, copyWidth);
                CropToActiveInto(u2Ref, samplesPerLine, sampleRate, _u2Active!, copyWidth);
                CropToActiveInto(v2Ref, samplesPerLine, sampleRate, _v2Active!, copyWidth);
                t1 = Stopwatch.GetTimestamp();
                frameCropTicks = t1 - t0;
                int activeWidth = copyWidth;

                // Convert each field (use activeWidth as samplesPerLine)
                int fieldHeight = activeLinesPerField; // lum1Active lines
                EnsureRgbBuffer(ref _rgbField1Buffer, fieldHeight, activeWidth);
                EnsureRgbBuffer(ref _rgbField2Buffer, fieldHeight, activeWidth);
                var rgbField1 = _rgbField1Buffer!;
                var rgbField2 = _rgbField2Buffer!;
                t0 = Stopwatch.GetTimestamp();
                ConvertYUVToRGB_BT601_OptimizedInto(_lum1Active!, _u1Active!, _v1Active!, activeWidth, rgbField1);
                ConvertYUVToRGB_BT601_OptimizedInto(_lum2Active!, _u2Active!, _v2Active!, activeWidth, rgbField2);
                t1 = Stopwatch.GetTimestamp();
                frameRgbTicks = t1 - t0;

                // Interleave fields for display into reusable frame buffer
                int finalHeight = Math.Min(PAL_VISIBLE_LINES, fieldHeight * 2);
                EnsureRgbBuffer(ref _rgbInterleavedFrameBuffer, finalHeight, activeWidth);
                t0 = Stopwatch.GetTimestamp();
                InterleaveFieldsInto(rgbField1, rgbField2, _rgbInterleavedFrameBuffer!);
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
                continue; // short read or error
            }

            DisplayVideoFrame(rgbFrame);

            if (!isFile && rgbFrame == null)
            {
                Console.WriteLine("Continuous stream ended.");
                break;
            }
        }
    }

    // Streaming reader + AM magnitude demodulation into provided buffer.
    // Returns false if EOF before filling buffer.
    private bool ReadAndDemodFrame(Stream source, int samplesPerFrame, double[] frameBuffer)
    {
        int filled = 0;
        double sum = 0;
        while (filled < samplesPerFrame)
        {
            int need = samplesPerFrame - filled;
            if (source is FileStream fsrc)
            {
                var seg = IQWavReader.ReadIQIntoRingOptimized(fsrc, need);
                if (seg.IsEmpty) break; // EOF
                var first = seg.First;
                for (int i = 0; i < first.Length; i++)
                {
                    double mag = first[i].Magnitude;
                    frameBuffer[filled++] = mag;
                    sum += mag;
                }
                var second = seg.Second;
                for (int i = 0; i < second.Length && filled < samplesPerFrame; i++)
                {
                    double mag = second[i].Magnitude;
                    frameBuffer[filled++] = mag;
                    sum += mag;
                }
            }
            else
            {
                // Fallback continuous stream reader: assume int8 I,Q interleaved.
                const int bytesPerSample = 2; // I + Q (signed 8-bit)
                int bytesNeeded = need * bytesPerSample;
                byte[] tmp = new byte[Math.Min(bytesNeeded, 8192)];
                int read = source.Read(tmp, 0, tmp.Length);
                if (read <= 0) break;
                int samplesRead = read / bytesPerSample;
                for (int s = 0; s < samplesRead && filled < samplesPerFrame; s++)
                {
                    double iVal = (sbyte)tmp[2 * s] / 128.0;
                    double qVal = (sbyte)tmp[2 * s + 1] / 128.0;
                    double mag = Math.Sqrt(iVal * iVal + qVal * qVal);
                    frameBuffer[filled++] = mag;
                    sum += mag;
                }
            }
        }
        if (filled < samplesPerFrame) return false;
        double dc = sum / samplesPerFrame;
        for (int i = 0; i < samplesPerFrame; i++) frameBuffer[i] -= dc;
        return true;
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

    // In-place variant: copies cropped active region into provided destination (dest length must be >= lines * copyWidth)
    private void CropToActiveInto(double[] signal, int samplesPerLine, int sampleRate, double[] dest, int copyWidth)
    {
        int lines = signal.Length / samplesPerLine;
        // Original logic calculates desiredActiveCol but then uses: copyWidth = Max(activeWidth, Max(0, samplesPerLine - desiredActiveCol))
        // We receive copyWidth precomputed; perform straight per-line copy from column 0 for now (mirrors earlier implemented behavior).
        for (int ln = 0; ln < lines; ln++)
        {
            int src = ln * samplesPerLine;
            int dst = ln * copyWidth;
            // Safe length guard
            if (src + copyWidth > signal.Length) break;
            Buffer.BlockCopy(signal, src * sizeof(double), dest, dst * sizeof(double), copyWidth * sizeof(double));
        }
    }

    // auto-find the horizontal start by detecting the H-sync pulse per line and computing the active-video start from timing 
    // (4.7 µs sync + 5.8 µs back porch). 
    private static int EstimateHorizontalOffset(double[] videoSignal, int frameStart, int samplesPerLine, int sampleRate, int linesToUse = 24)
    {
        const double HSYNC_US = 4.7e-6;
        const double BACK_PORCH_US = 5.8e-6;
        const double SEARCH_US = 30e-6; // include front porch, sync, back porch
        const int LPF_TAPS = 51;        // must match LightLowPass()
        int gd = (LPF_TAPS - 1) / 2;    // FIR group delay in samples

        int searchCols = Math.Min(samplesPerLine, (int)Math.Round(SEARCH_US * sampleRate));
        int hsyncSamples = Math.Max(8, (int)Math.Round(HSYNC_US * sampleRate));
        int desiredActiveCol = (int)Math.Round((HSYNC_US + BACK_PORCH_US) * sampleRate); // ≈ 10.5 µs

        List<int> activeStarts = new();

        for (int ln = 0; ln < linesToUse; ln++)
        {
            int lineStart = frameStart + ln * samplesPerLine;
            if (lineStart + searchCols > videoSignal.Length) break;

            // Window and light LPF to kill chroma ripple
            double[] win = new double[searchCols];
            Array.Copy(videoSignal, lineStart, win, 0, searchCols);
            win = LightLowPass(win);

            // Auto polarity per line (ensure sync is negative)
            (double wMin, double wMax) = MinMax(win);
            if (Math.Abs(wMax) > Math.Abs(wMin))
            {
                for (int i = 0; i < win.Length; i++) win[i] = -win[i];
                (wMin, wMax) = (-wMax, -wMin);
            }

            // Rectangular correlation with a 4.7 µs window → find deepest segment
            int bestStart = -1;
            double sum = 0, bestSum = double.PositiveInfinity;
            int L = hsyncSamples;
            int N = win.Length;

            if (L <= N)
            {
                for (int i = 0; i < L; i++) sum += win[i];
                bestSum = sum; bestStart = 0;
                for (int i = L; i < N; i++)
                {
                    sum += win[i] - win[i - L];
                    if (sum < bestSum)
                    {
                        bestSum = sum;
                        bestStart = i - L + 1;
                    }
                }
            }

            bool added = false;
            if (bestStart >= 0)
            {
                // Refine with robust threshold
                double med = Median(win);
                double mad = Median(win.Select(v => Math.Abs(v - med)).ToArray());
                if (mad <= 1e-12) mad = (wMax - wMin) / 50.0;
                double thr = med - 2.5 * mad;

                int rs = bestStart;
                while (rs > 0 && win[rs - 1] < thr) rs--;
                int re = Math.Min(N - 1, bestStart + L - 1);
                while (re + 1 < N && win[re + 1] < thr) re++;

                int len = re - rs + 1;
                if (len >= (int)(0.5 * hsyncSamples) && len <= (int)(2.0 * hsyncSamples))
                {
                    // Compensate LPF group delay
                    int syncStartCol = Math.Max(0, rs - gd);
                    int activeCol = syncStartCol + (int)Math.Round((HSYNC_US + BACK_PORCH_US) * sampleRate);
                    activeStarts.Add(activeCol);
                    added = true;
                }
            }

            // Fallback: largest negative edge (sync leading edge)
            if (!added)
            {
                int edgeIdx = -1;
                double minDer = double.PositiveInfinity;
                for (int i = 1; i < N; i++)
                {
                    double d = win[i] - win[i - 1];
                    if (d < minDer) { minDer = d; edgeIdx = i; }
                }
                if (edgeIdx > 0)
                {
                    int syncEdge = Math.Max(0, edgeIdx - gd); // compensate LPF delay
                    int activeCol = syncEdge + (int)Math.Round((HSYNC_US + BACK_PORCH_US) * sampleRate);
                    activeStarts.Add(activeCol);
                }
            }
        }

        if (activeStarts.Count == 0)
        {
            Console.WriteLine("HSYNC not found; ");
            return 0;
        }

        int medianActiveCol = activeStarts.OrderBy(x => x).ElementAt(activeStarts.Count / 2);
        int offset = desiredActiveCol - medianActiveCol; // positive ⇒ shift right; negative ⇒ shift left
        offset = Math.Max(-samplesPerLine / 2, Math.Min(samplesPerLine / 2, offset));

        Console.WriteLine($"Auto horizontal offset = {offset} (median={medianActiveCol}, desired={desiredActiveCol}, gd={gd}, lines={activeStarts.Count})");
        return offset;
    }

    private static (double min, double max) MinMax(ReadOnlySpan<double> s)
    {
        if (s.IsEmpty) return (double.NaN, double.NaN);
        double min = s[0], max = s[0];
        for (int i = 1; i < s.Length; i++)
        {
            double v = s[i];
            if (v < min) min = v;
            if (v > max) max = v;
        }
        return (min, max);
    }

    // Works even if sync polarity is flipped by AM magnitude demod.
    // Threshold adapts to your signal distribution (median/MAD), not a fixed “min*0.5”.
    // Broad-pulse run-length avoids missing due to single-sample noise.
    // Low-pass reduces 4.43 MHz ripple so runs are cleaner.
    private int FindFrameStart(double[] videoSignal, int sampleRate, int samplesPerLine)
    {
        Console.WriteLine("Searching for frame start (robust VBI detection)...");
        if (videoSignal == null || videoSignal.Length < samplesPerLine * 50)
        {
            Console.WriteLine("Signal too short for VBI search, fallback to 0.");
            return 0;
        }

        int searchLength = Math.Min(videoSignal.Length, Math.Max(sampleRate / 4, samplesPerLine * PAL_LINES_PER_FRAME));
        var segment = LightLowPass(videoSignal.Take(searchLength).ToArray());

        // auto polarity
        (double segMin, double segMax) = MinMax(segment);
        bool inverted = Math.Abs(segMax) > Math.Abs(segMin);
        if (inverted)
        {
            for (int j = 0; j < segment.Length; j++) segment[j] = -segment[j];
            (segMin, segMax) = MinMax(segment);
        }

        // robust threshold
        double median = Median(segment);
        double mad = Median(segment.Select(x => Math.Abs(x - median)).ToArray());
        if (mad <= 1e-12) mad = (segMax - segMin) / 50.0;
        double syncThreshold = median - 3.0 * mad;

        int broadMin = (int)Math.Round(20e-6 * sampleRate); // broad pulse ≥20 µs
        int i = 0;

        while (i < segment.Length - samplesPerLine)
        {
            while (i < segment.Length && !(segment[i] < syncThreshold)) i++;
            if (i >= segment.Length) break;

            int runStart = i;
            while (i < segment.Length && segment[i] < syncThreshold) i++;
            int runLen = i - runStart;

            if (runLen >= broadMin)
            {
                // confirm a second broad pulse nearby
                int j = i;
                int endSearch = Math.Min(segment.Length, runStart + 4 * samplesPerLine);
                bool secondBroadFound = false;
                while (j < endSearch)
                {
                    while (j < endSearch && !(segment[j] < syncThreshold)) j++;
                    int r2s = j;
                    while (j < endSearch && segment[j] < syncThreshold) j++;
                    int r2len = j - r2s;
                    if (r2len >= broadMin) { secondBroadFound = true; break; }
                }

                if (secondBroadFound)
                {
                    // SPEC default (kept as fallback): int approx = runStart + (int)Math.Round(22.5 * samplesPerLine);
                    // Measurement-based refinement: find first line with valid burst after V-sync
                    int refined = RefineToFirstActiveLine(videoSignal, runStart, sampleRate, samplesPerLine);
                    Console.WriteLine($"VBI detected (inv={inverted}). Field1 active ≈ {refined}");
                    return refined;
                }
            }
            // continue scanning
        }

        Console.WriteLine("Could not find VBI pattern, using default of 0.");
        return 0;
    }

    // Find the first line after V-sync that has a strong 4.4336 MHz burst, and return its active start index.
    private int RefineToFirstActiveLine(double[] video, int vsyncFirstBroadStart, int sampleRate, int samplesPerLine)
    {
        // timing (PAL)
        int hsyncUs = (int)Math.Round(4.7e-6 * sampleRate);
        int breezewayUs = (int)Math.Round(0.6e-6 * sampleRate);
        int burstLen = (int)Math.Round(2.25e-6 * sampleRate);
        int backPorchToBurst = hsyncUs + breezewayUs;
        int desiredActiveCol = (int)Math.Round((4.7e-6 + 5.8e-6) * sampleRate);

        // Examine a few dozen lines after V-sync and pick the first with strong burst
        int firstLineIdx = Math.Max(0, (vsyncFirstBroadStart / samplesPerLine) - 1); // start a tad early
        int lastLineIdx = Math.Min(firstLineIdx + 60, (video.Length / samplesPerLine) - 1);

        // Establish a robust threshold by measuring burst energy over the candidate set
        List<(int line, double power)> candidates = new();
        for (int ln = firstLineIdx; ln <= lastLineIdx; ln++)
        {
            int lineStart = ln * samplesPerLine;
            if (lineStart + backPorchToBurst + burstLen >= video.Length) break;

            // Gate on expected burst window relative to the line start
            int burstStart = lineStart + backPorchToBurst;
            double p = GoertzelPower(video, burstStart, burstLen, PAL_COLOR_CARRIER_FREQ, sampleRate);
            candidates.Add((ln, p));
        }
        if (candidates.Count == 0) return vsyncFirstBroadStart + (int)Math.Round(22.5 * samplesPerLine);

        // Robust threshold: median + 3*MAD of burst power
        double[] powers = candidates.Select(t => t.power).ToArray();
        double m = Median(powers);
        double md = Median(powers.Select(v => Math.Abs(v - m)).ToArray());
        if (md <= 1e-20) md = (powers.Max() - powers.Min()) / 50.0;
        double thr = m + 3.0 * md;

        // First line with power above threshold (and preferably 3 consecutive hits)
        for (int k = 0; k < candidates.Count; k++)
        {
            bool strong = candidates[k].power > thr;
            bool nextStrong = (k + 1 < candidates.Count) && (candidates[k + 1].power > thr);
            bool next2Strong = (k + 2 < candidates.Count) && (candidates[k + 2].power > thr);
            if (strong && (nextStrong || next2Strong))
            {
                int ln = candidates[k].line;
                return ln * samplesPerLine + desiredActiveCol; // return active-video start
            }
        }

        // Fallback to the strongest burst line
        var best = candidates.OrderByDescending(t => t.power).First();
        return best.line * samplesPerLine + desiredActiveCol;
    }

    // Goertzel detector for burst power at subcarrier frequency in a short window
    private static double GoertzelPower(double[] x, int offset, int length, double targetFreq, int sampleRate)
    {
        double w = 2 * Math.PI * targetFreq / sampleRate;
        double cosw = Math.Cos(w);
        double coeff = 2 * cosw;
        double s0 = 0, s1 = 0, s2 = 0;

        int end = Math.Min(offset + length, x.Length);
        for (int n = offset; n < end; n++)
        {
            s0 = x[n] + coeff * s1 - s2;
            s2 = s1;
            s1 = s0;
        }
        double real = s1 - s2 * cosw;
        double imag = s2 * Math.Sin(w);
        return real * real + imag * imag;
    }

    private static double Median(double[] a)
    {
        var b = (double[])a.Clone();
        Array.Sort(b);
        int n = b.Length;
        return (n % 2 == 1) ? b[n / 2] : 0.5 * (b[n / 2 - 1] + b[n / 2]);
    }

    private static double[] LightLowPass(double[] x)
    {
        // 51-tap Hamming moving-average-like FIR (simple, fast)
        const int N = 51;
        if (x.Length < N) return x;
        double[] w = new double[N];
        for (int i = 0; i < N; i++) w[i] = 0.54 - 0.46 * Math.Cos(2 * Math.PI * i / (N - 1));
        double sumw = w.Sum();
        for (int i = 0; i < N; i++) w[i] /= sumw;

        double[] y = new double[x.Length];
        for (int n = 0; n < x.Length; n++)
        {
            double acc = 0;
            for (int k = 0; k < N; k++)
            {
                int idx = n - k;
                if (idx >= 0) acc += x[idx] * w[k];
            }
            y[n] = acc;
        }
        return y;
    }

    private (double[] luminance, double[] chrominance) SeparateLumaChromaPooled(double[] videoSignal, int sampleRate)
    {
        int len = videoSignal.Length;
        EnsureBuffer(ref _luminanceBuffer, len);
        EnsureBuffer(ref _chrominanceBuffer, len);
        // We previously performed two independent FIR passes (LPF for luma, BPF for chroma).
        // To reduce memory bandwidth and halve passes through the large video buffer we
        // compute both convolutions in a single streaming/vectorized loop.
        double lumaCutoff = Math.Min(_profile.LumaCutoffHz, 0.45 * sampleRate);
        var lumaFilter = CreateLowPassFilter(lumaCutoff, sampleRate, LUMA_LPF_TAPS);
        double chromaLow = _profile.ChromaLowHz;
        double chromaHigh = _profile.ChromaHighHz;
        var chromaFilter = CreateBandPassFilter(chromaLow, chromaHigh, sampleRate, CHROMA_SEPARATION_TAPS);
        ApplyTwoFiltersStaticToDest(videoSignal, lumaFilter, chromaFilter, _luminanceBuffer!, _chrominanceBuffer!);
        return (_luminanceBuffer!, _chrominanceBuffer!);
    }

    private (double[] uComponent, double[] vComponent) DecodeChroma(double[] chrominance, int sampleRate, int samplesPerLine, int startLineOffset)
    {
        int len = chrominance.Length;
        EnsureBuffer(ref _uScratch, len);   // demod intermediate (can later be removed if we fuse demod+FIR)
        EnsureBuffer(ref _vScratch, len);
        EnsureBuffer(ref _uFiltered, len);
        EnsureBuffer(ref _vFiltered, len);
        var pll = new ColorPll(sampleRate);
        for (int i = 0; i < len; i++)
        {
            int lineInBuffer = i / samplesPerLine;
            int sampleInLine = i % samplesPerLine;
            int absoluteLine = startLineOffset + lineInBuffer; // cross-field absolute line number
            bool isVInverted = (absoluteLine % 2 != 0);
            var referenceCarrier = pll.GetReference(new Complex(chrominance[i], 0), sampleInLine, isVInverted);
            var demod = new Complex(chrominance[i], 0) * Complex.Conjugate(referenceCarrier);
            _uScratch![i] = demod.Real;
            _vScratch![i] = isVInverted ? -demod.Imaginary : demod.Imaginary;
        }
        // Apply the SAME low-pass filter to U & V simultaneously (shared memory reads of filter and partial SIMD dot products)
        double chromaCutoff = 1.3e6; // PAL chroma baseband width
        var chromaLPF = CreateLowPassFilter(chromaCutoff, sampleRate, CHROMA_BASEBAND_LPF_TAPS);
        ApplySameFilterTwoSignalsStatic(_uScratch!, _vScratch!, chromaLPF, _uFiltered!, _vFiltered!);
        return (_uFiltered!, _vFiltered!);
    }

    // In-place variant: writes into provided destination buffer (height = y.Length / samplesPerLine, width = provided dest width)
    private void ConvertYUVToRGB_BT601_OptimizedInto(double[] y, double[] u, double[] v, int samplesPerLine, byte[,,] dest)
    {
        int destHeight = dest.GetLength(0);
        int destWidth = dest.GetLength(1);
        int height = Math.Min(destHeight, y.Length / samplesPerLine);
        int width = Math.Min(destWidth, samplesPerLine);

        const double c_rv = 1.402;
        const double c_gu = -0.344;
        const double c_gv = -0.714;
        const double c_bu = 1.772;

        for (int row = 0; row < height; row++)
        {
            int baseIdx = row * samplesPerLine;
            for (int col = 0; col < width; col++)
            {
                int index = baseIdx + col;
                double yVal = Clamp(y[index], 0, 1);
                double uVal = Clamp(u[index], -0.5, 0.5);
                double vVal = Clamp(v[index], -0.5, 0.5);

                double r = yVal + c_rv * vVal;
                double g = yVal + c_gu * uVal + c_gv * vVal;
                double b = yVal + c_bu * uVal;

                dest[row, col, 0] = (byte)(Clamp(r, 0, 1) * 255);
                dest[row, col, 1] = (byte)(Clamp(g, 0, 1) * 255);
                dest[row, col, 2] = (byte)(Clamp(b, 0, 1) * 255);
            }
        }
    }

    // In-place interleave into destination frame (avoids allocating a new 3D frame array)
    private void InterleaveFieldsInto(byte[,,] field1, byte[,,] field2, byte[,,] dest)
    {
        int fieldLines = field1.GetLength(0);
        int width = Math.Min(field1.GetLength(1), field2.GetLength(1));
        int destHeight = dest.GetLength(0);
        int height = Math.Min(destHeight, fieldLines * 2);
        int outLine = 0;
        for (int i = 0; i < fieldLines && outLine < height - 1; i++)
        {
            for (int j = 0; j < width; j++)
            {
                if (_fieldOrder == FieldOrder.BottomFieldFirst)
                {
                    dest[outLine, j, 0] = field2[i, j, 0];
                    dest[outLine, j, 1] = field2[i, j, 1];
                    dest[outLine, j, 2] = field2[i, j, 2];

                    dest[outLine + 1, j, 0] = field1[i, j, 0];
                    dest[outLine + 1, j, 1] = field1[i, j, 1];
                    dest[outLine + 1, j, 2] = field1[i, j, 2];
                }
                else
                {
                    dest[outLine, j, 0] = field1[i, j, 0];
                    dest[outLine, j, 1] = field1[i, j, 1];
                    dest[outLine, j, 2] = field1[i, j, 2];

                    dest[outLine + 1, j, 0] = field2[i, j, 0];
                    dest[outLine + 1, j, 1] = field2[i, j, 1];
                    dest[outLine + 1, j, 2] = field2[i, j, 2];
                }
            }
            outLine += 2;
        }
    }

    // Thread-safe shared filter cache (keyed by cutoff/sampleRate/tap count)
    private readonly ConcurrentDictionary<string, double[]> _filterCache = new();

    private double[] CreateLowPassFilter(double cutoffFreq, int sampleRate, int filterLength = 101)
    {
        string filterKey = $"LPF_{cutoffFreq}_{sampleRate}_{filterLength}";
        return _filterCache.GetOrAdd(filterKey, static key =>
        {
            // key format: LPF_{cutoff}_{sampleRate}_{length}
            var parts = key.Split('_');
            double cutoff = double.Parse(parts[1], System.Globalization.CultureInfo.InvariantCulture);
            int sr = int.Parse(parts[2], System.Globalization.CultureInfo.InvariantCulture);
            int fl = int.Parse(parts[3], System.Globalization.CultureInfo.InvariantCulture);
            double[] filter = new double[fl];
            double fc = cutoff / sr;
            for (int i = 0; i < fl; i++)
            {
                int n = i - fl / 2;
                double v = (n == 0) ? 2 * fc : Math.Sin(2 * Math.PI * fc * n) / (Math.PI * n);
                v *= 0.54 - 0.46 * Math.Cos(2 * Math.PI * i / (fl - 1)); // Hamming window
                filter[i] = v;
            }
            return filter;
        });
    }

    private double[] CreateBandPassFilter(double lowFreq, double highFreq, int sampleRate, int filterLength = 101)
    {
        string key = $"BPF_{lowFreq}_{highFreq}_{sampleRate}_{filterLength}";
        return _filterCache.GetOrAdd(key, _ =>
        {
            var lpf1 = CreateLowPassFilter(highFreq, sampleRate, filterLength);
            var lpf2 = CreateLowPassFilter(lowFreq, sampleRate, filterLength);
            double[] bandPass = new double[lpf1.Length];
            for (int i = 0; i < bandPass.Length; i++) bandPass[i] = lpf1[i] - lpf2[i];
            return bandPass;
        });
    }

    // Fused chroma demod + low-pass filtering (single pass per field)
    private void FusedChromaDemodAndFilter(double[] chroma, int sampleRate, int samplesPerLine, int startLineOffset,
        double[] lpf, double[] lpfRev, double[] uScratch, double[] vScratch, double[] uOut, double[] vOut)
    {
        int len = chroma.Length;
        int taps = lpf.Length;
        int warm = taps - 1;
        int maxWarm = Math.Min(warm, len - 1);
        var pll = new ColorPll(sampleRate);
        int simdWidth = Vector<double>.Count;
        bool useSimd = Vector.IsHardwareAccelerated && taps >= simdWidth * 2;
        for (int i = 0; i < len; i++)
        {
            int lineInBuffer = i / samplesPerLine;
            int sampleInLine = i % samplesPerLine;
            int absLine = startLineOffset + lineInBuffer;
            bool invertV = (absLine % 2 != 0);
            double c = chroma[i];
            var refCarrier = pll.GetReference(new Complex(c, 0), sampleInLine, invertV);
            var demod = new Complex(c, 0) * Complex.Conjugate(refCarrier);
            double u = demod.Real; double v = invertV ? -demod.Imaginary : demod.Imaginary;
            uScratch[i] = u; vScratch[i] = v;
            double sumU = 0, sumV = 0;
            if (i <= maxWarm)
            {
                int maxTap = Math.Min(taps - 1, i);
                for (int k = 0; k <= maxTap; k++) { double fk = lpf[k]; sumU += uScratch[i - k] * fk; sumV += vScratch[i - k] * fk; }
            }
            else
            {
                int start = i - taps + 1;
                if (useSimd)
                {
                    int k = 0; int limit = taps - simdWidth; var vaccU = Vector<double>.Zero; var vaccV = Vector<double>.Zero;
                    while (k <= limit)
                    {
                        var vUS = new Vector<double>(uScratch, start + k);
                        var vVS = new Vector<double>(vScratch, start + k);
                        var vF = new Vector<double>(lpfRev, k);
                        vaccU += vUS * vF; vaccV += vVS * vF; k += simdWidth;
                    }
                    for (int s = 0; s < simdWidth; s++) { sumU += vaccU[s]; sumV += vaccV[s]; }
                    for (; k < taps; k++) { double fk = lpfRev[k]; sumU += uScratch[start + k] * fk; sumV += vScratch[start + k] * fk; }
                }
                else
                {
                    for (int k = 0; k < taps; k++) { double fk = lpfRev[k]; sumU += uScratch[start + k] * fk; sumV += vScratch[start + k] * fk; }
                }
            }
            uOut[i] = sumU; vOut[i] = sumV;
        }
    }

    private static void ApplyFilterStaticToDest(double[] signal, double[] filter, double[] dest)
    {
        // Vectorized convolution (causal) with reversed filter trick for steady-state region.
        int n = signal.Length;
        int taps = filter.Length;
        if (n == 0) return;

        // Warm-up (edges) where full taps not yet available – scalar
        int warm = taps - 1;
        int maxWarm = Math.Min(warm, n - 1);
        for (int i = 0; i <= maxWarm; i++)
        {
            double sum = 0;
            int maxTap = Math.Min(taps - 1, i);
            for (int j = 0; j <= maxTap; j++) sum += signal[i - j] * filter[j];
            dest[i] = sum;
        }

        if (maxWarm == n - 1) return; // all done (short signal)

        // Prepare reversed filter for steady region so we can do forward dot product over contiguous slice
        double[] filtRev = GetReversedFilter(filter);
        int simdWidth = Vector<double>.Count;
        bool useSimd = Vector.IsHardwareAccelerated && taps >= simdWidth * 2; // require at least two vectors

        for (int i = warm; i < n; i++)
        {
            // Full window available: samples [i - taps + 1 .. i]
            int start = i - taps + 1;
            double sum;
            if (useSimd)
            {
                sum = 0;
                int k = 0;
                int limit = taps - simdWidth;
                Vector<double> vacc = Vector<double>.Zero;
                // Slice is contiguous forward segment
                // We avoid allocating by manual load
                while (k <= limit)
                {
                    var vSig = new Vector<double>(signal, start + k);
                    var vFlt = new Vector<double>(filtRev, k);
                    vacc += vSig * vFlt;
                    k += simdWidth;
                }
                sum = 0;
                for (int s = 0; s < simdWidth; s++) sum += vacc[s];
                for (; k < taps; k++) sum += signal[start + k] * filtRev[k];
            }
            else
            {
                sum = 0;
                for (int k = 0; k < taps; k++) sum += signal[start + k] * filtRev[k];
            }
            dest[i] = sum;
        }
    }

    // NEW: Compute two independent FIR outputs (signal * filterA, signal * filterB) in a single pass.
    // This halves memory traffic vs calling ApplyFilter twice when both outputs are required.
    private static void ApplyTwoFiltersStaticToDest(double[] signal, double[] filterA, double[] filterB, double[] destA, double[] destB)
    {
        int n = signal.Length;
        int tapsA = filterA.Length;
        int tapsB = filterB.Length;
        if (n == 0) return;
        // For simplicity require equal length; if not equal we fallback to independent calls.
        if (tapsA != tapsB)
        {
            ApplyFilterStaticToDest(signal, filterA, destA);
            ApplyFilterStaticToDest(signal, filterB, destB);
            return;
        }
        int taps = tapsA;
        int warm = taps - 1;
        int maxWarm = Math.Min(warm, n - 1);
        double[] revA = GetReversedFilter(filterA);
        double[] revB = ReferenceEquals(filterA, filterB) ? revA : GetReversedFilter(filterB);
        int simdWidth = Vector<double>.Count;
        bool useSimd = Vector.IsHardwareAccelerated && taps >= simdWidth * 2;

        // Warm-up edge (scalar)
        for (int i = 0; i <= maxWarm; i++)
        {
            double sumA = 0, sumB = 0; int maxTap = Math.Min(taps - 1, i);
            for (int k = 0; k <= maxTap; k++)
            {
                double s = signal[i - k];
                sumA += s * filterA[k];
                sumB += s * filterB[k];
            }
            destA[i] = sumA; destB[i] = sumB;
        }
        if (maxWarm == n - 1) return;
        for (int i = warm; i < n; i++)
        {
            int start = i - taps + 1;
            double sumA, sumB;
            if (useSimd)
            {
                int k = 0; int limit = taps - simdWidth; var vaccA = Vector<double>.Zero; var vaccB = Vector<double>.Zero;
                while (k <= limit)
                {
                    var vSig = new Vector<double>(signal, start + k);
                    var vA = new Vector<double>(revA, k);
                    var vB = new Vector<double>(revB, k);
                    vaccA += vSig * vA;
                    vaccB += vSig * vB;
                    k += simdWidth;
                }
                sumA = 0; sumB = 0;
                for (int s = 0; s < simdWidth; s++) { sumA += vaccA[s]; sumB += vaccB[s]; }
                for (; k < taps; k++)
                {
                    double sig = signal[start + k];
                    sumA += sig * revA[k];
                    sumB += sig * revB[k];
                }
            }
            else
            {
                sumA = 0; sumB = 0;
                for (int k = 0; k < taps; k++)
                {
                    double sig = signal[start + k];
                    sumA += sig * revA[k];
                    sumB += sig * revB[k];
                }
            }
            destA[i] = sumA; destB[i] = sumB;
        }
    }

    // NEW: Apply the same filter to two different source signals simultaneously.
    private static void ApplySameFilterTwoSignalsStatic(double[] signal1, double[] signal2, double[] filter, double[] dest1, double[] dest2)
    {
        int n = signal1.Length;
        int taps = filter.Length;
        if (n == 0) return;
        int warm = taps - 1;
        int maxWarm = Math.Min(warm, n - 1);
        double[] rev = GetReversedFilter(filter);
        int simdWidth = Vector<double>.Count;
        bool useSimd = Vector.IsHardwareAccelerated && taps >= simdWidth * 2;
        for (int i = 0; i <= maxWarm; i++)
        {
            double sum1 = 0, sum2 = 0; int maxTap = Math.Min(taps - 1, i);
            for (int k = 0; k <= maxTap; k++)
            {
                double s1 = signal1[i - k]; double s2 = signal2[i - k]; double fk = filter[k];
                sum1 += s1 * fk; sum2 += s2 * fk;
            }
            dest1[i] = sum1; dest2[i] = sum2;
        }
        if (maxWarm == n - 1) return;
        for (int i = warm; i < n; i++)
        {
            int start = i - taps + 1; double sum1, sum2;
            if (useSimd)
            {
                int k = 0; int limit = taps - simdWidth; var vacc1 = Vector<double>.Zero; var vacc2 = Vector<double>.Zero;
                while (k <= limit)
                {
                    var vS1 = new Vector<double>(signal1, start + k);
                    var vS2 = new Vector<double>(signal2, start + k);
                    var vF = new Vector<double>(rev, k);
                    vacc1 += vS1 * vF; vacc2 += vS2 * vF; k += simdWidth;
                }
                sum1 = 0; sum2 = 0; for (int s = 0; s < simdWidth; s++) { sum1 += vacc1[s]; sum2 += vacc2[s]; }
                for (; k < taps; k++) { double fk = rev[k]; sum1 += signal1[start + k] * fk; sum2 += signal2[start + k] * fk; }
            }
            else
            {
                sum1 = 0; sum2 = 0;
                for (int k = 0; k < taps; k++) { double fk = rev[k]; sum1 += signal1[start + k] * fk; sum2 += signal2[start + k] * fk; }
            }
            dest1[i] = sum1; dest2[i] = sum2;
        }
    }

    private static double Clamp(double value, double min, double max)
    {
        return Math.Max(min, Math.Min(max, value));
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

// This class implements a Phase-Locked Loop (PLL) to recover a stable color subcarrier
// that is phase-aligned with the incoming signal's color burst.
public class ColorPll
{
    private readonly double _sampleRate;
    private readonly double _burstStartTime; // Time after H-sync where burst starts
    private readonly double _burstDuration;  // Duration of the color burst

    // PLL state variables
    private double _frequency = PALDecoder.PAL_COLOR_CARRIER_FREQ;
    private double _phaseErrorIntegrator = 0.0;
    private double _phaseIncrement;
    private Complex _osc = Complex.One; // running oscillator (unit magnitude)
    private Complex _rot;               // per-sample rotation factor exp(j*phaseIncrement)
    private int _renormCounter;

    // PLL loop filter gains (these may need tuning for noisy signals)
    private const double ProportionalGain = 0.1;
    private const double IntegralGain = 0.005;

    private static readonly Complex PhasePos = Complex.FromPolarCoordinates(1, -Math.PI / 4);   // -45°
    private static readonly Complex PhaseNeg = Complex.FromPolarCoordinates(1, -3 * Math.PI / 4); // -135°

    public ColorPll(int sampleRate)
    {
        _sampleRate = sampleRate;
        // PAL spec: H-sync (4.7µs) + Breezeway (0.6µs) = 5.3µs
        _burstStartTime = 5.3e-6;
        _burstDuration = 2.25e-6;
        UpdatePhaseIncrement();
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void UpdatePhaseIncrement()
    {
        _phaseIncrement = 2 * Math.PI * _frequency / _sampleRate;
        _rot = Complex.FromPolarCoordinates(1, _phaseIncrement);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void StepOscillator()
    {
        _osc *= _rot;
        // Periodic renormalization to avoid floating drift
        if ((++_renormCounter & 1023) == 0)
        {
            double mag = _osc.Magnitude;
            _osc /= mag;
        }
    }

    // Generates one sample of the reference carrier and updates the PLL state
    public Complex GetReference(Complex chromaSample, int sampleIndexInLine, bool isVInverted)
    {
        double timeInLine = sampleIndexInLine / _sampleRate;

        // --- Phase Detection (only during the color burst) ---
        if (timeInLine >= _burstStartTime && timeInLine < _burstStartTime + _burstDuration)
        {
            // Reference with expected V-axis phase swing
            var burstReference = _osc * (isVInverted ? PhaseNeg : PhasePos);
            double phaseError = (chromaSample * Complex.Conjugate(burstReference)).Phase;

            // --- Loop Filter ---
            // Update the integrator (I-term)
            _phaseErrorIntegrator += phaseError * IntegralGain;
            // Update frequency
            _frequency = PALDecoder.PAL_COLOR_CARRIER_FREQ + phaseError * ProportionalGain + _phaseErrorIntegrator;
            UpdatePhaseIncrement();
        }
        StepOscillator();
        return _osc;
    }
}
