using ScottPlot;
using System.Numerics;
using Microsoft.UI.Dispatching;

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

    public PALDecoder(Plot plot, DispatcherQueue dispatcherQueue, TvSystem system = TvSystem.PAL_DK, FieldOrder fieldOrder = FieldOrder.BottomFieldFirst)
    {
        _plot = plot;
        _dispatcherQueue = dispatcherQueue;
        _profile = SystemProfile.For(system);
        _fieldOrder = fieldOrder;
    }

    public void DecodePALSignal(Complex[] iqSamples, int sampleRate)
    {
        Task.Factory.StartNew(() =>
           {
               int samplesPerLine = (int)(PAL_LINE_DURATION * sampleRate); // 640
               int samplesPerFrame = samplesPerLine * PAL_LINES_PER_FRAME; //400000

               double[] videoSignal = AMDemodulation(iqSamples);
               int frameStart = FindFrameStart(videoSignal, sampleRate, samplesPerLine);

               if (frameStart + samplesPerFrame > videoSignal.Length)
               {
                   Console.WriteLine("Not enough samples for a complete PAL frame after sync");
                   samplesPerFrame = videoSignal.Length - frameStart;
               }

               double[] frameData = new double[samplesPerFrame];
               int autoHOffset = EstimateHorizontalOffset(videoSignal, frameStart, samplesPerLine, sampleRate);
               var skipUntil = frameStart + autoHOffset;
               int numberOfFrames = (videoSignal.Length - skipUntil) / samplesPerFrame;

               var nonVideoData = (int)Math.Round((1.5 + 4.7 + 5.8) * 1e-6 * sampleRate);
               var delta = (samplesPerLine - nonVideoData) / 2 + nonVideoData;


               for (int i = 0; i < numberOfFrames; i++)
               {
                   Array.Copy(videoSignal, skipUntil + i * samplesPerFrame + delta, frameData, 0, samplesPerFrame);

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

                   double[] field1 = new double[fieldLines * samplesPerLine];
                   double[] field2 = new double[fieldLines * samplesPerLine];

                   for (int j = 0; j < fieldLinesVis; j++)
                   {
                       // Field 1: contiguous lines
                       Array.Copy(
                           frameData, (field1StartLine + j) * samplesPerLine,
                           field1, j * samplesPerLine,
                           samplesPerLine
                       );

                       // Field 2: next contiguous block (half-frame later)
                       Array.Copy(
                           frameData, (field2StartLine + j) * samplesPerLine,
                           field2, j * samplesPerLine,
                           samplesPerLine
                       );
                   }
                   // Process each field
                   var (lum1, chr1) = SeparateLumaChroma(field1, sampleRate, samplesPerLine);
                   var (lum2, chr2) = SeparateLumaChroma(field2, sampleRate, samplesPerLine);

                   // Use absolute line offsets for V-axis alternation (preserves 8-field parity across fields)
                   var (u1, v1) = DecodeChroma(chr1, sampleRate, samplesPerLine, startLineOffset: field1StartLine);
                   var (u2, v2) = DecodeChroma(chr2, sampleRate, samplesPerLine, startLineOffset: field2StartLine);


                   // OPTIONAL: crop Y/U/V after decode (safe; burst already used)
                   //    var activeWidth = samplesPerLine;
                   (lum1, var activeWidth) = CropToActive(lum1, samplesPerLine, sampleRate);
                   (u1, _) = CropToActive(u1, samplesPerLine, sampleRate);
                   (v1, _) = CropToActive(v1, samplesPerLine, sampleRate);

                   (lum2, _) = CropToActive(lum2, samplesPerLine, sampleRate);
                   (u2, _) = CropToActive(u2, samplesPerLine, sampleRate);
                   (v2, _) = CropToActive(v2, samplesPerLine, sampleRate);

                   // Convert each field (use activeWidth as samplesPerLine)
                   byte[,,] rgbField1 = ConvertYUVToRGB_BT601_Optimized(lum1, u1, v1, activeWidth);
                   byte[,,] rgbField2 = ConvertYUVToRGB_BT601_Optimized(lum2, u2, v2, activeWidth);

                   // Interleave fields for display
                   byte[,,] rgbFrame = InterleaveFields(rgbField1, rgbField2);
                   DisplayVideoFrame(rgbFrame);
               }
           }, TaskCreationOptions.LongRunning);
    }

    private (double[] data, int width) CropToActive(double[] signal, int samplesPerLine, int sampleRate)
    {
        int lines = signal.Length / samplesPerLine;
        int activeWidth = (int)Math.Round(52e-6 * sampleRate);
        int desiredActiveCol = (int)Math.Round((4.7 + 5.8) * 1e-6 * sampleRate);

        // Bounds safety
        int copyWidth = Math.Max(activeWidth, Math.Max(0, samplesPerLine - desiredActiveCol));
        double[] outSig = new double[lines * copyWidth];

        for (int ln = 0; ln < lines; ln++)
        {
            int src = ln * samplesPerLine;
            int dst = ln * copyWidth;
            Array.Copy(signal, src, outSig, dst, copyWidth);
        }
        return (outSig, copyWidth);
    }

    // auto-find the horizontal start by detecting the H-sync pulse per line and computing the active-video start from timing 
    // (4.7 µs sync + 5.8 µs back porch). 
    private int EstimateHorizontalOffset(double[] videoSignal, int frameStart, int samplesPerLine, int sampleRate, int linesToUse = 24)
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
            double wMin = win.Min(), wMax = win.Max();
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

    // Interleave two fields into a single frame (PAL interlacing)
    private byte[,,] InterleaveFields(byte[,,] field1, byte[,,] field2)
    {
        int fieldLines = field1.GetLength(0);
        int width = field1.GetLength(1);
        int height = Math.Min(PAL_VISIBLE_LINES, fieldLines * 2);

        byte[,,] frame = new byte[height, width, 3];

        int outLine = 0;
        for (int i = 0; i < fieldLines && outLine < height - 1; i++)
        {
            for (int j = 0; j < width; j++)
            {
                for (int c = 0; c < 3; c++)
                {
                    if (_fieldOrder == FieldOrder.BottomFieldFirst)
                    {
                        // BFF: top line from Field 2, next line from Field 1
                        frame[outLine, j, c] = field2[i, j, c];
                        frame[outLine + 1, j, c] = field1[i, j, c];
                    }
                    else
                    {
                        // TFF: top line from Field 1, next line from Field 2
                        frame[outLine, j, c] = field1[i, j, c];
                        frame[outLine + 1, j, c] = field2[i, j, c];
                    }
                }
            }
            outLine += 2;
        }
        return frame;
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

        int searchLength = Math.Min(videoSignal.Length, Math.Max(sampleRate / 4, samplesPerLine * 625));
        var segment = LightLowPass(videoSignal.Take(searchLength).ToArray());

        // auto polarity
        double segMin = segment.Min(), segMax = segment.Max();
        bool inverted = Math.Abs(segMax) > Math.Abs(segMin);
        if (inverted)
        {
            for (int j = 0; j < segment.Length; j++) segment[j] = -segment[j];
            (segMin, segMax) = (segment.Min(), segment.Max());
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

    private double[] AMDemodulation(Complex[] iqSamples)
    {
        // AM demodulation: extract magnitude of complex signal
        double[] demodulated = new double[iqSamples.Length];

        for (int i = 0; i < iqSamples.Length; i++)
        {
            demodulated[i] = iqSamples[i].Magnitude;
        }

        // Remove DC component
        double dcOffset = demodulated.Average();
        for (int i = 0; i < demodulated.Length; i++)
        {
            demodulated[i] -= dcOffset;
        }

        return demodulated;
    }


    private (double[] luminance, double[] chrominance) SeparateLumaChroma(double[] videoSignal, int sampleRate, int samplesPerLine)
    {
        // Luma LPF (System profile)
        double lumaCutoff = Math.Min(_profile.LumaCutoffHz, 0.45 * sampleRate); // safety vs Nyquist
        var lumaFilter = CreateLowPassFilter(lumaCutoff, sampleRate);
        double[] luminance = ApplyFilter(videoSignal, lumaFilter);

        // Chroma BPF (tighter around 4.4336 MHz to avoid sound at +6.5 MHz in PAL-D/K)
        double chromaLow = _profile.ChromaLowHz;
        double chromaHigh = _profile.ChromaHighHz;
        var chromaFilter = CreateBandPassFilter(chromaLow, chromaHigh, sampleRate);
        double[] chrominance = ApplyFilter(videoSignal, chromaFilter);

        return (luminance, chrominance);
    }

    private (double[] uComponent, double[] vComponent) DecodeChroma(double[] chrominance, int sampleRate, int samplesPerLine, int startLineOffset)
    {
        double[] uComponent = new double[chrominance.Length];
        double[] vComponent = new double[chrominance.Length];

        var pll = new ColorPll(sampleRate);

        for (int i = 0; i < chrominance.Length; i++)
        {
            int lineInBuffer = i / samplesPerLine;
            int sampleInLine = i % samplesPerLine;

            // Absolute line index across the full frame
            int absoluteLine = startLineOffset + lineInBuffer;

            // PAL V-axis alternates every line; absolute parity must be used across fields
            bool isVInverted = (absoluteLine % 2 != 0);

            // Get the stable, phase-locked reference carrier from the PLL
            var referenceCarrier = pll.GetReference(new Complex(chrominance[i], 0), sampleInLine, isVInverted);
            // Demodulate using the PLL's reference
            var demodulatedSample = new Complex(chrominance[i], 0) * Complex.Conjugate(referenceCarrier);

            // U demodulation (in-phase with subcarrier)
            uComponent[i] = demodulatedSample.Real;

            // For V, we must respect the PAL switch. We flip the sign back on inverted lines.
            vComponent[i] = isVInverted ? -demodulatedSample.Imaginary : demodulatedSample.Imaginary;
        }

        // Demodulated U & V Components:
        // |--0Hz----1.3MHz----|
        // |   Baseband U & V  |
        // The 1.3 MHz represents the baseband chroma bandwidth in PAL:
        // - U component bandwidth: ~1.3 MHz
        // - V component bandwidth: ~1.3 MHz

        // Result after multiplication:
        // - Desired: Baseband signal (0 - 1.3 MHz) ← What we want
        // - Unwanted: High frequency components around 8.86 MHz (4.43 + 4.43) ← Noise. To remove this we apply a low-pass filter.

        // Low-pass filter the demodulated components
        double chromaCutoff = 1.3e6;
        var chromaLPF = CreateLowPassFilter(chromaCutoff, sampleRate);
        uComponent = ApplyFilter(uComponent, chromaLPF);
        vComponent = ApplyFilter(vComponent, chromaLPF);

        return (uComponent, vComponent);
    }


    private byte[,,] ConvertYUVToRGB_BT709(double[] y, double[] u, double[] v, int samplesPerLine)
    {
        int width = Math.Min(720, samplesPerLine); // Standard PAL width
        int height = PAL_VISIBLE_LINES;

        byte[,,] rgbFrame = new byte[height, width, 3]; // RGB

        // BT.709 YUV to RGB conversion matrix
        // [R]   [1.0000   0.0000   1.5748] [Y]
        // [G] = [1.0000  -0.1873  -0.4681] [U]
        // [B]   [1.0000   1.8556   0.0000] [V]

        double[,] bt709Matrix = new double[3, 3]
        {
        { 1.0000,  0.0000,  1.5748 }, // R coefficients
        { 1.0000, -0.1873, -0.4681 }, // G coefficients  
        { 1.0000,  1.8556,  0.0000 }  // B coefficients
        };

        for (int row = 0; row < height; row++)
        {
            for (int col = 0; col < width; col++)
            {
                int index = row * samplesPerLine + col;

                if (index < y.Length)
                {
                    // Normalize and clamp input values
                    double yVal = Clamp(y[index], 0, 1);
                    double uVal = Clamp(u[index], -0.5, 0.5);
                    double vVal = Clamp(v[index], -0.5, 0.5);

                    // Create YUV vector
                    double[] yuvVector = { yVal, uVal, vVal };

                    // Matrix multiplication: RGB = Matrix × YUV
                    double r = bt709Matrix[0, 0] * yuvVector[0] +
                              bt709Matrix[0, 1] * yuvVector[1] +
                              bt709Matrix[0, 2] * yuvVector[2];

                    double g = bt709Matrix[1, 0] * yuvVector[0] +
                              bt709Matrix[1, 1] * yuvVector[1] +
                              bt709Matrix[1, 2] * yuvVector[2];

                    double b = bt709Matrix[2, 0] * yuvVector[0] +
                              bt709Matrix[2, 1] * yuvVector[1] +
                              bt709Matrix[2, 2] * yuvVector[2];

                    // Clamp and convert to byte values
                    rgbFrame[row, col, 0] = (byte)(Clamp(r, 0, 1) * 255); // R
                    rgbFrame[row, col, 1] = (byte)(Clamp(g, 0, 1) * 255); // G
                    rgbFrame[row, col, 2] = (byte)(Clamp(b, 0, 1) * 255); // B
                }
            }
        }

        return rgbFrame;
    }

    private byte[,,] ConvertYUVToRGB_BT601_Optimized(double[] y, double[] u, double[] v, int samplesPerLine)
    {
        int width = Math.Min(720, samplesPerLine);
        int height = y.Length / samplesPerLine; // derive from buffer length

        byte[,,] rgbFrame = new byte[height, width, 3];

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

                rgbFrame[row, col, 0] = (byte)(Clamp(r, 0, 1) * 255);
                rgbFrame[row, col, 1] = (byte)(Clamp(g, 0, 1) * 255);
                rgbFrame[row, col, 2] = (byte)(Clamp(b, 0, 1) * 255);
            }
        }
        return rgbFrame;
    }

    private double[] CreateLowPassFilter(double cutoffFreq, int sampleRate)
    {
        // Simple FIR low-pass filter
        int filterLength = 101;
        double[] filter = new double[filterLength];
        double fc = cutoffFreq / sampleRate;

        for (int i = 0; i < filterLength; i++)
        {
            int n = i - filterLength / 2;
            if (n == 0)
                filter[i] = 2 * fc;
            else
                filter[i] = Math.Sin(2 * Math.PI * fc * n) / (Math.PI * n);

            // Apply Hamming window
            filter[i] *= 0.54 - 0.46 * Math.Cos(2 * Math.PI * i / (filterLength - 1));
        }

        return filter;
    }

    private double[] CreateBandPassFilter(double lowFreq, double highFreq, int sampleRate)
    {
        // Create band-pass as difference of two low-pass filters
        var lpf1 = CreateLowPassFilter(highFreq, sampleRate);
        var lpf2 = CreateLowPassFilter(lowFreq, sampleRate);

        double[] bandPass = new double[lpf1.Length];
        for (int i = 0; i < bandPass.Length; i++)
        {
            bandPass[i] = lpf1[i] - lpf2[i];
        }

        return bandPass;
    }

    private double[] ApplyFilter(double[] signal, double[] filter)
    {
        int signalLength = signal.Length;
        int filterLength = filter.Length;
        double[] filtered = new double[signalLength];

        for (int i = 0; i < signalLength; i++)
        {
            double sum = 0;
            for (int j = 0; j < filterLength; j++)
            {
                int index = i - j;
                if (index >= 0 && index < signalLength)
                {
                    sum += signal[index] * filter[j];
                }
            }
            filtered[i] = sum;
        }

        return filtered;
    }

    private double Clamp(double value, double min, double max)
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
    private readonly double _lineDuration;
    private readonly double _burstStartTime; // Time after H-sync where burst starts
    private readonly double _burstDuration;  // Duration of the color burst

    // PLL state variables
    private double _phase = 0.0;
    private double _frequency = PALDecoder.PAL_COLOR_CARRIER_FREQ;
    private double _phaseErrorIntegrator = 0.0;

    // PLL loop filter gains (these may need tuning for noisy signals)
    private const double ProportionalGain = 0.1;
    private const double IntegralGain = 0.005;

    public ColorPll(int sampleRate)
    {
        _sampleRate = sampleRate;
        _lineDuration = PALDecoder.PAL_LINE_DURATION;
        // PAL spec: H-sync (4.7µs) + Breezeway (0.6µs) = 5.3µs
        _burstStartTime = 5.3e-6;
        _burstDuration = 2.25e-6;
    }

    // Generates one sample of the reference carrier and updates the PLL state
    public Complex GetReference(Complex chromaSample, int sampleIndexInLine, bool isVInverted)
    {
        double timeInLine = sampleIndexInLine / _sampleRate;

        // --- Phase Detection (only during the color burst) ---
        if (timeInLine >= _burstStartTime && timeInLine < _burstStartTime + _burstDuration)
        {
            // Generate the PLL's current estimate of the carrier
            var pllReference = Complex.FromPolarCoordinates(1, _phase);

            // The expected phase of the burst depends on the V-axis switch
            double expectedBurstPhase = isVInverted ? -Math.PI * 3 / 4 : -Math.PI / 4; // -135° or -45° for PAL swing
            var burstReference = Complex.FromPolarCoordinates(1, _phase + expectedBurstPhase);

            // Calculate phase error: how far off is our PLL from the actual burst?
            double phaseError = (chromaSample * Complex.Conjugate(burstReference)).Phase;

            // --- Loop Filter ---
            // Update the integrator (I-term)
            _phaseErrorIntegrator += phaseError * IntegralGain;

            // Update the frequency based on the P and I terms
            _frequency = PALDecoder.PAL_COLOR_CARRIER_FREQ + (phaseError * ProportionalGain) + _phaseErrorIntegrator;
        }

        // --- Numerically Controlled Oscillator (NCO) ---
        // Advance the phase for the next sample using the (potentially updated) frequency
        _phase += 2 * Math.PI * _frequency / _sampleRate;

        // Wrap phase to keep it within [-PI, PI]
        if (_phase > Math.PI) _phase -= 2 * Math.PI;
        if (_phase < -Math.PI) _phase += 2 * Math.PI;

        // Return the final, stable reference carrier for this sample
        return Complex.FromPolarCoordinates(1, _phase);
    }
}