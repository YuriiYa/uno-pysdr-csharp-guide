using ScottPlot;
using System.Numerics;
using Microsoft.UI.Dispatching;

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

public class PALDecoder
{
    private readonly Plot _plot;
    private readonly DispatcherQueue _dispatcherQueue;

    // PAL I specifications
    public const double PAL_FRAME_RATE = 25.0; // 25 fps
    public const int PAL_LINES_PER_FRAME = 625;
    public const int PAL_VISIBLE_LINES = 576;
    public const double PAL_LINE_DURATION = 64e-6; // 64 microseconds per line
    public const double PAL_COLOR_CARRIER_FREQ = 4433618.75; // Hz
    public const double PAL_VIDEO_BANDWIDTH = 5.5e6; // 5.5 MHz

    public PALDecoder(Plot plot, DispatcherQueue dispatcherQueue)
    {
        _plot = plot;
        _dispatcherQueue = dispatcherQueue;
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
               var skipUntil = frameStart + 80;
               int numberOfFrames = (videoSignal.Length - skipUntil) / samplesPerFrame;
               for (int i = 0; i < numberOfFrames; i++)
               {
                   Array.Copy(videoSignal, skipUntil + i * samplesPerFrame, frameData, 0, samplesPerFrame);

                   // Separate fields
                   int fieldLines = PAL_VISIBLE_LINES / 2; // 288
                                                           // After you fill frameData with exactly one PAL frame (625 lines worth) OR one visible frame window,
                                                           // split fields as contiguous blocks, not every-other-line.
                   int linesPerFieldAll = PAL_LINES_PER_FRAME / 2;   // 312
                   int fieldLinesVis = PAL_VISIBLE_LINES / 2;        // 288

                   // If frameStart already points to Field-1 active top, use 0 and +linesPerFieldAll.
                   // Otherwise add per-field VBI offsets (~22–25 lines) to land on active video:
                   const int VBI_LINES_FIELD1 = 23;
                   const int VBI_LINES_FIELD2 = 23;

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

                   var (u1, v1) = DecodeChroma(chr1, sampleRate);
                   var (u2, v2) = DecodeChroma(chr2, sampleRate);

                   // Interleave fields for display
                   byte[,,] rgbFrame = InterleaveFields(
                       ConvertYUVToRGB_BT601_Optimized(lum1, u1, v1, samplesPerLine),
                       ConvertYUVToRGB_BT601_Optimized(lum2, u2, v2, samplesPerLine)
                   );

                   DisplayVideoFrame(rgbFrame);

               }
           }, TaskCreationOptions.LongRunning);
    }

    // Interleave two fields into a single frame (PAL interlacing)
    private byte[,,] InterleaveFields(byte[,,] field1, byte[,,] field2)
    {
        int fieldLines = field1.GetLength(0);
        int width = field1.GetLength(1);
        int height = PAL_VISIBLE_LINES; // Only visible lines
        //int height = (PAL_VISIBLE_LINES / 4 + 10) * 2;
        byte[,,] frame = new byte[height, width, 3];

        int outLine = 0;
        for (int i = 0; i < fieldLines && outLine < height - 1; i++)
        {
            for (int j = 0; j < width; j++)
            {
                for (int c = 0; c < 3; c++)
                {
                    frame[outLine, j, c] = field1[i, j, c];     // Even lines
                    frame[outLine + 1, j, c] = field2[i, j, c]; // Odd lines
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

        // Search up to ~0.25 s (enough for multiple frames)
        int searchLength = Math.Min(videoSignal.Length, Math.Max(sampleRate / 4, samplesPerLine * 625));
        var segment = videoSignal.Take(searchLength).ToArray();

        // 1) Optional light low-pass (reduces chroma ripple). 51-tap Hamming MA.
        segment = LightLowPass(segment);

        // 2) Auto-polarity: ensure sync pulses are negative
        double segMin = segment.Min();
        double segMax = segment.Max();
        bool inverted = Math.Abs(segMax) > Math.Abs(segMin); // peaks stronger than troughs → invert
        if (inverted)
        {
            for (int j = 0; j < segment.Length; j++) segment[j] = -segment[j];
            (segMin, segMax) = (segment.Min(), segment.Max());
        }

        // 3) Robust threshold from median and MAD
        double median = Median(segment);
        double mad = Median(segment.Select(x => Math.Abs(x - median)).ToArray());
        if (mad <= 1e-12) mad = (segMax - segMin) / 50.0; // fallback
                                                          // Sync pulses are well below blanking; go several MADs below median
        double syncThreshold = median - 3.0 * mad;

        int broadMin = (int)Math.Round(20e-6 * sampleRate); // ≥20 µs → broad pulse
        int i = 0;

        while (i < segment.Length - samplesPerLine)
        {
            // Find start of a run below threshold
            while (i < segment.Length && !(segment[i] < syncThreshold)) i++;
            if (i >= segment.Length) break;

            int runStart = i;
            while (i < segment.Length && segment[i] < syncThreshold) i++;
            int runLen = i - runStart;

            if (runLen >= broadMin)
            {
                // Validate there is another broad run shortly after (within a couple of lines)
                int j = i;
                int endSearch = Math.Min(segment.Length, runStart + 4 * samplesPerLine);
                bool secondBroadFound = false;
                while (j < endSearch)
                {
                    // skip above threshold
                    while (j < endSearch && !(segment[j] < syncThreshold)) j++;
                    int r2s = j;
                    while (j < endSearch && segment[j] < syncThreshold) j++;
                    int r2len = j - r2s;

                    if (r2len >= broadMin)
                    {
                        secondBroadFound = true;
                        break;
                    }
                }

                if (secondBroadFound)
                {
                    int frameStart = runStart + (int)Math.Round(22.5 * samplesPerLine);
                    Console.WriteLine($"VBI detected (inv={inverted}, thr={syncThreshold:E3}). Field1 active ≈ {frameStart}");
                    return frameStart;
                }
            }
            // continue scanning
        }

        Console.WriteLine("Could not find VBI pattern, using default of 0.");
        return 0;
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
        // Low-pass filter for luminance (cutoff ~3 MHz)
        double lumaCutoff = 3e6;
        var lumaFilter = CreateLowPassFilter(lumaCutoff, sampleRate);
        double[] luminance = ApplyFilter(videoSignal, lumaFilter);

        // Band-pass filter for chrominance (around 4.43 MHz)
        double chromaLow = 3.5e6;
        double chromaHigh = 5.5e6;
        var chromaFilter = CreateBandPassFilter(chromaLow, chromaHigh, sampleRate);
        double[] chrominance = ApplyFilter(videoSignal, chromaFilter);

        return (luminance, chrominance);
    }

    private (double[] uComponent, double[] vComponent) DecodeChroma(double[] chrominance, int sampleRate)
    {
        double[] uComponent = new double[chrominance.Length];
        double[] vComponent = new double[chrominance.Length];

        int samplesPerLine = (int)(PAL_LINE_DURATION * sampleRate);
        var pll = new ColorPll(sampleRate);

        for (int i = 0; i < chrominance.Length; i++)
        {
            int lineNumber = i / samplesPerLine;
            int sampleInLine = i % samplesPerLine;

            // The V-axis phase is inverted on every other line
            bool isVInverted = (lineNumber % 2 != 0);

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
        double chromaCutoff = 1.3e6; // Chroma bandwidth
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
        int height = PAL_VISIBLE_LINES;

        byte[,,] rgbFrame = new byte[height, width, 3];

        // BT.601 coefficients (precomputed for performance)
        // Based on: Kr = 0.299, Kg = 0.587, Kb = 0.114
        const double c_rv = 1.402;    // R from V coefficient
        const double c_gu = -0.344;   // G from U coefficient  
        const double c_gv = -0.714;   // G from V coefficient
        const double c_bu = 1.772;    // B from U coefficient

        for (int row = 0; row < height; row++)
        {
            for (int col = 0; col < width; col++)
            {
                int index = row * samplesPerLine + col;

                if (index < y.Length)
                {
                    double yVal = Clamp(y[index], 0, 1);
                    double uVal = Clamp(u[index], -0.5, 0.5);
                    double vVal = Clamp(v[index], -0.5, 0.5);

                    // BT.601 conversion with precomputed coefficients
                    double r = yVal + c_rv * vVal;
                    double g = yVal + c_gu * uVal + c_gv * vVal;
                    double b = yVal + c_bu * uVal;

                    rgbFrame[row, col, 0] = (byte)(Clamp(r, 0, 1) * 255); // R
                    rgbFrame[row, col, 1] = (byte)(Clamp(g, 0, 1) * 255); // G
                    rgbFrame[row, col, 2] = (byte)(Clamp(b, 0, 1) * 255); // B
                }
            }
        }

        return rgbFrame;
    }

    private byte[,,] ConvertYUVToRGB(double[] y, double[] u, double[] v, int samplesPerLine)
    {
        int width = Math.Min(720, samplesPerLine); // Standard PAL width
        int height = PAL_VISIBLE_LINES;

        byte[,,] rgbFrame = new byte[height, width, 3]; // RGB

        for (int row = 0; row < height; row++)
        {
            for (int col = 0; col < width; col++)
            {
                int index = row * samplesPerLine + col;

                if (index < y.Length)
                {
                    // YUV to RGB conversion
                    double yVal = Clamp(y[index], 0, 1);
                    double uVal = Clamp(u[index], -0.5, 0.5);
                    double vVal = Clamp(v[index], -0.5, 0.5);

                    // ITU-R BT.601 conversion
                    double r = yVal + 1.402 * vVal;
                    double g = yVal - 0.344 * uVal - 0.714 * vVal;
                    double b = yVal + 1.772 * uVal;

                    rgbFrame[row, col, 0] = (byte)(Clamp(r, 0, 1) * 255); // R
                    rgbFrame[row, col, 1] = (byte)(Clamp(g, 0, 1) * 255); // G
                    rgbFrame[row, col, 2] = (byte)(Clamp(b, 0, 1) * 255); // B
                }
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