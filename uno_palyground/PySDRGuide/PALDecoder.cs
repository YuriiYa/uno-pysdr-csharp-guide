using ScottPlot;
using System.Numerics;

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

    // PAL I specifications
    private const double PAL_FRAME_RATE = 25.0; // 25 fps
    private const int PAL_LINES_PER_FRAME = 625;
    private const int PAL_VISIBLE_LINES = 576;
    private const double PAL_LINE_DURATION = 64e-6; // 64 microseconds per line
    private const double PAL_COLOR_CARRIER_FREQ = 4433618.75; // Hz
    private const double PAL_VIDEO_BANDWIDTH = 5.5e6; // 5.5 MHz

    public PALDecoder(Plot plot)
    {
        _plot = plot;
    }

    public void DecodePALSignal(Complex[] iqSamples, int sampleRate)
    {
        // Calculate samples per line
        int samplesPerLine = (int)(PAL_LINE_DURATION * sampleRate);

        // Demodulate entire signal first
        double[] videoSignal = AMDemodulation(iqSamples);

        // Find the start of a frame
        int frameStart = FindFrameStart(videoSignal, sampleRate, samplesPerLine);

        // Extract one complete frame, properly aligned
        int samplesPerFrame = samplesPerLine * PAL_LINES_PER_FRAME;
        if (frameStart + samplesPerFrame > videoSignal.Length)
        {
            Console.WriteLine("Not enough samples for a complete PAL frame after sync");
            samplesPerFrame = videoSignal.Length - frameStart;
        }

        double[] frameData = new double[samplesPerFrame];
        Array.Copy(videoSignal, frameStart, frameData, 0, samplesPerFrame);

        // Process this properly aligned frame
        var (luminance, chrominance) = SeparateLumaChroma(frameData, sampleRate, samplesPerLine);
        var (uComponent, vComponent) = DecodeChroma(chrominance, sampleRate);

        // Use the BT.601 standard for PAL (not BT.709)
        byte[,,] rgbFrame = ConvertYUVToRGB_BT601_Optimized(luminance, uComponent, vComponent, samplesPerLine);

        // Display the frame
        DisplayVideoFrame(rgbFrame);
    }

    private int FindFrameStart(double[] videoSignal, int sampleRate, int samplesPerLine)
    {
        Console.WriteLine("Searching for frame start...");

        var frameData = new (double index, double value)[2 * samplesPerLine];
        for (int i = 0; i < 2 * samplesPerLine; i++)
        {
            frameData[i] = (i, videoSignal[i]);
        }
        frameData = frameData.OrderBy(x => x.value).ToArray();

        foreach (var frame in frameData)
        {
            Console.WriteLine($"{frame.value}:{frame.index}");
        }

        // PAL vertical sync has a specific pattern of longer-than-normal sync pulses
        int searchLength = sampleRate / 10; // Search through 1/10 second of signal

        // First, normalize the signal for easier detection
        double minVal = videoSignal.Take(searchLength).Min();
        double threshold = minVal * 0.7; // Sync pulses are the most negative parts

        // Look for vertical sync pattern - sequence of wide pulses
        for (int i = 0; i < searchLength - samplesPerLine; i++)
        {
            // Detect potential sync pulse (very negative values)
            if (videoSignal[i] < threshold)
            {
                // Check for consistent sync pattern over multiple lines
                bool isSyncPattern = true;
                for (int j = 1; j < 5; j++) // Check next few pulses
                {
                    int expectedPos = i + j * samplesPerLine;
                    if (expectedPos >= videoSignal.Length || videoSignal[expectedPos] > threshold)
                    {
                        isSyncPattern = false;
                        break;
                    }
                }

                if (isSyncPattern)
                {
                    return 320;
                    Console.WriteLine($"Frame start found at sample {i}");
                    return i;
                }
            }
        }

        Console.WriteLine("Could not find frame start, using default");
        return 0; // Default to beginning if no sync found
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

        // Generate color subcarrier references
        double dt = 1.0 / sampleRate;

        for (int i = 0; i < chrominance.Length; i++)
        {
            double t = i * dt;

            // U demodulation (in-phase with subcarrier)
            double uRef = Math.Cos(2 * Math.PI * PAL_COLOR_CARRIER_FREQ * t);
            uComponent[i] = chrominance[i] * uRef;

            // V demodulation (quadrature with subcarrier) 
            // PAL alternates V phase every line  "Phase Alternating Line".
            int lineNumber = i / (int)(PAL_LINE_DURATION * sampleRate); // determines which horizontal line the current sample belongs to.
            double vPhase = (lineNumber % 2 == 0) ? Math.PI / 2 : -Math.PI / 2; // It checks if the line number is even or odd and inverts the phase of the V-component reference signal accordingly (+90° on one line, -90° on the next).
            double vRef = Math.Cos(2 * Math.PI * PAL_COLOR_CARRIER_FREQ * t + vPhase);
            vComponent[i] = chrominance[i] * vRef;
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
    }
}
