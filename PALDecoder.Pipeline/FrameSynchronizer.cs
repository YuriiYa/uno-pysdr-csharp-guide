namespace uno_palyground.PySDRGuide.Pipeline;

public class FrameSynchronizer
{
    public int FindFrameStart(double[] videoSignal, int sampleRate, int samplesPerLine)
    {
        Console.WriteLine("Searching for frame start (robust VBI detection)...");
        if (videoSignal == null || videoSignal.Length < samplesPerLine * 50)
        {
            Console.WriteLine("Signal too short for VBI search, fallback to 0.");
            return 0;
        }

        int searchLength = Math.Min(videoSignal.Length, Math.Max(sampleRate / 4, samplesPerLine * PalConstants.PAL_LINES_PER_FRAME));
        var segment = LightLowPass(videoSignal.Take(searchLength).ToArray());

        (double segMin, double segMax) = MinMax(segment);
        bool inverted = Math.Abs(segMax) > Math.Abs(segMin);
        if (inverted)
        {
            for (int j = 0; j < segment.Length; j++) segment[j] = -segment[j];
            (segMin, segMax) = MinMax(segment);
        }

        double median = Median(segment);
        double mad = Median(segment.Select(x => Math.Abs(x - median)).ToArray());
        if (mad <= PalConstants.MadNearZeroThreshold) mad = (segMax - segMin) / PalConstants.MadFallbackRangeDivisor;
        // TODO: check the correctness: VSyncMadMultiplier threshold is empirically tuned
        double syncThreshold = median - PalConstants.VSyncMadMultiplier * mad;

        int broadMin = (int)Math.Round(PalConstants.BroadPulseMinSeconds * sampleRate);
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
                    int refined = RefineToFirstActiveLine(videoSignal, runStart, sampleRate, samplesPerLine);
                    Console.WriteLine($"VBI detected (inv={inverted}). Field1 active ≈ {refined}");
                    return refined;
                }
            }
        }

        Console.WriteLine("Could not find VBI pattern, using default of 0.");
        return 0;
    }

    private int RefineToFirstActiveLine(double[] video, int vsyncFirstBroadStart, int sampleRate, int samplesPerLine)
    {
        int hsyncSamples = (int)Math.Round(PalConstants.HSyncSeconds * sampleRate);
        int breezewaySamples = (int)Math.Round(PalConstants.BreezewaySeconds * sampleRate);
        int burstLenSamples = (int)Math.Round(PalConstants.BurstDurationSeconds * sampleRate);
        int backPorchToBurst = hsyncSamples + breezewaySamples;
        int desiredActiveCol = (int)Math.Round((PalConstants.HSyncSeconds + PalConstants.BackPorchSeconds) * sampleRate);

        int firstLineIdx = Math.Max(0, (vsyncFirstBroadStart / samplesPerLine) - 1);
        int lastLineIdx = Math.Min(firstLineIdx + PalConstants.BurstSearchLinesAfterVSync, (video.Length / samplesPerLine) - 1);

        List<(int line, double power)> candidates = new();
        for (int ln = firstLineIdx; ln <= lastLineIdx; ln++)
        {
            int lineStart = ln * samplesPerLine;
            if (lineStart + backPorchToBurst + burstLenSamples >= video.Length) break;

            int burstStart = lineStart + backPorchToBurst;
            double p = GoertzelPower(video, burstStart, burstLenSamples, PalConstants.PAL_COLOR_CARRIER_FREQ, sampleRate);
            candidates.Add((ln, p));
        }
        if (candidates.Count == 0) return vsyncFirstBroadStart + (int)Math.Round(PalConstants.VbiLinesToActiveVideo * samplesPerLine);

        double[] powers = candidates.Select(t => t.power).ToArray();
        double m = Median(powers);
        double md = Median(powers.Select(v => Math.Abs(v - m)).ToArray());
        if (md <= PalConstants.BurstMadNearZeroThreshold) md = (powers.Max() - powers.Min()) / PalConstants.MadFallbackRangeDivisor;
        double thr = m + PalConstants.VSyncMadMultiplier * md;

        for (int k = 0; k < candidates.Count; k++)
        {
            bool strong = candidates[k].power > thr;
            bool nextStrong = (k + 1 < candidates.Count) && (candidates[k + 1].power > thr);
            bool next2Strong = (k + 2 < candidates.Count) && (candidates[k + 2].power > thr);
            if (strong && (nextStrong || next2Strong))
            {
                int ln = candidates[k].line;
                return ln * samplesPerLine + desiredActiveCol;
            }
        }

        var best = candidates.OrderByDescending(t => t.power).First();
        return best.line * samplesPerLine + desiredActiveCol;
    }

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

    internal static double[] LightLowPass(double[] x)
    {
        int N = PalConstants.LightLpfTaps;
        if (x.Length < N) return x;
        double[] w = new double[N];
        for (int i = 0; i < N; i++) w[i] = PalConstants.HammingAlpha - PalConstants.HammingBeta * Math.Cos(2 * Math.PI * i / (N - 1));
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

    internal static (double min, double max) MinMax(ReadOnlySpan<double> s)
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
}
