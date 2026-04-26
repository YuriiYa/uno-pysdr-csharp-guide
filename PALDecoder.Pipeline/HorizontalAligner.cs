namespace uno_palyground.PySDRGuide.Pipeline;

public class HorizontalAligner
{
    public int EstimateOffset(double[] videoSignal, int frameStart, int samplesPerLine, int sampleRate, int linesToUse = 24)
    {
        // LPF group delay: (taps-1)/2 samples; must match LightLpfTaps in PalConstants.
        int gd = (PalConstants.LightLpfTaps - 1) / 2;

        int searchCols = Math.Min(samplesPerLine, (int)Math.Round(PalConstants.HSyncSearchSeconds * sampleRate));
        int hsyncSamples = Math.Max(8, (int)Math.Round(PalConstants.HSyncSeconds * sampleRate));
        int desiredActiveCol = (int)Math.Round((PalConstants.HSyncSeconds + PalConstants.BackPorchSeconds) * sampleRate);

        List<int> activeStarts = new();

        for (int ln = 0; ln < linesToUse; ln++)
        {
            int lineStart = frameStart + ln * samplesPerLine;
            if (lineStart + searchCols > videoSignal.Length) break;

            double[] win = new double[searchCols];
            Array.Copy(videoSignal, lineStart, win, 0, searchCols);
            win = LightLowPass(win);

            (double wMin, double wMax) = MinMax(win);
            if (Math.Abs(wMax) > Math.Abs(wMin))
            {
                for (int i = 0; i < win.Length; i++) win[i] = -win[i];
                (wMin, wMax) = (-wMax, -wMin);
            }

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
                double med = Median(win);
                double mad = Median(win.Select(v => Math.Abs(v - med)).ToArray());
                if (mad <= PalConstants.MadNearZeroThreshold) mad = (wMax - wMin) / PalConstants.MadFallbackRangeDivisor;
                // TODO: check the correctness: HSyncMadMultiplier threshold is empirically tuned, not derived from spec
                double thr = med - PalConstants.HSyncMadMultiplier * mad;

                int rs = bestStart;
                while (rs > 0 && win[rs - 1] < thr) rs--;
                int re = Math.Min(N - 1, bestStart + L - 1);
                while (re + 1 < N && win[re + 1] < thr) re++;

                int len = re - rs + 1;
                if (len >= (int)(0.5 * hsyncSamples) && len <= (int)(2.0 * hsyncSamples))
                {
                    int syncStartCol = Math.Max(0, rs - gd);
                    int activeCol = syncStartCol + (int)Math.Round((PalConstants.HSyncSeconds + PalConstants.BackPorchSeconds) * sampleRate);
                    activeStarts.Add(activeCol);
                    added = true;
                }
            }

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
                    int syncEdge = Math.Max(0, edgeIdx - gd);
                    int activeCol = syncEdge + (int)Math.Round((PalConstants.HSyncSeconds + PalConstants.BackPorchSeconds) * sampleRate);
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
        int offset = desiredActiveCol - medianActiveCol;
        offset = Math.Max(-samplesPerLine / 2, Math.Min(samplesPerLine / 2, offset));

        Console.WriteLine($"Auto horizontal offset = {offset} (median={medianActiveCol}, desired={desiredActiveCol}, gd={gd}, lines={activeStarts.Count})");
        return offset;
    }

    internal static double Median(double[] a)
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
