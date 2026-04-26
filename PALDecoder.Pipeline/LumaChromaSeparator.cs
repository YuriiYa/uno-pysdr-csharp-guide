using System.Collections.Concurrent;
using System.Runtime.CompilerServices;
using System.Numerics;

namespace uno_palyground.PySDRGuide.Pipeline;

public class LumaChromaSeparator
{
    private readonly SystemProfile _profile;
    private readonly int _sampleRate;
    private readonly int _lumaLpfTaps;
    private readonly int _chromaSeparationTaps;
    private static readonly ConcurrentDictionary<string, double[]> _filterCache = new();

    private static class FilterReverseCache
    {
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

    public LumaChromaSeparator(SystemProfile profile, int sampleRate, int lumaLpfTaps = 81, int chromaSeparationTaps = 81)
    {
        _profile = profile;
        _sampleRate = sampleRate;
        _lumaLpfTaps = lumaLpfTaps;
        _chromaSeparationTaps = chromaSeparationTaps;
    }

    public void Separate(double[] field, double[] lumaOut, double[] chromaOut)
    {
        double lumaCutoff = Math.Min(_profile.LumaCutoffHz, PalConstants.LumaCutoffNyquistMargin * _sampleRate);
        var lumaFilter = CreateLowPassFilter(lumaCutoff, _sampleRate, _lumaLpfTaps);
        var chromaFilter = CreateBandPassFilter(_profile.ChromaLowHz, _profile.ChromaHighHz, _sampleRate, _chromaSeparationTaps);
        ApplyTwoFiltersStaticToDest(field, lumaFilter, chromaFilter, lumaOut, chromaOut);
    }

    private double[] CreateLowPassFilter(double cutoffFreq, int sampleRate, int filterLength = 101)
    {
        string filterKey = $"LPF_{cutoffFreq}_{sampleRate}_{filterLength}";
        return _filterCache.GetOrAdd(filterKey, static key =>
        {
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
                v *= PalConstants.HammingAlpha - PalConstants.HammingBeta * Math.Cos(2 * Math.PI * i / (fl - 1));
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

    private static void ApplyTwoFiltersStaticToDest(double[] signal, double[] filterA, double[] filterB, double[] destA, double[] destB)
    {
        int n = signal.Length;
        int tapsA = filterA.Length;
        int tapsB = filterB.Length;
        if (n == 0) return;
        if (tapsA != tapsB)
        {
            ApplyFilterStaticToDest(signal, filterA, destA);
            ApplyFilterStaticToDest(signal, filterB, destB);
            return;
        }
        int taps = tapsA;
        int warm = taps - 1;
        int maxWarm = Math.Min(warm, n - 1);
        double[] revA = FilterReverseCache.Get(filterA);
        double[] revB = ReferenceEquals(filterA, filterB) ? revA : FilterReverseCache.Get(filterB);
        int simdWidth = Vector<double>.Count;
        bool useSimd = Vector.IsHardwareAccelerated && taps >= simdWidth * 2;

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

    private static void ApplyFilterStaticToDest(double[] signal, double[] filter, double[] dest)
    {
        int n = signal.Length;
        int taps = filter.Length;
        if (n == 0) return;
        int warm = taps - 1;
        int maxWarm = Math.Min(warm, n - 1);
        for (int i = 0; i <= maxWarm; i++)
        {
            double sum = 0;
            int maxTap = Math.Min(taps - 1, i);
            for (int j = 0; j <= maxTap; j++) sum += signal[i - j] * filter[j];
            dest[i] = sum;
        }
        if (maxWarm == n - 1) return;
        double[] filtRev = FilterReverseCache.Get(filter);
        int simdWidth = Vector<double>.Count;
        bool useSimd = Vector.IsHardwareAccelerated && taps >= simdWidth * 2;
        for (int i = warm; i < n; i++)
        {
            int start = i - taps + 1;
            double sum;
            if (useSimd)
            {
                sum = 0;
                int k = 0;
                int limit = taps - simdWidth;
                Vector<double> vacc = Vector<double>.Zero;
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
}
