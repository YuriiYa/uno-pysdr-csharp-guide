using System.Numerics;
using System.Runtime.CompilerServices;
using System.Collections.Concurrent;

namespace uno_palyground.PySDRGuide.Pipeline;

public class ChromaDecoder
{
    private readonly int _sampleRate;
    private readonly int _chromaBasebandLpfTaps;
    private static readonly ConcurrentDictionary<string, double[]> _filterCache = new();

    private double[]? _uScratch;
    private double[]? _vScratch;

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

    private static void EnsureBuffer(ref double[]? buf, int len)
    {
        if (buf == null || buf.Length < len) buf = new double[len];
    }

    public ChromaDecoder(int sampleRate, int chromaBasebandLpfTaps = 63)
    {
        _sampleRate = sampleRate;
        _chromaBasebandLpfTaps = chromaBasebandLpfTaps;
    }

    public void Decode(double[] chroma, int sampleRate, int samplesPerLine, int startLineOffset, double[] uOut, double[] vOut)
    {
        int len = chroma.Length;
        EnsureBuffer(ref _uScratch, len);
        EnsureBuffer(ref _vScratch, len);

        var chromaLPF = CreateLowPassFilter(1.3e6, _sampleRate, _chromaBasebandLpfTaps);
        var revLPF = FilterReverseCache.Get(chromaLPF);

        FusedChromaDemodAndFilter(chroma, _sampleRate, samplesPerLine, startLineOffset, chromaLPF, revLPF, _uScratch!, _vScratch!, uOut, vOut);
    }

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
}
