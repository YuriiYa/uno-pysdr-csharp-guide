using uno_palyground.PySDRGuide.Pipeline;
using Xunit;

namespace PALDecoder.Tests;

public class LumaChromaSeparatorTests
{
    private const int SampleRate = 10_000_000; // 10 MHz
    private static SystemProfile DefaultProfile => SystemProfile.For(TvSystem.PAL_DK);

    // 5.8 DC-only input → chroma output near zero, luma output near DC value
    [Fact]
    public void Separate_DcInput_ChromaNearZeroLumaPassesThrough()
    {
        const int length = 2000;
        double[] field = new double[length];
        double[] lumaOut = new double[length];
        double[] chromaOut = new double[length];

        const double dcValue = 0.5;
        Array.Fill(field, dcValue);

        var sep = new LumaChromaSeparator(DefaultProfile, SampleRate);
        sep.Separate(field, lumaOut, chromaOut);

        // After FIR settling (skip warm-up) chroma should be near zero
        int skip = 100; // skip FIR warm-up samples
        double maxChroma = chromaOut.Skip(skip).Max(Math.Abs);
        Assert.InRange(maxChroma, 0, 0.01); // chroma << DC level

        // Luma should converge to ≈ dcValue
        double avgLuma = lumaOut.Skip(skip).Average();
        Assert.InRange(avgLuma, dcValue - 0.05, dcValue + 0.05);
    }

    // 5.9 4.43 MHz tone → chroma amplitude significantly larger than luma amplitude
    [Fact]
    public void Separate_ChromaTone_ChromaLargerThanLuma()
    {
        const int length = 4000;
        double[] field = new double[length];
        double[] lumaOut = new double[length];
        double[] chromaOut = new double[length];

        double chromaFreq = PalConstants.PAL_COLOR_CARRIER_FREQ; // 4.43361875 MHz
        for (int i = 0; i < length; i++)
            field[i] = Math.Sin(2 * Math.PI * chromaFreq / SampleRate * i);

        var sep = new LumaChromaSeparator(DefaultProfile, SampleRate);
        sep.Separate(field, lumaOut, chromaOut);

        // Skip FIR warm-up
        int skip = 200;
        double chromaRms = Math.Sqrt(chromaOut.Skip(skip).Average(x => x * x));
        double lumaRms = Math.Sqrt(lumaOut.Skip(skip).Average(x => x * x));

        // Chroma at 4.43 MHz should pass the BPF while the luma LPF attenuates it near its cutoff
        // (at 10 MHz sample rate luma cutoff = 4.5 MHz, so 4.43 MHz is right at the edge)
        Assert.True(chromaRms > lumaRms,
            $"Expected chromaRms ({chromaRms:F4}) > lumaRms ({lumaRms:F4})");
    }
}
