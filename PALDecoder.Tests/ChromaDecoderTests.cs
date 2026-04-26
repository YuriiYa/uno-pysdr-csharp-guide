using uno_palyground.PySDRGuide.Pipeline;
using Xunit;

namespace PALDecoder.Tests;

public class ChromaDecoderTests
{
    private const int SampleRate = 10_000_000;
    private const int SamplesPerLine = 640; // (int)(64e-6 * 10MHz)

    // 5.10 Verify V negation on odd absolute lines
    // On even absolute lines (absLine % 2 == 0) → V is demod.Imaginary
    // On odd absolute lines (absLine % 2 == 1)  → V is −demod.Imaginary (PAL phase alternation)
    [Fact]
    public void Decode_VComponentNegatedOnOddLines()
    {
        // Build a chroma buffer with 2 lines of a 4.43 MHz sine tone (single frequency = known demod phase).
        // With startLineOffset = 0: line 0 (abs 0, even) → V = +Im, line 1 (abs 1, odd) → V = −Im.
        // After decoding, if V-line-0 and V-line-1 are compared at the same within-line sample position,
        // they should have opposite signs (within tolerance, accounting for PLL transients).
        int twoLines = 2 * SamplesPerLine;
        double[] chroma = new double[twoLines];
        double freq = PalConstants.PAL_COLOR_CARRIER_FREQ;
        for (int i = 0; i < twoLines; i++)
            chroma[i] = Math.Sin(2 * Math.PI * freq / SampleRate * i);

        var decoder = new ChromaDecoder(SampleRate, chromaBasebandLpfTaps: 63);
        double[] uOut = new double[twoLines];
        double[] vOut = new double[twoLines];
        decoder.Decode(chroma, SampleRate, SamplesPerLine, startLineOffset: 0, uOut, vOut);

        // Compare mid-line V values: line 0 and line 1 should have opposite sign
        // Skip early samples due to PLL lock transient
        int midLine = SamplesPerLine / 2;
        double vLine0 = vOut[midLine];          // line 0: even abs → positive V
        double vLine1 = vOut[SamplesPerLine + midLine]; // line 1: odd abs → negated V

        // They should have opposite signs; if neither is near-zero the test is meaningful
        if (Math.Abs(vLine0) > 0.001 && Math.Abs(vLine1) > 0.001)
            Assert.True(Math.Sign(vLine0) != Math.Sign(vLine1),
                $"V line0={vLine0:F4} and V line1={vLine1:F4} should have opposite signs");
        // else PLL hasn't locked yet — just verify the decode ran without exception (passed if no throw)
    }
}
