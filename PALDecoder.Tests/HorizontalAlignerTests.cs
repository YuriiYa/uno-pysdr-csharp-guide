using uno_palyground.PySDRGuide.Pipeline;
using Xunit;

namespace PALDecoder.Tests;

public class HorizontalAlignerTests
{
    // Use 20 MHz so hsyncSamples (94) >> gd (25): LPF-smeared pulse width stays within the 2x limit
    private const int SampleRate = 20_000_000;
    private static int SamplesPerLine => (int)(64e-6 * SampleRate); // 1280

    // 5.5 Synthetic line with a 4.7 us negative H-sync at a known column
    [Fact]
    public void EstimateOffset_SyntheticHSync_ReturnsOffsetWithinTwoSamples()
    {
        int samplesPerLine = SamplesPerLine;
        int hsyncSamples = (int)Math.Round(4.7e-6 * SampleRate);  // 94

        // Build 30 lines with H-sync at column 0 and blanking at 0.0
        // Using -1.0 / 0.0 ensures |wMin| >> |wMax| so no polarity inversion occurs
        int linesToUse = 30;
        double[] signal = new double[samplesPerLine * (linesToUse + 2)];
        for (int ln = 0; ln < linesToUse; ln++)
        {
            int lineStart = ln * samplesPerLine;
            for (int s = 0; s < hsyncSamples; s++)
                signal[lineStart + s] = -1.0;
        }

        var aligner = new HorizontalAligner();
        int offset = aligner.EstimateOffset(signal, frameStart: 0, samplesPerLine: samplesPerLine,
            sampleRate: SampleRate, linesToUse: linesToUse);

        // desiredActiveCol = round((4.7+5.8)us * 20MHz) = 210
        // Sync at col 0, detected activeCol ~= 210, offset ~= 0
        Assert.InRange(offset, -2, 2);
    }
}
