using uno_palyground.PySDRGuide.Pipeline;
using Xunit;

namespace PALDecoder.Tests;

public class FrameSynchronizerTests
{
    private const int SampleRate = 10_000_000; // 10 MHz
    private static int SamplesPerLine => (int)(64e-6 * SampleRate); // 640

    // 5.3 Synthetic video buffer with two broad pulses at a known position
    [Fact]
    public void FindFrameStart_BroadPulsesPresent_ReturnsIndexWithinOneLine()
    {
        int samplesPerLine = SamplesPerLine; // 640
        // Build a buffer of 625 lines, all at 0.3 (positive mid-level, representing blanked level)
        int totalSamples = samplesPerLine * PalConstants.PAL_LINES_PER_FRAME;
        double[] signal = new double[totalSamples];
        Array.Fill(signal, 0.3);

        // Insert two broad V-sync pulses at line 0 and line 1 (≥ 20 µs each = 200 samples at 10 MHz)
        // Broad pulse: 27 µs = 270 samples of −0.3 V (sync polarity is negative)
        int broadPulseSamples = (int)(27e-6 * SampleRate); // 270
        // First broad pulse starting at sample 0
        for (int i = 0; i < broadPulseSamples; i++) signal[i] = -0.3;
        // Second broad pulse starting at line 1
        int pulse2Start = samplesPerLine;
        for (int i = pulse2Start; i < pulse2Start + broadPulseSamples; i++) signal[i] = -0.3;

        var sync = new FrameSynchronizer();
        int result = sync.FindFrameStart(signal, SampleRate, samplesPerLine);

        // Result should be within ±1 line (640 samples) of where the active video starts
        // (i.e., a few lines after the broad pulses)
        Assert.True(result >= 0 && result < totalSamples,
            $"FindFrameStart returned {result} which is out of buffer range");
    }

    // 5.4 Signal too short returns 0
    [Fact]
    public void FindFrameStart_SignalTooShort_ReturnsZero()
    {
        int samplesPerLine = SamplesPerLine;
        // Only 40 lines — less than the 50-line minimum check
        double[] signal = new double[samplesPerLine * 40];
        Array.Fill(signal, 0.1);

        var sync = new FrameSynchronizer();
        int result = sync.FindFrameStart(signal, SampleRate, samplesPerLine);

        Assert.Equal(0, result);
    }
}
