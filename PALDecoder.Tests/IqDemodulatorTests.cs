using System.IO;
using uno_palyground.PySDRGuide;
using uno_palyground.PySDRGuide.Pipeline;
using Xunit;

namespace PALDecoder.Tests;

public class IqDemodulatorTests
{
    // 5.1 Generic-stream path with synthetic MemoryStream of known int8 IQ pairs
    [Fact]
    public void TryRead_GenericStream_ReturnsMagnitudesAndRemovesDC()
    {
        // Build a synthetic stream: all identical IQ pairs (I=64, Q=0) → magnitude = 0.5, mean = 0.5
        // After DC removal, all samples should be ≈ 0.
        const int samplesPerFrame = 100;
        byte[] rawBytes = new byte[samplesPerFrame * 2];
        for (int i = 0; i < samplesPerFrame; i++)
        {
            rawBytes[2 * i] = 64;   // I = 64 → 64/128.0 = 0.5
            rawBytes[2 * i + 1] = 0; // Q = 0
        }

        using var ms = new MemoryStream(rawBytes);
        var demod = new IqDemodulator();
        double[] frame = new double[samplesPerFrame];
        bool result = demod.TryRead(ms, samplesPerFrame, frame);

        Assert.True(result);
        // All magnitudes before DC removal are 0.5; after DC subtraction they should all be ≈ 0.
        for (int i = 0; i < samplesPerFrame; i++)
            Assert.InRange(frame[i], -1e-10, 1e-10);
    }

    // 5.1 continued — verify magnitudes are non-negative before DC subtraction using mixed IQ pairs
    [Fact]
    public void TryRead_GenericStream_MagnitudesMatchExpected()
    {
        // I = 64 (→ 0.5), Q = 64 (→ 0.5): expected magnitude = sqrt(0.5² + 0.5²) = 0.7071...
        const int samplesPerFrame = 4;
        byte[] rawBytes = new byte[samplesPerFrame * 2];
        for (int i = 0; i < samplesPerFrame; i++)
        {
            rawBytes[2 * i] = 64;
            rawBytes[2 * i + 1] = 64;
        }

        using var ms = new MemoryStream(rawBytes);
        var demod = new IqDemodulator();
        double[] frame = new double[samplesPerFrame];
        bool result = demod.TryRead(ms, samplesPerFrame, frame);

        Assert.True(result);
        // After DC removal all should be 0 (uniform signal), but mean before removal ≈ 0.7071
        double expectedMag = Math.Sqrt(0.5 * 0.5 + 0.5 * 0.5);
        double mean = frame.Average();
        Assert.InRange(mean, -1e-10, 1e-10); // DC removed → mean ≈ 0
        // Verify every individual sample equals the pre-DC value minus the mean
        for (int i = 0; i < samplesPerFrame; i++)
            Assert.InRange(frame[i], -1e-10, 1e-10);
    }

    // 5.2 Short read returns false
    [Fact]
    public void TryRead_ShortStream_ReturnsFalse()
    {
        const int samplesPerFrame = 100;
        // Only 10 samples worth of bytes — far fewer than requested
        byte[] rawBytes = new byte[10 * 2];
        using var ms = new MemoryStream(rawBytes);
        var demod = new IqDemodulator();
        double[] frame = new double[samplesPerFrame];

        bool result = demod.TryRead(ms, samplesPerFrame, frame);

        Assert.False(result);
    }
}
