using uno_palyground.PySDRGuide.Pipeline;
using Xunit;

namespace PALDecoder.Tests;

public class ActiveAreaCropperTests
{
    // 5.11 Verify per-line copy copies exactly copyWidth samples from each line
    [Fact]
    public void CropInto_CopiesExactlyCopyWidthSamplesPerLine()
    {
        const int samplesPerLine = 100;
        const int copyWidth = 60;
        const int lines = 10;
        double[] signal = new double[samplesPerLine * lines];

        // Fill each sample with a unique value: sample index
        for (int i = 0; i < signal.Length; i++)
            signal[i] = i;

        double[] dest = new double[copyWidth * lines];

        var cropper = new ActiveAreaCropper();
        cropper.CropInto(signal, samplesPerLine, dest, copyWidth);

        // Verify: dest[ln * copyWidth + s] should equal signal[ln * samplesPerLine + s]
        for (int ln = 0; ln < lines; ln++)
            for (int s = 0; s < copyWidth; s++)
                Assert.Equal(signal[ln * samplesPerLine + s], dest[ln * copyWidth + s]);
    }

    // Also verify that samples past copyWidth in the source line are NOT copied
    [Fact]
    public void CropInto_DoesNotCopyBeyondCopyWidth()
    {
        const int samplesPerLine = 100;
        const int copyWidth = 60;
        const int lines = 5;
        double[] signal = new double[samplesPerLine * lines];
        for (int i = 0; i < signal.Length; i++)
            signal[i] = i + 1; // no zeros in signal

        double[] dest = new double[copyWidth * lines]; // initialized to 0

        var cropper = new ActiveAreaCropper();
        cropper.CropInto(signal, samplesPerLine, dest, copyWidth);

        // Destination should only have copyWidth values per line; verify line boundaries
        for (int ln = 0; ln < lines; ln++)
        {
            // First sample of each line in dest == first sample of that line in signal
            Assert.Equal(signal[ln * samplesPerLine], dest[ln * copyWidth]);
            // Last copied sample == signal[ln * samplesPerLine + copyWidth - 1]
            Assert.Equal(signal[ln * samplesPerLine + copyWidth - 1], dest[ln * copyWidth + copyWidth - 1]);
        }
    }
}
