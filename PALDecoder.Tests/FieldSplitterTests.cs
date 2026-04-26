using uno_palyground.PySDRGuide.Pipeline;
using Xunit;

namespace PALDecoder.Tests;

public class FieldSplitterTests
{
    private const int SamplesPerLine = 100;

    // 5.6 Frame buffer filled with line-index markers; assert field boundaries
    [Fact]
    public void TrySplit_FrameBuffer_FieldBoundariesAreCorrect()
    {
        // Strategy A: field1 starts at line 0, field2 starts at line PAL_LINES_PER_FRAME/2 = 312
        int field1StartLine = 0;
        int field2StartLine = PalConstants.PAL_LINES_PER_FRAME / 2; // 312
        int fieldLines = PalConstants.PAL_VISIBLE_LINES / 2;        // 288

        // Total lines needed: field2StartLine + fieldLines = 312 + 288 = 600 lines
        int totalLines = field2StartLine + fieldLines;
        int frameLength = totalLines * SamplesPerLine;

        double[] frame = new double[frameLength];
        // Fill each sample with its line index so we can verify copies
        for (int ln = 0; ln < totalLines; ln++)
            for (int s = 0; s < SamplesPerLine; s++)
                frame[ln * SamplesPerLine + s] = ln;

        double[] field1 = new double[fieldLines * SamplesPerLine];
        double[] field2 = new double[fieldLines * SamplesPerLine];

        var splitter = new FieldSplitter();
        bool ok = splitter.TrySplit(frame, SamplesPerLine, field1, field2);

        Assert.True(ok);
        // field1[0..SamplesPerLine-1] should come from line field1StartLine = 0
        Assert.Equal(field1StartLine, (int)field1[0]);
        // field2[0..SamplesPerLine-1] should come from line field2StartLine = 312
        Assert.Equal(field2StartLine, (int)field2[0]);
    }

    // 5.7 Insufficient lines returns false
    [Fact]
    public void TrySplit_InsufficientLines_ReturnsFalse()
    {
        // Buffer too small (only 100 lines) — need 312 + 288 = 600
        double[] frame = new double[100 * SamplesPerLine];
        double[] field1 = new double[288 * SamplesPerLine];
        double[] field2 = new double[288 * SamplesPerLine];

        var splitter = new FieldSplitter();
        bool ok = splitter.TrySplit(frame, SamplesPerLine, field1, field2);

        Assert.False(ok);
    }
}
