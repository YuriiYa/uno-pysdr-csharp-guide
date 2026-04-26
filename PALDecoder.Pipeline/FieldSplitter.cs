namespace uno_palyground.PySDRGuide.Pipeline;

public class FieldSplitter
{
    public bool TrySplit(double[] frameData, int samplesPerLine, double[] field1, double[] field2)
    {
        // TODO: check the correctness: Strategy A assumes frameData starts at Field-1 active line 0 with no VBI; verify this holds for all stream sources
        int field1StartLine = 0;
        int field2StartLine = PalConstants.PAL_LINES_PER_FRAME / 2;  // 312
        int fieldLinesVis = PalConstants.PAL_VISIBLE_LINES / 2;      // 288

        int availableLines = frameData.Length / samplesPerLine;
        if (availableLines < (field2StartLine + fieldLinesVis)) return false;

        int srcField1OffsetSamples = field1StartLine * samplesPerLine;
        int srcField2OffsetSamples = field2StartLine * samplesPerLine;
        int copySamples = fieldLinesVis * samplesPerLine;
        int bytesToCopy = copySamples * sizeof(double);

        if (srcField1OffsetSamples + copySamples <= frameData.Length)
            Buffer.BlockCopy(frameData, srcField1OffsetSamples * sizeof(double), field1, 0, bytesToCopy);
        if (srcField2OffsetSamples + copySamples <= frameData.Length)
            Buffer.BlockCopy(frameData, srcField2OffsetSamples * sizeof(double), field2, 0, bytesToCopy);
        return true;
    }
}
