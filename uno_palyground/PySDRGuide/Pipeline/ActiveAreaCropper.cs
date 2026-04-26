namespace uno_palyground.PySDRGuide.Pipeline;

public class ActiveAreaCropper
{
    public void CropInto(double[] signal, int samplesPerLine, double[] dest, int copyWidth)
    {
        int lines = signal.Length / samplesPerLine;
        for (int ln = 0; ln < lines; ln++)
        {
            int src = ln * samplesPerLine;
            int dst = ln * copyWidth;
            if (src + copyWidth > signal.Length) break;
            Buffer.BlockCopy(signal, src * sizeof(double), dest, dst * sizeof(double), copyWidth * sizeof(double));
        }
    }
}
