namespace uno_palyground.PySDRGuide.Pipeline;

public class YuvToRgbConverter
{
    public void Convert(double[] y, double[] u, double[] v, int samplesPerLine, byte[,,] dest)
    {
        int destHeight = dest.GetLength(0);
        int destWidth = dest.GetLength(1);
        int height = Math.Min(destHeight, y.Length / samplesPerLine);
        int width = Math.Min(destWidth, samplesPerLine);

        // BT.601 YUV→RGB matrix coefficients — see PalConstants.Bt601C*
        double c_rv = PalConstants.Bt601Crv;
        double c_gu = PalConstants.Bt601Cgu;
        double c_gv = PalConstants.Bt601Cgv;
        double c_bu = PalConstants.Bt601Cbu;

        for (int row = 0; row < height; row++)
        {
            int baseIdx = row * samplesPerLine;
            for (int col = 0; col < width; col++)
            {
                int index = baseIdx + col;
                double yVal = Clamp(y[index], 0, 1);
                double uVal = Clamp(u[index], -0.5, 0.5);
                double vVal = Clamp(v[index], -0.5, 0.5);

                double r = yVal + c_rv * vVal;
                double g = yVal + c_gu * uVal + c_gv * vVal;
                double b = yVal + c_bu * uVal;

                dest[row, col, 0] = (byte)(Clamp(r, 0, 1) * 255);
                dest[row, col, 1] = (byte)(Clamp(g, 0, 1) * 255);
                dest[row, col, 2] = (byte)(Clamp(b, 0, 1) * 255);
            }
        }
    }

    private static double Clamp(double value, double min, double max) =>
        Math.Max(min, Math.Min(max, value));
}
