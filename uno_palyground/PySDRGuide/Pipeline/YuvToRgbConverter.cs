namespace uno_palyground.PySDRGuide.Pipeline;

public class YuvToRgbConverter
{
    public void Convert(double[] y, double[] u, double[] v, int samplesPerLine, byte[,,] dest)
    {
        int destHeight = dest.GetLength(0);
        int destWidth = dest.GetLength(1);
        int height = Math.Min(destHeight, y.Length / samplesPerLine);
        int width = Math.Min(destWidth, samplesPerLine);

        // BT.601 YUV→RGB matrix coefficients (full-range, R′G′B′ = Y′ + matrix × [U, V])
        // c_rv: V → R contribution  (R = Y + 1.402 * V)
        // c_gu: U → G contribution  (G = Y − 0.344 * U − 0.714 * V)
        // c_gv: V → G contribution
        // c_bu: U → B contribution  (B = Y + 1.772 * U)
        const double c_rv = 1.402;
        const double c_gu = -0.344;
        const double c_gv = -0.714;
        const double c_bu = 1.772;

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
