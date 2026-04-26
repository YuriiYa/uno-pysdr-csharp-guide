using uno_palyground.PySDRGuide.Pipeline;
using Xunit;

namespace PALDecoder.Tests;

public class YuvToRgbConverterTests
{
    // 5.12 White (Y=1, U=0, V=0) → (255, 255, 255)
    [Fact]
    public void Convert_White_AllChannels255()
    {
        const int width = 4;
        const int height = 2;
        double[] y = new double[width * height];
        double[] u = new double[width * height];
        double[] v = new double[width * height];
        Array.Fill(y, 1.0); // white
        // u and v remain 0.0

        byte[,,] dest = new byte[height, width, 3];
        var converter = new YuvToRgbConverter();
        converter.Convert(y, u, v, width, dest);

        for (int row = 0; row < height; row++)
            for (int col = 0; col < width; col++)
            {
                Assert.Equal(255, dest[row, col, 0]); // R
                Assert.Equal(255, dest[row, col, 1]); // G
                Assert.Equal(255, dest[row, col, 2]); // B
            }
    }

    // 5.12 Black (Y=0, U=0, V=0) → (0, 0, 0)
    [Fact]
    public void Convert_Black_AllChannels0()
    {
        const int width = 4;
        const int height = 2;
        double[] y = new double[width * height]; // all 0.0
        double[] u = new double[width * height];
        double[] v = new double[width * height];

        byte[,,] dest = new byte[height, width, 3];
        var converter = new YuvToRgbConverter();
        converter.Convert(y, u, v, width, dest);

        for (int row = 0; row < height; row++)
            for (int col = 0; col < width; col++)
            {
                Assert.Equal(0, dest[row, col, 0]);
                Assert.Equal(0, dest[row, col, 1]);
                Assert.Equal(0, dest[row, col, 2]);
            }
    }

    // 5.12 Out-of-range Y=1.5 → clamped to 255
    [Fact]
    public void Convert_OverrangeY_ClampedTo255()
    {
        const int width = 2;
        const int height = 1;
        double[] y = new double[width * height];
        double[] u = new double[width * height];
        double[] v = new double[width * height];
        Array.Fill(y, 1.5); // above maximum — should clamp

        byte[,,] dest = new byte[height, width, 3];
        var converter = new YuvToRgbConverter();
        converter.Convert(y, u, v, width, dest);

        for (int col = 0; col < width; col++)
        {
            Assert.Equal(255, dest[0, col, 0]);
            Assert.Equal(255, dest[0, col, 1]);
            Assert.Equal(255, dest[0, col, 2]);
        }
    }
}
