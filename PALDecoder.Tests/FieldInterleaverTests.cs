using uno_palyground.PySDRGuide.Pipeline;
using Xunit;

namespace PALDecoder.Tests;

public class FieldInterleaverTests
{
    private const int FieldLines = 4;
    private const int Width = 8;

    private static byte[,,] MakeField(byte value)
    {
        byte[,,] f = new byte[FieldLines, Width, 3];
        for (int r = 0; r < FieldLines; r++)
            for (int c = 0; c < Width; c++)
                for (int ch = 0; ch < 3; ch++)
                    f[r, c, ch] = value;
        return f;
    }

    // 5.13 BottomFieldFirst → row 0 = field2, row 1 = field1
    [Fact]
    public void Interleave_BottomFieldFirst_EvenRowsAreField2OddRowsAreField1()
    {
        byte[,,] field1 = MakeField(1);
        byte[,,] field2 = MakeField(2);
        byte[,,] dest = new byte[FieldLines * 2, Width, 3];

        var interleaver = new FieldInterleaver();
        interleaver.Interleave(field1, field2, dest, FieldOrder.BottomFieldFirst);

        // Even rows (0, 2, 4, ...) come from field2 (value = 2)
        for (int r = 0; r < FieldLines * 2; r += 2)
            Assert.Equal(2, dest[r, 0, 0]);

        // Odd rows (1, 3, 5, ...) come from field1 (value = 1)
        for (int r = 1; r < FieldLines * 2; r += 2)
            Assert.Equal(1, dest[r, 0, 0]);
    }

    // 5.14 TopFieldFirst → row 0 = field1, row 1 = field2
    [Fact]
    public void Interleave_TopFieldFirst_EvenRowsAreField1OddRowsAreField2()
    {
        byte[,,] field1 = MakeField(1);
        byte[,,] field2 = MakeField(2);
        byte[,,] dest = new byte[FieldLines * 2, Width, 3];

        var interleaver = new FieldInterleaver();
        interleaver.Interleave(field1, field2, dest, FieldOrder.TopFieldFirst);

        // Even rows come from field1 (value = 1)
        for (int r = 0; r < FieldLines * 2; r += 2)
            Assert.Equal(1, dest[r, 0, 0]);

        // Odd rows come from field2 (value = 2)
        for (int r = 1; r < FieldLines * 2; r += 2)
            Assert.Equal(2, dest[r, 0, 0]);
    }
}
