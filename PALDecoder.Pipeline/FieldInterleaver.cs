namespace uno_palyground.PySDRGuide.Pipeline;

public class FieldInterleaver
{
    public void Interleave(byte[,,] field1, byte[,,] field2, byte[,,] dest, FieldOrder order)
    {
        int fieldLines = field1.GetLength(0);
        int width = Math.Min(field1.GetLength(1), field2.GetLength(1));
        int destHeight = dest.GetLength(0);
        int height = Math.Min(destHeight, fieldLines * 2);
        int outLine = 0;
        for (int i = 0; i < fieldLines && outLine < height - 1; i++)
        {
            for (int j = 0; j < width; j++)
            {
                if (order == FieldOrder.BottomFieldFirst)
                {
                    dest[outLine, j, 0] = field2[i, j, 0];
                    dest[outLine, j, 1] = field2[i, j, 1];
                    dest[outLine, j, 2] = field2[i, j, 2];

                    dest[outLine + 1, j, 0] = field1[i, j, 0];
                    dest[outLine + 1, j, 1] = field1[i, j, 1];
                    dest[outLine + 1, j, 2] = field1[i, j, 2];
                }
                else
                {
                    dest[outLine, j, 0] = field1[i, j, 0];
                    dest[outLine, j, 1] = field1[i, j, 1];
                    dest[outLine, j, 2] = field1[i, j, 2];

                    dest[outLine + 1, j, 0] = field2[i, j, 0];
                    dest[outLine + 1, j, 1] = field2[i, j, 1];
                    dest[outLine + 1, j, 2] = field2[i, j, 2];
                }
            }
            outLine += 2;
        }
    }
}
