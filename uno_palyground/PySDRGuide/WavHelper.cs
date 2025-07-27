using System.Diagnostics;
using System.Runtime.InteropServices;

public class WavHelper
{
    public static TType? ConvertByteArraytoType<TType>(Int16[] int16Header)
    {
        byte[] byteHeader = new byte[int16Header.Length * sizeof(Int16)];
        Buffer.BlockCopy(int16Header, 0, byteHeader, 0, byteHeader.Length);
        return ConvertByteArraytoType<TType>(byteHeader);
    }

    public static TType? ConvertByteArraytoType<TType>(byte[] header)
    {
        IntPtr ptPoit = Marshal.AllocHGlobal(header.Length);
        Marshal.Copy(header, 0, ptPoit, header.Length);
        TType? structureFromByte = Marshal.PtrToStructure<TType>(ptPoit);
        Marshal.FreeHGlobal(ptPoit);
        return structureFromByte;
    }

    public byte[] ReadWaveHeader(string flePath)
    {
        byte[] buffer = new byte[44]; // header consist of 44 bytes
        try
        {
            using (FileStream fs = new FileStream(flePath, FileMode.Open, FileAccess.Read))
            {
                var bytes_read = fs.Read(buffer, 0, buffer.Length);
                fs.Close();

                if (bytes_read != buffer.Length)
                {
                    // Couldn't read 512 bytes
                }
            }
        }
        catch (System.UnauthorizedAccessException ex)
        {
            Debug.Print(ex.Message);
        }
        return buffer;
    }
}
