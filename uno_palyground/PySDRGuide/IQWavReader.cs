using System.Numerics;
using System.IO;
using System.Runtime.InteropServices; // For MemoryMarshal.AsBytes

namespace uno_palyground.PySDRGuide;

/// <summary>
/// Utility to read an SDRSharp-style IQ WAV file (44-byte PCM header + interleaved 16-bit IQ samples)
/// without relying on Numpy.NET. Loads the entire file into memory similarly to the previous numpy code.
/// </summary>
public static class IQWavReader
{
    // TODO:  Returns spans into a reusable ring buffer (no Complex[] allocation), or
    // Writes into a caller-provided Span<Complex> target,
    // I can add those too.
    /// <summary>
    /// Reads all interleaved IQ samples (Int16 I,Q pairs) from disk.
    /// Scales I/Q to [-1, +1] range (double) and returns as Complex[].
    /// </summary>
    /// <param name="fs">FileStream for IQ WAV file</param>
    /// <param name="length">Number of bytes to read</param>
    /// <returns>Tuple of parsed WavHeader and Complex[] IQ samples</returns>
    public static Complex[] ReadIQWavNext(FileStream fs, long length)
    {
        // Read only a portion of the remaining stream (length bytes), converting directly to Complex on-the-fly.
        // No intermediate short[] or large byte[] heap allocations; only the final Complex[] is allocated.

        long remaining = fs.Length - fs.Position;
        if (remaining <= 0 || length <= 0)
            return Array.Empty<Complex>();

        long targetBytes = Math.Min(remaining, length);
        // We need complete IQ pairs: 4 bytes (I int16 + Q int16) per Complex sample.
        long usableBytes = targetBytes - (targetBytes % 4); // floor to multiple of 4
        if (usableBytes <= 0)
            return Array.Empty<Complex>();

        long totalPairsLong = usableBytes / 4;
        if (totalPairsLong > int.MaxValue)
            throw new NotSupportedException("Requested segment too large to materialize into a single Complex[]");
        int totalPairs = (int)totalPairsLong;

        Complex[] result = new Complex[totalPairs];
        double scale = 1.0 / short.MaxValue;

        const int ChunkPairs = 4096; // tuneable; 4096 pairs = 16 KB of raw bytes
        Span<byte> buffer = stackalloc byte[ChunkPairs * 4];

        int written = 0;
        while (written < totalPairs)
        {
            int pairsNeeded = totalPairs - written;
            if (pairsNeeded > ChunkPairs) pairsNeeded = ChunkPairs;
            int bytesNeeded = pairsNeeded * 4;

            int readTotal = 0;
            while (readTotal < bytesNeeded)
            {
                int n = fs.Read(buffer.Slice(readTotal, bytesNeeded - readTotal));
                if (n == 0) break; // EOF encountered mid-chunk
                readTotal += n;
            }

            if (readTotal == 0)
            {
                // EOF before reading any more data; shrink if partial earlier
                if (written < totalPairs)
                {
                    if (written == 0) return Array.Empty<Complex>();
                    Array.Resize(ref result, written);
                }
                break;
            }

            int usable = readTotal - (readTotal % 4); // ensure whole pairs
            int pairsRead = usable / 4;
            int b = 0;
            for (int p = 0; p < pairsRead; p++)
            {
                short iVal = (short)(buffer[b] | (buffer[b + 1] << 8));
                short qVal = (short)(buffer[b + 2] | (buffer[b + 3] << 8));
                b += 4;
                result[written + p] = new Complex(iVal * scale, qVal * scale);
            }
            written += pairsRead;

            if (pairsRead < pairsNeeded) // hit EOF or partial chunk
            {
                if (written < totalPairs)
                    Array.Resize(ref result, written);
                break;
            }
        }

        return result;
    }

    public static void SkipBytes(FileStream fs, long length)
    {
        byte[] buffer = new byte[length];
        long totalSkipped = 0;

        while (totalSkipped < length)
        {
            int n = fs.Read(buffer, (int)totalSkipped, (int)(length - totalSkipped));
            if (n == 0) break; // EOF encountered mid-chunk
            totalSkipped += n;
        }
    }
}
