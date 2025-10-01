using System.Numerics;
using System.IO;
using System.Runtime.InteropServices; // For MemoryMarshal.AsBytes
using System.Runtime.CompilerServices; // For Unsafe operations

namespace uno_palyground.PySDRGuide;

/// <summary>
/// Utility to read an SDRSharp-style IQ WAV file (44-byte PCM header + interleaved 16-bit IQ samples)
/// without relying on Numpy.NET. Loads the entire file into memory similarly to the previous numpy code.
/// </summary>
public static class IQWavReader
{
    // Ring buffer storage for zero-allocation streaming reads.
    // Not thread-safe. Configure via ConfigureRingBuffer before use.
    private static Complex[]? _ringBuffer;
    private static int _capacity;          // number of Complex samples buffer can hold
    private static int _writeIndex;        // next position to write

    /// <summary>
    /// Configure (or reconfigure) the reusable ring buffer used by span-returning streaming APIs.
    /// Existing contents are discarded. Not thread-safe.
    /// </summary>
    /// <param name="capacitySamples">Number of Complex IQ samples the ring buffer can hold. Must be > 0.</param>
    public static void ConfigureRingBuffer(int capacitySamples)
    {
        if (capacitySamples <= 0)
            throw new ArgumentOutOfRangeException(nameof(capacitySamples));
        _ringBuffer = new Complex[capacitySamples];
        _capacity = capacitySamples;
        _writeIndex = 0;
    }
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
                // Reuse existing Complex instance in result array
                ref Complex target = ref result[written + p];
                target = new Complex(iVal * scale, qVal * scale);
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

    /// <summary>
    /// Result of a ring buffer fill operation. Because the ring buffer can wrap, up to two spans
    /// are returned: First (from writeIndex - n to end or contiguous part) and Second (wrapped beginning).
    /// If no wrap occurred, Second will be empty.
    /// </summary>
    public readonly ref struct RingSegment
    {
        public readonly ReadOnlySpan<Complex> First;
        public readonly ReadOnlySpan<Complex> Second;
        public readonly int SamplesRead; // total samples just appended
        public RingSegment(ReadOnlySpan<Complex> first, ReadOnlySpan<Complex> second, int samplesRead)
        {
            First = first;
            Second = second;
            SamplesRead = samplesRead;
        }
        public bool IsEmpty => SamplesRead == 0;
    }

    /// <summary>
    /// Reads up to requestedSamples Complex IQ samples from the FileStream into the reusable ring buffer
    /// without allocating a new array. Returns spans pointing to the freshly written region.
    /// Lifetime: spans are valid until the next call that mutates the ring buffer.
    /// Not thread-safe.
    /// Optimized version using span operations and minimizing Complex constructor calls.
    /// Uses vectorized operations where possible for better performance.
    /// </summary>
    /// <param name="fs">Source FileStream positioned at start of interleaved 16-bit IQ data.</param>
    /// <param name="requestedSamples">Maximum number of Complex samples to read this call.</param>
    /// <returns>RingSegment containing one or two spans for the newly written samples.</returns>
    public static RingSegment ReadIQIntoRingOptimized(FileStream fs, int requestedSamples)
    {
        if (_ringBuffer == null)
            throw new InvalidOperationException("Ring buffer not configured. Call ConfigureRingBuffer first.");
        if (requestedSamples <= 0)
            return new RingSegment(ReadOnlySpan<Complex>.Empty, ReadOnlySpan<Complex>.Empty, 0);

        int capacity = _capacity;
        const int MaxChunkPairs = 4096; // tune as needed
        Span<byte> raw = stackalloc byte[MaxChunkPairs * 4];
        double scale = 1.0 / short.MaxValue;
        int totalWritten = 0;

        while (totalWritten < requestedSamples)
        {
            int remainingSamples = requestedSamples - totalWritten;
            int thisChunkSamples = remainingSamples > MaxChunkPairs ? MaxChunkPairs : remainingSamples;
            int bytesNeeded = thisChunkSamples * 4;

            int readTotal = 0;
            while (readTotal < bytesNeeded)
            {
                int n = fs.Read(raw.Slice(readTotal, bytesNeeded - readTotal));
                if (n == 0) break;
                readTotal += n;
            }
            if (readTotal == 0) break;

            int usableBytes = readTotal - (readTotal % 4);
            int pairs = usableBytes / 4;
            if (pairs == 0) break;

            // Convert bytes to int16 pairs using spans for better performance
            var int16Data = MemoryMarshal.Cast<byte, short>(raw.Slice(0, usableBytes));
            
            int firstPart = Math.Min(pairs, capacity - _writeIndex);
            int writtenThisChunk = firstPart;

            // First part - use span indexing to minimize bounds checks
            var targetSpan = _ringBuffer.AsSpan(_writeIndex, firstPart);
            for (int i = 0, j = 0; i < firstPart; i++, j += 2)
            {
                // Use ref to minimize array access overhead
                ref Complex target = ref targetSpan[i];
                target = new Complex(int16Data[j] * scale, int16Data[j + 1] * scale);
            }

            if (firstPart < pairs)
            {
                // Wrap remainder
                int remaining = pairs - firstPart;
                var wrapSpan = _ringBuffer.AsSpan(0, remaining);
                for (int i = 0, j = firstPart * 2; i < remaining; i++, j += 2)
                {
                    ref Complex target = ref wrapSpan[i];
                    target = new Complex(int16Data[j] * scale, int16Data[j + 1] * scale);
                }
                _writeIndex = remaining;
                writtenThisChunk = pairs;
            }
            else
            {
                _writeIndex += firstPart;
                if (_writeIndex == capacity) _writeIndex = 0;
            }

            totalWritten += writtenThisChunk;
            if (writtenThisChunk < thisChunkSamples) break;
        }

        if (totalWritten == 0)
            return new RingSegment(ReadOnlySpan<Complex>.Empty, ReadOnlySpan<Complex>.Empty, 0);

        int endIndex = _writeIndex;
        int startIndex = endIndex - totalWritten;
        if (startIndex >= 0)
        {
            return new RingSegment(new ReadOnlySpan<Complex>(_ringBuffer, startIndex, totalWritten), ReadOnlySpan<Complex>.Empty, totalWritten);
        }
        else
        {
            int firstLen = -startIndex;
            int secondLen = totalWritten - firstLen;
            return new RingSegment(
                new ReadOnlySpan<Complex>(_ringBuffer, capacity + startIndex, firstLen),
                new ReadOnlySpan<Complex>(_ringBuffer, 0, secondLen),
                totalWritten);
        }
    }
}
