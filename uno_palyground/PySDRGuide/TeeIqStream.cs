using System;
using System.Collections.Concurrent;
using System.IO;
using System.Threading;

namespace HackRF.Namespace
{
    // Stream that supplies a continuous sequence of captured IQ bytes added by producer via AddChunk.
    // Each chunk is raw interleaved int8 I,Q samples. Consumer (PALDecoder) reads sequentially.
    internal sealed class TeeIqStream : Stream
    {
        private readonly ConcurrentQueue<byte[]> _queue = new();
        private byte[]? _current;
        private int _offset;
        private bool _completed;
        private readonly ManualResetEventSlim _dataEvent = new(false);
        private readonly object _sync = new();

        public void AddChunk(byte[] buffer, int validLength)
        {
            if (validLength <= 0 || buffer.Length == 0 || _completed) return;
            if (validLength != buffer.Length)
            {
                var trimmed = new byte[validLength];
                Buffer.BlockCopy(buffer, 0, trimmed, 0, validLength);
                _queue.Enqueue(trimmed);
            }
            else
            {
                _queue.Enqueue(buffer);
            }
            _dataEvent.Set();
        }

        public void Complete()
        {
            _completed = true;
            _dataEvent.Set();
        }

        // Discard queued data older than the last frameBytes worth; keeps only most recent segment to reduce latency.
        // frameBytes: number of raw interleaved I/Q bytes representing one PAL frame (samplesPerFrame * 2).
        public int DrainToLatestFrame(int frameBytes)
        {
            if (frameBytes <= 0) return 0;
            lock (_sync)
            {
                if (_queue.IsEmpty) return 0;
                // Aggregate all queued buffers
                int total = 0;
                foreach (var b in _queue) total += b.Length;
                if (total <= frameBytes) return total; // nothing to drain; we keep all
                // Build combined tail slice of size frameBytes
                byte[] combined = new byte[frameBytes];
                int copyStart = total - frameBytes; // offset into concatenated stream
                int copied = 0;
                int traversed = 0;
                foreach (var b in _queue)
                {
                    int nextTraversed = traversed + b.Length;
                    if (nextTraversed > copyStart && copied < frameBytes)
                    {
                        int srcOffset = Math.Max(0, copyStart - traversed);
                        int available = b.Length - srcOffset;
                        int toCopy = Math.Min(available, frameBytes - copied);
                        Buffer.BlockCopy(b, srcOffset, combined, copied, toCopy);
                        copied += toCopy;
                    }
                    traversed = nextTraversed;
                }
                // Replace queue with single latest frame buffer
                while (_queue.TryDequeue(out _)) { }
                _queue.Enqueue(combined);
                _current = null; _offset = 0; // reset current cursor so reader starts at new buffer
                return frameBytes;
            }
        }

        public override int Read(byte[] dest, int offset, int count)
        {
            if (dest == null) throw new ArgumentNullException(nameof(dest));
            if (offset < 0 || count < 0 || (offset + count) > dest.Length) throw new ArgumentOutOfRangeException();
            int written = 0;
            while (written < count)
            {
                if (_current == null || _offset >= _current.Length)
                {
                    if (!_queue.TryDequeue(out _current))
                    {
                        if (_completed) break;
                        _dataEvent.Wait(50);
                        _dataEvent.Reset();
                        continue;
                    }
                    _offset = 0;
                }
                int toCopy = Math.Min(count - written, _current.Length - _offset);
                Buffer.BlockCopy(_current, _offset, dest, offset + written, toCopy);
                _offset += toCopy;
                written += toCopy;
            }
            return written;
        }

        public override bool CanRead => true;
        public override bool CanSeek => false;
        public override bool CanWrite => false;
        public override long Length => throw new NotSupportedException();
        public override long Position { get => throw new NotSupportedException(); set => throw new NotSupportedException(); }
        public override void Flush() { }
        public override long Seek(long offset, SeekOrigin origin) => throw new NotSupportedException();
        public override void SetLength(long value) => throw new NotSupportedException();
        public override void Write(byte[] buffer, int offset, int count) => throw new NotSupportedException();
        protected override void Dispose(bool disposing)
        {
            if (disposing)
            {
                Complete();
                _dataEvent.Dispose();
            }
            base.Dispose(disposing);
        }
    }
}
