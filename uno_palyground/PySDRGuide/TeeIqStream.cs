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
