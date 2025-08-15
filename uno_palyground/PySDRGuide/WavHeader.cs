using System.Runtime.InteropServices;

[StructLayout(LayoutKind.Sequential, Pack = 1)]
    public struct WavHeader
    {
        //	Marks the file as a riff file. Characters are each 1 byte long.
        [MarshalAs(UnmanagedType.ByValTStr, SizeConst = 4)]
        public string FileMark;
        //Size of the overall file - 8 bytes, in bytes (32-bit integer). Typically, you’d fill this in after creation.
        public Int32 Size;
        // File Type Header. For our purposes, it always equals “WAVE”.
        [MarshalAs(UnmanagedType.ByValTStr, SizeConst = 4)]
        public string FileFormat;
        // Format chunk marker. Includes trailing null
        [MarshalAs(UnmanagedType.ByValTStr, SizeConst = 4)]
        public string ChunkMarker;
        // Length of format data as listed above
        public Int32 ChankLength;
        //	Type of format (1 is PCM) - 2 byte integer
        public Int16 AudioFormat;
        // Number of Channels - 2 byte integer
        public Int16 Channels;
        // Sample Rate - 32 bit integer. Common values are 44100 (CD), 48000 (DAT). Sample Rate = Number of Samples per second, or Hertz.
        public Int32 SampleRate;
        // (Sample Rate * BitsPerSample * Channels) / 8.
        public Int32 AverageBytePerSecond;
        // (BitsPerSample * Channels) / 8.1 - 8 bit mono2 - 8 bit stereo/16 bit mono4 - 16 bit stereo
        public Int16 BlockAlign;
        // Bits per sample
        public Int16 AudioBitPerSample;
        //“data” chunk header. Marks the beginning of the data section.
        public Int32 Data;
        // Size of the data section.
        public Int32 DataFileSize;
    }
