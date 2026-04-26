# uno-pysdr-csharp-guide

C# reimplementation of [pysdr.org](https://pysdr.org) SDR examples using the [Uno Platform](https://platform.uno) desktop framework, with standalone Python scripts for comparison.

## Solution structure

```
uno_palyground.sln
├── uno_palyground/              # Uno Platform desktop app (net9.0-desktop, Uno.Sdk)
│   ├── Presentation/            # MVVM: *Page.xaml + *Page.xaml.cs + *ViewModel.cs triads
│   ├── PySDRGuide/              # C# reimplementations of pysdr.org chapters
│   └── Services/Endpoints/      # Dependency-injected services
├── PALDecoder.Pipeline/         # .NET 9 class library – PAL video decode pipeline stages
│   ├── IqDemodulator.cs         # IQ → real baseband demodulation
│   ├── FrameSynchronizer.cs     # V-sync detection (Goertzel)
│   ├── HorizontalAligner.cs     # H-sync alignment via LPF + median
│   ├── FieldSplitter.cs         # Interlaced field separation
│   ├── LumaChromaSeparator.cs   # LPF/BPF luma–chroma separation
│   ├── ChromaDecoder.cs         # PAL colour PLL + U/V demodulation
│   ├── ActiveAreaCropper.cs     # Active video line extraction
│   ├── YuvToRgbConverter.cs     # BT.601 YUV → RGB conversion
│   ├── FieldInterleaver.cs      # Interlaced field recombination
│   ├── ColorPll.cs              # PAL colour sub-carrier PLL
│   ├── IQWavReader.cs           # WAV file reader for IQ samples
│   ├── PalConstants.cs          # PAL timing and frequency constants
│   └── PalTypes.cs              # TvSystem, SystemProfile, FieldOrder types
├── PALDecoder.Tests/            # xUnit test project (net9.0) – 18 unit tests
│   └── *Tests.cs                # One file per pipeline stage
├── NetHackRF/                   # .NET 9 class library – P/Invoke wrapper for hackrf.dll
└── python/                      # Standalone Python scripts (separate venv)
```

## Build & run

```sh
# Build the entire solution
dotnet build uno_palyground.sln

# Run the desktop app
dotnet run --project uno_palyground/uno_palyground.csproj

# Run unit tests
dotnet test PALDecoder.Tests/PALDecoder.Tests.csproj

# Build a single library
dotnet build PALDecoder.Pipeline/PALDecoder.Pipeline.csproj
```

## Python environment

```sh
# One-time setup
python -m venv venv
./venv/Scripts/activate          # Windows
# source venv/bin/activate       # Linux/macOS
pip install git+https://github.com/GvozdevLeonid/python_hackrf.git
pip install scipy matplotlib

# Run a script
py python/HackRFReadSamples.py

# Deactivate when done
deactivate
```

## Key libraries

| Library | Purpose |
|---------|---------|
| [Uno Platform](https://platform.uno) | Cross-platform WinUI desktop UI |
| [FftSharp](https://github.com/swharden/FftSharp) | FFT utilities |
| [ScottPlot](https://scottplot.net) | Signal plots |
| [Spectrogram](https://github.com/swharden/Spectrogram) | Waterfall spectrogram |
| [Numpy (.NET)](https://github.com/SciSharp/Numpy.NET) | Array math |
| [CommunityToolkit.Mvvm](https://learn.microsoft.com/dotnet/communitytoolkit/mvvm/) | MVVM source generators |
| [HackRF](https://github.com/greatscottgadgets/hackrf) | SDR hardware (via P/Invoke) |

## HackRF IQ format

Raw samples from HackRF are interleaved `int8` bytes: `[I₀, Q₀, I₁, Q₁, …]`.  
Convert to `Complex`:

```csharp
new Complex((sbyte)buf[2*i], (sbyte)buf[2*i+1]) / 128.0
```

## Package management

Package versions are centrally managed in `Directory.Packages.props` (CPM).  
Add new packages there first, then reference without a version in `.csproj`.