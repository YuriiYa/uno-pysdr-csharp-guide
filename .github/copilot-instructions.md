# Copilot Instructions

## Project Overview

This repo implements [pysdr.org](https://pysdr.org) SDR examples in C# with an [Uno Platform](https://platform.uno) desktop UI, plus standalone Python scripts for comparison. The goal is a C# counterpart to pysdr.org's Python-based guide.

## Solution Structure

```
uno_palyground.sln
├── uno_palyground/          # Uno Platform desktop app (net9.0-desktop, Uno.Sdk 6.0.96)
│   ├── Presentation/        # MVVM: *Page.xaml + *Page.xaml.cs + *ViewModel.cs triads
│   ├── PySDRGuide/          # C# reimplementations of pysdr.org chapters
│   └── Services/Endpoints/  # Dependency-injected services
├── NetHackRF/               # .NET 9 class library – P/Invoke wrapper for hackrf.dll
└── python/                  # Standalone Python scripts (separate venv)
```

## Build & Run

```sh
# Build the solution
dotnet build uno_palyground.sln

# Run the desktop app
dotnet run --project uno_palyground/uno_palyground.csproj

# Build a single project
dotnet build NetHackRF/NetHackrf.csproj
```

There are no automated tests in this repository.

## Python Environment

```sh
# One-time setup
python -m venv venv
./venv/Scripts/activate          # Windows
pip install git+https://github.com/GvozdevLeonid/python_hackrf.git
pip install scipy matplotlib

# Run a script
py python/HackRFReadSamples.py
```

## Key Conventions

### MVVM / Navigation

- Each screen is a triad: `*Page.xaml`, `*Page.xaml.cs` (code-behind), `*ViewModel.cs`.
- ViewModels use `[ObservableProperty]` / `[RelayCommand]` from **CommunityToolkit.Mvvm** (partial classes).
- Navigation is registered in `App.xaml.cs` via `ViewMap` / `RouteMap` using Uno Toolkit navigation. Adding a new page requires entries in both `views.Register(...)` and `routes.Register(...)`.
- `GlobalUsings.cs` declares solution-wide `global using` statements; avoid repeating those imports in individual files.

### PySDRGuide Implementations

- Each file in `uno_palyground/PySDRGuide/` maps to one pysdr.org chapter. The originating URL is noted in a comment at the top of each file (e.g., `// https://pysdr.org/content/frequency_domain.html`).
- Signal processing uses **Numpy** (C# NuGet port) for array math, **FftSharp** for FFT utilities, **ScottPlot** for plotting, and **Spectrogram** for waterfall visuals.
- IQ samples from HackRF are raw **interleaved int8** bytes: `[I₀, Q₀, I₁, Q₁, …]`. Convert to `Complex` via `new Complex((sbyte)buf[2*i], (sbyte)buf[2*i+1]) / 128.0`.

### NetHackRF Library

- Wraps `hackrf.dll` (native binary, must be present alongside the executable) via `libhackrf.cs` P/Invoke declarations.
- Uses `unsafe` code — **Debug** build targets `x64`, **Release** targets `x86`.
- `NetHackrf` is `IDisposable`; always `Dispose()` or `using`-scope it to stop the radio and free resources. Never leave a stream open when switching modes.
- HackRF device properties use setters only (read back from hardware is not supported for most settings): `CarrierFrequencyMHz`, `SampleFrequencyMHz`, `LNAGainDb` (0–40 dB, 8 dB steps), `VGAGainDb` (0–62 dB, 2 dB steps), `FilterBandwidthMHz`.
- `TeeIqStream` (in `PySDRGuide/`) multiplexes the single hardware RX stream to multiple consumers (e.g., spectrogram + PAL decoder) without duplicating device handles.

### Package Management

- Package **versions** are centrally managed in `Directory.Packages.props` (CPM). Add new packages there first, then reference them without a version in `.csproj`.
- Uno SDK version is updated only in `global.json` (`msbuild-sdks` → `Uno.Sdk`), not in project files.
- `UnoFeatures` in the `.csproj` controls implicit Uno package references (Material, DSP, Navigation, etc.).

### SDR Abbreviations

A full abbreviation reference lives at the top of `uno_palyground/PySDRGuide/PALDecoder.cs` (ADC, LNA, VGA, IQ, FFT, FIR, etc.).

## OpenSpec Change Workflow

Design changes are tracked under `openspec/`. Use the skill shortcuts:
- **Propose** a change: `/opsx-propose` (`.github/prompts/opsx-propose.prompt.md`)
- **Explore** an idea: `/opsx-explore`
- **Apply** a change: `/opsx-apply`
- **Archive** when done: `/opsx-archive`
