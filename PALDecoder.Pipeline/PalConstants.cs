namespace uno_palyground.PySDRGuide.Pipeline;

/// <summary>
/// Central repository of PAL signal timing and filter constants.
/// All magic numbers used across pipeline stages should be defined here.
/// </summary>
public static class PalConstants
{
    // ── Line timing (PAL ITU-R BT.470) ──────────────────────────────────────
    /// <summary>Front porch: 1.5 µs silent period before H-sync pulse.</summary>
    public const double FrontPorchSeconds = 1.5e-6;

    /// <summary>Horizontal sync pulse duration: the −0.3 V "blacker-than-black" pulse used for line lock.</summary>
    public const double HSyncSeconds = 4.7e-6;

    /// <summary>Breezeway: 0.6 µs gap at 0 V between H-sync and color burst to isolate them.</summary>
    public const double BreezewaySeconds = 0.6e-6;

    /// <summary>Color burst duration: ~10 cycles of 4.43 MHz subcarrier used as phase reference.</summary>
    public const double BurstDurationSeconds = 2.25e-6;

    /// <summary>Back porch total: breezeway + burst + residual blanking (5.8 µs nominal).</summary>
    public const double BackPorchSeconds = 5.8e-6;

    /// <summary>
    /// Nominal offset (in lines) from the start of the first broad V-sync pulse
    /// to the first active video line. Used as a fallback when burst detection fails.
    /// </summary>
    public const double VbiLinesToActiveVideo = 22.5;

    // ── VBI / sync detection ─────────────────────────────────────────────────
    /// <summary>Minimum run-length (seconds) for a "broad" V-sync equalising/serration pulse (≥20 µs).</summary>
    public const double BroadPulseMinSeconds = 20e-6;

    /// <summary>Search window duration (seconds) used when scanning for H-sync pulses per line.</summary>
    public const double HSyncSearchSeconds = 30e-6;

    // ── Robust statistics ────────────────────────────────────────────────────
    /// <summary>
    /// MAD (Median Absolute Deviation) near-zero guard. When MAD falls below this value
    /// the signal is essentially flat; switch to a range-based fallback denominator.
    /// </summary>
    public const double MadNearZeroThreshold = 1e-12;

    /// <summary>
    /// Burst-power MAD near-zero guard (used in RefineToFirstActiveLine).
    /// Chosen tighter than MadNearZeroThreshold because burst energies are squared.
    /// </summary>
    public const double BurstMadNearZeroThreshold = 1e-20;

    /// <summary>
    /// Fallback range divisor when the signal range is used instead of MAD.
    /// Dividing by 50 gives a ~2 % of range per step, empirically reasonable.
    /// </summary>
    public const double MadFallbackRangeDivisor = 50.0;

    /// <summary>
    /// V-sync threshold: number of MADs below the median that a sample must reach
    /// to be classified as part of a sync pulse (empirically tuned).
    /// </summary>
    public const double VSyncMadMultiplier = 3.0;

    /// <summary>
    /// H-sync threshold: number of MADs below the median used when refining the
    /// sync-pulse run boundary during horizontal alignment (empirically tuned).
    /// </summary>
    public const double HSyncMadMultiplier = 2.5;

    /// <summary>
    /// Number of lines to examine after the V-sync broad pulse when scanning
    /// for lines with a valid color burst. Covers the VBI (≈22 lines) plus margin.
    /// </summary>
    public const int BurstSearchLinesAfterVSync = 60;

    // ── Light low-pass FIR (LightLowPass helper) ────────────────────────────
    /// <summary>Tap count for the lightweight pre-filter applied before sync detection.</summary>
    public const int LightLpfTaps = 51;

    /// <summary>Hamming window coefficient α (raised-cosine a₀ = 0.54).</summary>
    public const double HammingAlpha = 0.54;

    /// <summary>Hamming window coefficient β (raised-cosine a₁ = 0.46).</summary>
    public const double HammingBeta = 0.46;

    // ── Luma/chroma filter design ────────────────────────────────────────────
    /// <summary>
    /// Safety margin applied to the luma LPF cutoff as a fraction of the sample rate.
    /// Ensures the luma filter never aliases past the Nyquist boundary (0.5 × Fs).
    /// Using 0.45 (not 0.5) leaves a small guard band.
    /// </summary>
    public const double LumaCutoffNyquistMargin = 0.45;

    // ── PAL frame / line geometry ────────────────────────────────────────────
    /// <summary>PAL standard: 625 lines per full interlaced frame (two fields).</summary>
    public const int PAL_LINES_PER_FRAME = 625;

    /// <summary>PAL standard: 576 active (visible) lines per frame (288 per field).</summary>
    public const int PAL_VISIBLE_LINES = 576;

    /// <summary>PAL standard frame rate: 25 frames per second.</summary>
    public const double PAL_FRAME_RATE = 25.0;

    /// <summary>Nominal duration of one PAL line: 64 µs.</summary>
    public const double PAL_LINE_DURATION = 64e-6;

    /// <summary>PAL color subcarrier frequency in Hz (4.43361875 MHz).</summary>
    public const double PAL_COLOR_CARRIER_FREQ = 4_433_618.75;

    /// <summary>PAL composite video bandwidth: 5.5 MHz.</summary>
    public const double PAL_VIDEO_BANDWIDTH = 5.5e6;

    /// <summary>
    /// Active video duration per PAL line: 64 µs total − 12 µs blanking = 52 µs.
    /// Blanking = front porch (1.5 µs) + H-sync (4.7 µs) + back porch (5.8 µs) = 12 µs.
    /// </summary>
    public const double ActiveLineSeconds = 52e-6;

    // ── FIR filter tap lengths ────────────────────────────────────────────────
    /// <summary>
    /// Tap count for the luma low-pass FIR. Keep equal to <see cref="ChromaSeparationTaps"/>
    /// so the dual-filter single-pass path (<c>ApplyTwoFiltersStaticToDest</c>) is used.
    /// Raise to 91/101 for better stop-band rejection; lower to 73 if CPU-bound.
    /// </summary>
    public const int LumaLpfTaps = 81;

    /// <summary>
    /// Tap count for the chroma band-pass FIR. Must equal <see cref="LumaLpfTaps"/> for
    /// the single-pass dual-filter optimisation to activate.
    /// </summary>
    public const int ChromaSeparationTaps = 81;

    /// <summary>Tap count for the post-demodulation chroma baseband low-pass FIR.</summary>
    public const int ChromaBasebandLpfTaps = 63;

    // ── BT.601 YUV → RGB coefficients ────────────────────────────────────────
    /// <summary>BT.601: V → R channel contribution  (R = Y + 1.402·V).</summary>
    public const double Bt601Crv = 1.402;

    /// <summary>BT.601: U → G channel contribution  (G = Y − 0.344·U − 0.714·V).</summary>
    public const double Bt601Cgu = -0.344;

    /// <summary>BT.601: V → G channel contribution.</summary>
    public const double Bt601Cgv = -0.714;

    /// <summary>BT.601: U → B channel contribution  (B = Y + 1.772·U).</summary>
    public const double Bt601Cbu = 1.772;

    // ── ColorPLL ─────────────────────────────────────────────────────────────
    /// <summary>
    /// Time (seconds) from H-sync leading edge to the start of the color burst.
    /// Equal to H-sync (4.7 µs) + Breezeway (0.6 µs) = 5.3 µs.
    /// </summary>
    public const double BurstStartSeconds = 5.3e-6;

    /// <summary>
    /// PLL proportional (P) gain. Scales the instantaneous phase error directly
    /// into a frequency correction. Higher values give faster lock but more noise sensitivity.
    /// </summary>
    public const double PllProportionalGain = 0.1;

    /// <summary>
    /// PLL integral (I) gain. Accumulates phase error over time to remove steady-state
    /// frequency offsets. Lower than proportional gain to avoid instability.
    /// </summary>
    public const double PllIntegralGain = 0.005;

    /// <summary>
    /// Oscillator renormalization interval (samples). Every 1024 steps the complex
    /// oscillator phasor is renormalized to unit magnitude to prevent floating-point drift.
    /// Must be a power of two so the modulo can be replaced by a bitwise AND.
    /// </summary>
    public const int OscillatorRenormInterval = 1024;
}
