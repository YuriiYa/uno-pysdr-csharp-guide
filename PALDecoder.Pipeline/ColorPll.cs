using System.Numerics;
using System.Runtime.CompilerServices;

namespace uno_palyground.PySDRGuide.Pipeline;

public class ColorPll
{
    private readonly double _sampleRate;
    private readonly double _burstStartTime;
    private readonly double _burstDuration;

    private double _frequency = PalConstants.PAL_COLOR_CARRIER_FREQ;
    private double _phaseErrorIntegrator = 0.0;
    private double _phaseIncrement;
    private Complex _osc = Complex.One;
    private Complex _rot;
    private int _renormCounter;

    private static readonly Complex PhasePos = Complex.FromPolarCoordinates(1, -Math.PI / 4);
    private static readonly Complex PhaseNeg = Complex.FromPolarCoordinates(1, -3 * Math.PI / 4);

    public ColorPll(int sampleRate)
    {
        _sampleRate = sampleRate;
        _burstStartTime = PalConstants.BurstStartSeconds;
        _burstDuration = PalConstants.BurstDurationSeconds;
        UpdatePhaseIncrement();
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void UpdatePhaseIncrement()
    {
        _phaseIncrement = 2 * Math.PI * _frequency / _sampleRate;
        _rot = Complex.FromPolarCoordinates(1, _phaseIncrement);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void StepOscillator()
    {
        _osc *= _rot;
        if ((++_renormCounter & (PalConstants.OscillatorRenormInterval - 1)) == 0)
        {
            double mag = _osc.Magnitude;
            _osc /= mag;
        }
    }

    public Complex GetReference(Complex chromaSample, int sampleIndexInLine, bool isVInverted)
    {
        double timeInLine = sampleIndexInLine / _sampleRate;

        if (timeInLine >= _burstStartTime && timeInLine < _burstStartTime + _burstDuration)
        {
            var burstReference = _osc * (isVInverted ? PhaseNeg : PhasePos);
            double phaseError = (chromaSample * Complex.Conjugate(burstReference)).Phase;

            _phaseErrorIntegrator += phaseError * PalConstants.PllIntegralGain;
            _frequency = PalConstants.PAL_COLOR_CARRIER_FREQ + phaseError * PalConstants.PllProportionalGain + _phaseErrorIntegrator;
            UpdatePhaseIncrement();
        }
        StepOscillator();
        return _osc;
    }
}
