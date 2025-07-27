using System.Numerics;
using ScottPlot;

// https://pysdr.org/content/sampling.html
public class IQSampling
{
    private readonly Plot _plot;
    public IQSampling(Plot plot)
    {
        _plot = plot;
    }

    public void CalculatePowerDpectralDensity()
    {
        double Fs = 300; // Sampling frequency in Hz
        var Ts = 1 / Fs; // sample period
        int N = 2048; // Number of points to simulate, and our FFT size

        var t = Ts * Numpy.np.arange(N);
        //var x = Numpy.np.exp(1j * 2 * Numpy.np.pi * 50 * t); // simulates sinusoid at 50 Hz
        var s = new Complex[t.len];
        for (int i = 0; i < t.len; i++)
        {
            var tcurrent = t[i].real.Get<double>();
            double angle = 2 * Numpy.np.pi * 50 * tcurrent;
            var c = Complex.Exp(new Complex(0, angle));
            s[i] = c;
        }

        //var n = (Numpy.np.random.randn(N) + 1j*Numpy.np.random.randn(N))/ Numpy.np.sqrt(2); //complex noise with unity power
        var n = new Complex[N];
        var cmpn = Numpy.np.random.randn(N);
        for (int i = 0; i < cmpn.len; i++)
        {
            var tcurrent = cmpn[i].real.Get<double>();
            n[i] = new Complex(0, 1 * tcurrent / Math.Sqrt(2));
        }
        var noise_power = 2;
        var r = (Numpy.np.random.randn(N) / Math.Sqrt(2) + Numpy.np.array(n)) * Math.Sqrt(noise_power) + Numpy.np.array(s);
        r = r * Numpy.np.hamming(r.len); // apply a Hamming window optionally
        var PSD = Numpy.np.abs(Numpy.np.fft.fft_(r)).pow(2) / (N * Fs);
        var PSD_log = 10.0 * Numpy.np.log10(PSD);
        var PSD_shifted = Numpy.np.fft.fftshift(PSD_log);

        var f = Numpy.np.arange(Fs / -2.0, Fs / 2.0, Fs / N); // start, stop, step;
        _plot.Add.SignalXY(f.GetData<double>(), PSD_shifted.GetData<double>());
        _plot.YLabel("Magnitude [dB]");
        _plot.XLabel("Frequency [Hz]");
        _plot.ShowGrid();
    }
}