using System.Drawing;
using System.Numerics;
using ScottPlot;
using Spectrogram;

// https://pysdr.org/content/frequency_domain.html
public class FrequencyDomain : Page
{
    private readonly Plot _plot;
    public FrequencyDomain(Plot plot)
    {
        _plot = plot;
    }
    public void PopulateFrequency()
    {
        double Fs = 1;// Hz
        var N = 100;// number of points to simulate, and our FFT size

        var t = Numpy.np.arange(N);
        //var t = Numpy.np.linspace(0, 100, 100);
        //var s = Numpy.np.sin(Numpy.np.array(Enumerable.Range(1,360).ToArray())*1.0* Math.PI/180 ); // generate a sine wave
        var s = Numpy.np.sin(0.15 * 2 * Numpy.np.pi * t);
        s = s * Numpy.np.hamming(N); // apply a Hamming window to the signal
        //var x = Numpy.np.linspace(-Numpy.np.pi, Numpy.np.pi, 201);
        var S = Numpy.np.fft.fftshift(Numpy.np.fft.fft_(s));
        var f = Numpy.np.arange(Fs / -2, Fs / 2, Fs / N).GetData<double>();
        var S_mag = Numpy.np.abs(S);
        var S_phase = Numpy.np.angle(S);

        _plot.Add.SignalXY(f, s.GetData<double>().ToArray());
        _plot.Add.SignalXY(f, S_mag.GetData<double>().ToArray());
        _plot.Add.SignalXY(f, S_phase.GetData<double>().ToArray());
        _plot.Axes.AutoScale();
        _plot.YLabel("Magnitude/phase/Amplitude");
        _plot.XLabel("Frequency (Hz)");
    }

    // heatmap https://scottplot.net/cookbook/4.1/category/plottable-heatmap/
    public void DrawSpectogram()
    {
        var sample_rate = 1e6;

        // Generate tone plus noise
        var t = Numpy.np.arange(1024 * 1000) / sample_rate; // time vector
        var f = 50e3;// freq of tone
        var x = (Numpy.np.sin(2 * Numpy.np.pi * f * t) + 0.2 * Numpy.np.random.randn(t.len)).GetData<double>();

        var fft_size = 1024;
        var num_rows = (int)Math.Round(x.Length * 1.0 / fft_size);// fft_size # // is an integer division which rounds down
        var spectrogram = Numpy.np.zeros((num_rows, fft_size));


        // spectrogram[i, :] = 10 * 
        // np.log10
        // (
        //  np.abs(
        //      np.fft.fftshift(
        //          np.fft.fft( 
        //              x[i * fft_size:(i + 1) * fft_size]
        //              )
        //          )
        //       ) ** 2
        // )
        for (int i = 0; i < num_rows; i++)
        {
            // 1. Get segment
            double[] segment = new double[fft_size];
            Array.Copy(x, i * fft_size, segment, 0, fft_size);

            // 2. Compute FFT
            var fft = Numpy.np.fft.fft_(segment);

            // 3. Shift zero-frequency to center
            var spectrum = Numpy.np.fft.fftshift(fft);

            // 4. Compute magnitude squared (power)
            var power = Numpy.np.multiply(Numpy.np.abs(spectrum), Numpy.np.abs(spectrum));
            //power = Numpy.np.trim_zeros(power); // remove trailing zeros
            // 5. Convert to dB
            for (int j = 0; j < fft_size; j++)
                spectrogram[i, j] = 10 * Numpy.np.log10(power[j]); // add epsilon to avoid log(0)
        }

        int numRows = spectrogram.shape[0];
        int numCols = spectrogram.shape[1];
        double[] flat = spectrogram.GetData<double>();
        double[,] spectrogramArray = new double[numRows, numCols];
        for (int i = 0; i < numRows; i++)
            for (int j = 0; j < numCols; j++)
                spectrogramArray[i, j] = flat[i * numCols + j];

        var hm = _plot.Add.Heatmap(spectrogramArray);
        _plot.YLabel("Time (s)");
        _plot.XLabel("Frequency (Hz)");
        hm.Axes.XAxis.Min = sample_rate / -2 / 1e6;
        hm.Axes.XAxis.Max = sample_rate / 2 / 1e6;
        hm.Axes.YAxis.Min = 0;
        hm.Axes.YAxis.Max = x.Length / sample_rate;
        _plot.Axes.AutoScale();
    }

    public void FFTSimulation()
    {
        // Simulate a tone + noise
        var sample_rate = 1e6;
        var f_offset = 0.2e6; // 200 kHz offset from carrier
        var N = 1024;
        var t = Numpy.np.arange(N) / sample_rate;
        var s = new Complex[t.len];
        //s = Numpy.np.exp(new Complex(1, 2) * Numpy.np.pi * f_offset * t);

        for (int i = 0; i < t.len; i++)
        {
            var tcurrent = t[i].real.Get<double>();
            double angle = 2 * Math.PI * f_offset * tcurrent;
            var c = Complex.Exp(new Complex(0, angle));
            s[i] = c;
        }

        // Generate real and imaginary noise parts
        var n_real = Numpy.np.random.randn(N).GetData<double>();
        var n_imag = Numpy.np.random.randn(N).GetData<double>();
        var n = new Complex[N];
        for (int i = 0; i < N; i++)
            n[i] = new Complex(n_real[i], n_imag[i]) / Math.Sqrt(2);

        //var n = (Numpy.np.random.randn(N) + new Complex(0, 1) * Numpy.np.random.randn(N)) / Math.Sqrt(2); // unity complex noise
        var r = Numpy.np.array(s) + Numpy.np.array(n); // 0 dB SNR

        // Perform fft, fftshift, convert to dB

        // var X = Numpy.np.fft.fft_(r);
        var X = fft(r);
        var X_shifted = Numpy.np.fft.fftshift(X);//# equivalent to np.fft.fftshift
        var X_mag = 10 * Numpy.np.log10(Numpy.np.abs(X_shifted) * Numpy.np.abs(X_shifted));

        //Plot results
        var f = Numpy.np.linspace(sample_rate / -2, sample_rate / 2, N) / 1e6;// plt in MHz
        var hm = _plot.Add.SignalXY(f.GetData<double>(), X_mag.GetData<double>());
        _plot.Grid.IsVisible = true;
        //_plot.Axes.AutoScale();
        _plot.XLabel("Frequency [MHz]");
        _plot.YLabel("Magnitude [dB]");

    }

private Complex[] ToComplex(Numpy.NDarray x)
    {
        var data_r = x.real.GetData<double>();
        var data_im = x.imag.GetData<double>();
        var complexArray = new Complex[data_r.Length];
        for (int i = 0; i < data_r.Length; i++)
            complexArray[i] = new Complex(data_r[i], data_im[i]);
        return complexArray;
    }
    private Numpy.NDarray fft(Numpy.NDarray x)
    {
        var N = x.len;
        if (N == 1)
            return x;
        //var twiddle_factors = Numpy.np.exp(-new Complex(0, 2) * Numpy.np.pi * Numpy.np.arange(N / 2) / N);
        int halfN = N / 2;
        Complex[] twiddle_factors = new Complex[halfN];
        for (int k = 0; k < halfN; k++)
        {
            double angle = -2 * Math.PI * k / N;
            twiddle_factors[k] = Complex.Exp(new Complex(0, angle));
        }

        var x_even = fft(x["::2"]); // yay recursion!

        var x_odd = fft(x["1::2"]);
        return Numpy.np.concatenate([x_even + Numpy.np.array(twiddle_factors) * x_odd,
                           x_even - Numpy.np.array(twiddle_factors) * x_odd]);
    }
}