using System.Numerics;
using ScottPlot;

public class Noise
{
    private readonly Plot _plot;
    private readonly Plot _plotModulation;
    public Noise(Plot plot, Plot plotModulation)
    {
        _plot = plot;
        _plotModulation = plotModulation;
    }

    public void GaussianNoise()
    {
        int N = 1000; // number of points to simulate
        var x = Numpy.np.random.randn(N);
        // Create time vector
        var t = Numpy.np.arange(0, N, 1).real;


        float mean = 0.0f;
        float stddev = 1.0f;

        // Generate Gaussian noise
        var noise = Numpy.np.random.normal(mean, stddev, N).real.GetData<double>();
        // only look at positive frequencies.  remember // is just an integer divide
        var X = Numpy.np.fft.fftshift(Numpy.np.fft.fft_(x));
        X = X.real / Numpy.np.max(X);

        _plot.Add.Signal(x.real.GetData<double>(), 1, Color.FromColor(System.Drawing.Color.Red));
        _plot.Add.Signal(X.GetData<double>(), 1, Color.FromColor(System.Drawing.Color.Blue));
        _plot.Add.Signal(noise, 1, Color.FromColor(System.Drawing.Color.Green));
        _plot.YLabel("Amplitude");
        _plot.XLabel("Sample Index");
        _plot.Title("Gaussian Noise with FFT");
        _plot.Axes.AutoScale();

        _plot.PlotControl?.Refresh();

    }

    public void ComplexNoise()
    {
        int N = 1000; // number of points to simulate
        // Create time vector
        var guasian_noise = new Complex[N];
        var nsigq = Numpy.np.random.randn(N);
        var nsigi = Numpy.np.random.randn(N);
        // AWGN with unity power
        var noise_power = 0.1;

        for (var i = 0; i < N; i++)
        {
            guasian_noise[i] = new Complex(nsigq.GetData<double>()[i], nsigi.GetData<double>()[i]) * Math.Sqrt(noise_power) / Math.Sqrt(2);
        }
        
        var complexNoise = Numpy.np.array(guasian_noise);

        _plot.Add.Signal(complexNoise.real.GetData<double>(), N, Color.FromColor(System.Drawing.Color.Blue));
        _plot.Add.Signal(complexNoise.imag.GetData<double>(), N, Color.FromColor(System.Drawing.Color.Red));
        _plot.YLabel("Amplitude");
        _plot.XLabel("Sample Index");
        _plot.Title("Complex Noise");
        var realcolor = Color.FromColor(System.Drawing.Color.Blue);
        var imagcolor = Color.FromColor(System.Drawing.Color.Red);
        _plot.Legend = new Legend(_plot)
        {
            Alignment = Alignment.UpperRight,
            FontSize = 10,
            ManualItems = new  List<LegendItem>(new[]
            {
                new LegendItem(){ LabelText = "Real Part", MarkerColor = realcolor, LineColor  = realcolor, LabelBorderColor = realcolor, FillColor = realcolor, LabelFontColor = realcolor},
                new LegendItem(){ LabelText = "Imaginary Part", MarkerColor = imagcolor, LineColor  = imagcolor, LabelBorderColor = imagcolor, FillColor = imagcolor, LabelFontColor = imagcolor},
            })
        };

        _plot.Axes.AutoScale();

        _plot.PlotControl?.Refresh();

    }
}