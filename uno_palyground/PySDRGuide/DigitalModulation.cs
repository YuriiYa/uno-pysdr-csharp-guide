using System.Numerics;
using ScottPlot;
using ScottPlot.DataSources;

public class DigitalModulation
{
    private readonly Plot _plot;

    public DigitalModulation(Plot plot)
    {
        _plot = plot;
    }

    public void QPSK()
    {
        var num_symbols = 1000;
        var x_symbols = Numpy.np.random.randint(0, 4, new int[] { num_symbols }); // 0 to 3
        //var x_int = Numpy.np.array<double>(new double[] { 0, 1, 2, 3 });
        var x_degrees = x_symbols * 360 / 4.0 + 45; // 45, 135, 225, 315 degrees
        var x_radians = x_degrees * Numpy.np.pi / 180.0; // sin() and cos() takes in radians
        var img = Numpy.np.sin(x_radians);
        var real = Numpy.np.cos(x_radians);

        var signal = new Complex[x_symbols.len];
        // AWGN with unity power
        var noise_power = 0.1;
        var phase_noise = Numpy.np.random.randn(x_symbols.len) * noise_power; //adjust multiplier for "strength" of phase noise
                                                                              //var dd = x_symbols * Numpy.np.exp(1j * phase_noise);
        var iqphase_noise = new Complex[x_symbols.len];
        for (var i = 0; i < x_symbols.len; i++)
        {
            iqphase_noise[i] = new Complex(0, phase_noise.GetData<double>()[i]);
            signal[i] = new Complex(real.GetData<double>()[i], img.GetData<double>()[i]);
        }
        var npiq_phase_noise = Numpy.np.exp(Numpy.np.array(iqphase_noise));
        //var r = Numpy.np.array(sig);
        //var x_symbols = Numpy.np.cos(x_radians) + new Complex(0, 1) * Numpy.np.sin(x_radians); // this produces our QPSK complex symbols
        var guasian_noise = new Complex[num_symbols];
        var nsigq = Numpy.np.random.randn(num_symbols);
        var nsigi = Numpy.np.random.randn(num_symbols);

        for (var i = 0; i < num_symbols; i++)
        {
            guasian_noise[i] = new Complex(nsigq.GetData<double>()[i], nsigi.GetData<double>()[i]) * Math.Sqrt(noise_power) / Math.Sqrt(2);
        }

        var sig = Numpy.np.array(signal) * npiq_phase_noise + Numpy.np.array(guasian_noise);
        var x = sig.real.GetData<double>();
        var y = sig.imag.GetData<double>();
        var hm = _plot.Add.Markers(x, y);
        // _plot.Add.Signal(y);
        IColormap cmap = new ScottPlot.Colormaps.MellowRainbow();
        hm.Color = new Color(0, 255, 255, 90);

        //_plot.Grid.IsVisible = true;
        _plot.Axes.AutoScale();
        _plot.PlotControl.Refresh();
        _plot.XLabel("I");
        _plot.YLabel("Q");

    }
}