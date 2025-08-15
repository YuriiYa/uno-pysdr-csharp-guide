using ScottPlot;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.InteropServices;

// Numpy.NET documentation: https://github.com/SciSharp/Numpy.NET
public class VideoFromFile
{
    private readonly Plot _plot;
    private readonly Plot _plotModulation;
    public VideoFromFile(Plot plot, Plot plotModulation)
    {
        _plot = plot;
        _plotModulation = plotModulation;
    }

    public static Complex[] ConvertToIQ(Int16[] data)
    {
        var halfLength = data.Length / 2;
        Complex[] ret = new Complex[halfLength];
        for (int i = 0; i < ret.Length; i++)
        {
            ret[i] = new Complex(data[2 * i], data[2 * i + 1]) / Int16.MaxValue; // even are real and odd are img numbers and adjuct it to -1 to +1;
        }
        return ret;
    }

    public void DoShowVideoFile()
    {
        // loaded from here https://www.sigidwiki.com/wiki/PAL_Broadcast 
        var filename = @"C:\projects\sdr\RFData\SDRSharp_20170122_171736Z_179100000Hz_IQ.wav";
        var samples = Numpy.np.fromfile(filename, Numpy.np.int16);

        var int16Header = samples[":22"].real.GetData<Int16>();
        WavHeader header = WavHelper.ConvertByteArraytoType<WavHeader>(int16Header);

        samples = samples["22:"];

        for (int i = 0; i < 10; i++)
        {
            Console.WriteLine(samples[i]);
        }

        var iq_samples = Numpy.np.array(ConvertToIQ(samples.real.GetData<Int16>()));

        var sample_rate = header.SampleRate;
        var center_freq = 179100000;
        double CenterFreqinMHz = center_freq ;

        // Take one batch of samples equivalent to 8 MHz bandwidth
        // For proper frequency resolution, use a power of 2 that covers the bandwidth well
        var batch_size = 8192; // Good FFT size for 8 MHz bandwidth
        var batch_samples = iq_samples[$":{batch_size}"];
        //  Apply window function to reduce spectral leakage
        var windowed_samples = batch_samples * Numpy.np.hanning(batch_samples.len);
        //  Apply FFT to convert to frequency domain
        var fft_result = Numpy.np.fft.fft_(windowed_samples);
        var fft_shifted = Numpy.np.fft.fftshift(fft_result); // Center the spectrum
        var freqs = Numpy.np.fft.fftshift(Numpy.np.fft.fftfreq(batch_size, (float)1 / sample_rate));
        var actual_freqs = center_freq + freqs; // Actual RF frequencies

        // Convert to power spectral density (magnitude squared, in dB)
        var psd = 20 * Numpy.np.log10(Numpy.np.abs(fft_shifted) + 1e-12);  // Add small value to avoid log(0)
                                                                           // Define the heatmap's boundaries using its Extent
        var left = CenterFreqinMHz + sample_rate / (-2.0) ;
        var right = CenterFreqinMHz + sample_rate / 2.0 ;

        int num_rows = samples.len / batch_size;
        var spectrogram = Numpy.np.zeros((num_rows, batch_size));

        for (int i = 0; i < num_rows; i++)
            spectrogram[$"{i}, :"] = 10 * Numpy.np.log10(
                (
                    Numpy.np.abs(
                        Numpy.np.fft.fftshift(
                            Numpy.np.fft.fft_(
                                //  Apply window function to reduce spectral leakage
                                samples[$"{i * batch_size}:{(i + 1) * batch_size}"] * Numpy.np.hanning(batch_size)
                        )
                    ) + 1e-12 // Add small value to avoid log(0)
                )
                ).pow(2));

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

        // Define the heatmap's boundaries using its Extent
        var bottom = 0;
        var top = samples.len / ((double)sample_rate / 1e6);

        hm.Extent = new ScottPlot.CoordinateRect(
            left: left,
            right: right,
            bottom: bottom,
            top: top
         );

        _plotModulation.Add.SignalXY(actual_freqs.GetData<double>(), psd.real.GetData<double>());
        _plotModulation.YLabel("Power");
        _plotModulation.XLabel("Frequency (Hz)");
        _plotModulation.Axes.SetLimitsX(left, right);

        _plot.Axes.InvertY();
        _plot.Axes.SetLimitsX(left, right);
        //_plot.Axes.SetLimitsY(bottom, top);
        _plot.Axes.MarginsX(0.1);
        _plotModulation.Axes.MarginsX(0.1);

        PixelPadding padding = new(50, 20, 30, 5);
        _plot.Layout.Fixed(padding);
        _plotModulation.Layout.Fixed(padding);


        _plot.Axes.Link(_plotModulation, x: true, y: false);
        _plotModulation.Axes.Link(_plot, x: true, y: false);
        //_plotModulation.Axes.AutoScale();
        _plotModulation.PlotControl?.Refresh();
        _plot.Axes.AutoScale();
        _plot.PlotControl?.Refresh();


    }
}