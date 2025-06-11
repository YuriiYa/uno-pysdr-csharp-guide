using System.Drawing;
using Spectrogram;

// https://pysdr.org/content/frequency_domain.html
// TODO: DELETE
public class Spectograph
{
    public Spectograph()
    {
    }


    public string DrawSpectogram()
    {
        var sample_rate = 1e6;

        // Generate tone plus noise
        var t = Numpy.np.arange(1024 * 1024) / sample_rate; // time vector
        var f = 50e3;// freq of tone
        var x = Numpy.np.sin(2 * Numpy.np.pi * f * t) + 0.2 * Numpy.np.random.randn(t.len);

        var fft_size = 1024;
        var num_rows = x.len;// fft_size # // is an integer division which rounds down
        var spectrogram = Numpy.np.zeros((num_rows, fft_size));
        int stepSize = x.len / fft_size;

        double[] signal = FftSharp.SampleData.SampleAudio1();
        var window = new FftSharp.Windows.Hanning();
        window.ApplyInPlace(signal);
        System.Numerics.Complex[] spectrum = FftSharp.FFT.Forward(signal);
        double[] magnitude = FftSharp.FFT.Magnitude(spectrum);

        var sg = new SpectrogramGenerator(48000, fftSize: 1024, stepSize: 20, maxFreq: 1000);

        sg.Add(magnitude);
        var path = @"C:\projects\sdr\HackRFTutorial\uno_palyground\uno_palyground\bin\Debug\net9.0-desktop\tmp.png";
        sg.SaveImage(path, intensity: 5, dB: true);
        return path;
    }
}