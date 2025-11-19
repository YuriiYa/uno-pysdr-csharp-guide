using System.Numerics;
using Microsoft.UI.Dispatching;
using nethackrf;
using ScottPlot;

namespace HackRF.Namespace;

// https://pysdr.org/content/hackrf.html
// to find infromation about connected HackRF devices use command line tool
// 'hackrf_info'
public class HackRF
{
    private readonly Plot _plot;
    private readonly Plot _plotModulation;
    private readonly HackRFInteraction _hackrfInteraction;
    private HackRFConfiguration _configuration;
    private readonly DispatcherQueue _dispatcherQueue;

    private bool _stopRequested = false;
    private CancellationTokenSource? _cts;
    private Task? _readerTask;

    public HackRF(Plot plotSpectrogram, DispatcherQueue dispatcherQueue, Plot plotModulation, HackRFConfiguration configuration, HackRFInteraction hackRFInteraction)
    {
        _dispatcherQueue = dispatcherQueue;
        _plot = plotSpectrogram;
        _plotModulation = plotModulation;
        _configuration = configuration;
        _hackrfInteraction = hackRFInteraction;
    }

    public HackRFConfiguration GetDfaultConfig(HackRFInteraction hackRFInteraction)
    {
        //These settings should match the hackrf_transfer example used in the textbook, and the resulting waterfall should look about the same
        var recording_time = 1; // seconds
        var center_freq = 1000 * 1000 * 100;  // 100 MHz
        var sample_rate = 1000 * 1000 * 10; //10 MHz
        var numSamples = (int)recording_time * sample_rate;//(int)recording_time * sample_rate; // 10 million samples for 1 second at 10 MHz
        var baseband_filter = 7.5e6;
        var lna_gain = 30; // 0 to 40 dB in 8 dB steps
        var vga_gain = 50; // 0 to 62 dB in 2 dB steps

        _configuration = new HackRFConfiguration(
                   recordingTimeInSec: recording_time,
                   deviceSerialNumber: "0000000000000000930464dc242ea317",
                   centerFrequency: center_freq,
                   sampleFrequency: sample_rate,
                   numSamples: numSamples,
                   lnaGain: lna_gain, // dB
                   vgaGain: vga_gain, // dB
                   filterBandwidth: baseband_filter, // MHz
                   fftSize: 2048 * 16
               );
        return _configuration;
    }

    public void Start()
    {
        if (_readerTask != null) return; // already running
        _cts = new CancellationTokenSource();
        _stopRequested = false;
        _readerTask = Task.Run(() => RunCaptureLoop(_cts.Token), _cts.Token);
    }

    public async Task StopAsync()
    {
        _stopRequested = true;
        if (_cts != null)
        {
            _cts.Cancel();
        }
        if (_readerTask != null)
        {
            try { await _readerTask.ConfigureAwait(false); } catch (OperationCanceledException) { }
        }
        _readerTask = null;
        _cts?.Dispose();
        _cts = null;
        _hackrfInteraction.Stop();
    }

    private void RunCaptureLoop(CancellationToken ct)
    {
        var isStarted = _hackrfInteraction.Setup(
            _configuration
        );
        if (!isStarted || !_hackrfInteraction.Start())
        {
            var error ="Failed to start HackRF interaction.";
            System.Console.WriteLine(error);
            throw new EntryPointNotFoundException(error);
        }
        try
        {
            while (!ct.IsCancellationRequested && !_stopRequested)
            {
                _hackrfInteraction.ReadFromHackRF(OnDataReceived);
            }
        }
        catch (Exception ex)
        {
            System.Console.WriteLine($"HackRF capture loop error: {ex.Message}");
        }
    }


    private void OnDataReceived(Complex[] values)
    {
        //var samples = Numpy.np.array(values)["100000:"]; // get rid of the first 100k samples just to be safe, due to transients
        var samples = Numpy.np.array(values);

        var fft_size = _configuration.FFTSize;
        int num_rows = samples.len / fft_size;
        var spectrogram = Numpy.np.zeros((num_rows, fft_size));

        for (int i = 0; i < num_rows; i++)
            spectrogram[$"{i}, :"] = 10 * Numpy.np.log10(Numpy.np.abs(Numpy.np.fft.fftshift(Numpy.np.fft.fft_(samples[$"{i * fft_size}:{(i + 1) * fft_size}"]))).pow(2));

        double CenterFreqinMHz = _configuration.CenterFrequency / 1e6;

        var lengthValues = values.Length;
        var lengthAbsFreqs = samples.len;

        int numRows = spectrogram.shape[0];
        int numCols = spectrogram.shape[1];
        double[] flat = spectrogram.GetData<double>();
        double[,] spectrogramArray = new double[numRows, numCols];
        for (int i = 0; i < numRows; i++)
            for (int j = 0; j < numCols; j++)
                spectrogramArray[i, j] = flat[i * numCols + j];
        var power = spectrogram["-1"].real.GetData<double>();
        var left = CenterFreqinMHz + _configuration.SampleFrequency / (-2.0) / 1e6;
        var right = CenterFreqinMHz + _configuration.SampleFrequency / 2.0 / 1e6;
        var bottom = 0;
        var top = samples.len / ((double)_configuration.SampleFrequency / 1e6);
        double[] freq = Numpy.np.arange(left, right, (right - left) / spectrogram["-1"].len).GetData<double>();


        _ = _dispatcherQueue.TryEnqueue(() =>
        {
            _plot.Clear();
            _plotModulation.Clear();
            var hm = _plot.Add.Heatmap(spectrogramArray);
            _plot.YLabel("Time (s)");
            _plot.XLabel("Frequency (Hz)");

            // Define the heatmap's boundaries using its Extent
            hm.Extent = new ScottPlot.CoordinateRect(
                left: left,
                right: right,
                bottom: bottom,
                top: top
             );

            _plot.Axes.InvertY();
            _plot.Axes.SetLimitsX(left, right);
            _plot.Axes.SetLimitsY(bottom, top);

            _plot.Axes.MarginsX(0.1);

            _plotModulation.Add.SignalXY(freq, power);
            _plotModulation.YLabel("Power");
            _plotModulation.XLabel("Frequency (Hz)");
            _plotModulation.Axes.SetLimitsX(left, right);

            PixelPadding padding = new(50, 20, 30, 5);
            _plot.Layout.Fixed(padding);
            _plotModulation.Layout.Fixed(padding);

            _plot.Axes.Link(_plotModulation, x: true, y: false);
            _plotModulation.Axes.Link(_plot, x: true, y: false);

            _plotModulation.Axes.AutoScale();
            _plot.Axes.AutoScale();
            _plot.PlotControl?.Refresh();
            _plotModulation.PlotControl?.Refresh();
        });
    }
}

public class HackRFConfiguration
{
    public int RecordingTimeInSec { get; set; }
    public string DeviceSerialNumber { get; set; }
    public int CenterFrequency { get; set; }
    public int SampleFrequency { get; set; }
    public int NumSamples { get; set; }
    public int LnaGain { get; set; }
    public int VgaGain { get; set; }
    public double FilterBandwidth { get; set; }

    public int FFTSize { get; set; }

    public HackRFConfiguration(int recordingTimeInSec, string deviceSerialNumber, int centerFrequency, int sampleFrequency, int numSamples, int lnaGain, int vgaGain, double filterBandwidth, int fftSize)
    {
        RecordingTimeInSec = recordingTimeInSec;
        DeviceSerialNumber = deviceSerialNumber;
        CenterFrequency = centerFrequency;
        SampleFrequency = sampleFrequency;
        NumSamples = numSamples;
        LnaGain = lnaGain;
        VgaGain = vgaGain;
        FilterBandwidth = filterBandwidth;
        FFTSize = fftSize;
    }
}

public class HackRFInteraction
{
    private NetHackrf? _device;
    private HackRFConfiguration? _hackRFConfiguration;
    private Stream? _dataStream;
    public int NumSamples { get; private set; }
    //public int SampleFrequency { get; private set; }

    static int[] max2837_ft = [
            1750000,
            2500000,
            3500000,
            5000000,
            5500000,
            6000000,
            7000000,
            8000000,
            9000000,
            10000000,
            12000000,
            14000000,
            15000000,
            20000000,
            24000000,
            28000000];
    public bool Setup(HackRFConfiguration hackRFConfiguration)
    {
        _hackRFConfiguration = hackRFConfiguration;
        var devices = NetHackrf.HackrfDeviceList(); // get list of all connected hackrf transceivers
        if (devices.Length == 0) // if no hackrfs discovered
        {
            System.Console.WriteLine("No hackrf devices were found");
            return false;
        }

        var foundDeviceByName = devices.FirstOrDefault(device => device.serial_number == _hackRFConfiguration.DeviceSerialNumber);
        if (foundDeviceByName == null)
        {
            var notfoundError = $"No hackrf devices were found with the specified serial number {_hackRFConfiguration.DeviceSerialNumber} will be using the first one in the list";
            System.Console.WriteLine(notfoundError);
            return false;
        }

        _device = (foundDeviceByName ?? devices[0]).OpenDevice(); // connecting to the first transceiver in the list

        _device.CarrierFrequencyMHz = _hackRFConfiguration.CenterFrequency / 1e6;//  in  MHz
        NumSamples = _hackRFConfiguration.NumSamples;
        //SampleFrequency = (int)_hackRFConfiguration.SampleFrequency;
        _device.SampleFrequencyMHz = _hackRFConfiguration.SampleFrequency / 1e6;//  in  MHz
        _device.LNAGainDb = _hackRFConfiguration.LnaGain;
        _device.VGAGainDb = _hackRFConfiguration.VgaGain;
        _device.FilterBandwidthMHz = ComputeBasebandFilterBwRoundDownLt((int)(_hackRFConfiguration.FilterBandwidth * 1e6)) / 1e6;
        _device.AMPEnable = false;
        //_device.Reset(); // reset hackRF to apply all settings
        return true;
    }

    public bool Start(
    )
    {
        if (_device == null)
            return false;

        _dataStream = _device?.StartRX();

        return true;
    }

    public void Stop()
    {
        try
        {
            _dataStream?.Dispose();
            _dataStream = null;
            _device?.Reset();
        }
        catch (Exception ex)
        {
            System.Console.WriteLine($"Error stopping HackRF: {ex.Message}");
        }
    }

    public void ReadFromHackRF(Action<Complex[]> onDataRecieved)
    {

        byte[] buffer;
        buffer = new byte[NumSamples];

        try
        {
            if (_dataStream == null || !_dataStream.CanRead) return;
            var valid_length = _dataStream.Read(buffer, 0, buffer.Length); // reading interpolated IQ data from stream
            System.Console.WriteLine("demodulating...");
            var IQ_data = ConvertToIQ(buffer, valid_length); // converting interpolated IQ data to complex array
            onDataRecieved(IQ_data);

        }
        catch (Exception ex)
        {
            System.Console.WriteLine($"Error while reading from HackRF: {ex.Message}");
            _device?.Reset(); // reset hackRF if something goes wrong
        }
    }

    static Complex[] ConvertToIQ(byte[] buffer, int valid_length)
    {
        var halfLength = valid_length / 2;
        Complex[] ret = new Complex[halfLength];
        for (int i = 0; i < ret.Length; i++)
        {
            ret[i] = new Complex((sbyte)buffer[2 * i], (sbyte)buffer[2 * i + 1]) / 128; // even are real and odd are img numbers and adjuct it to -1 to +1;
        }
        return ret;
    }


    static double[] coefs = { -0.006052, -0.005539, -0.007277, -0.008615, -0.009164, -0.008511, -0.006271, -0.002146, 0.004031, 0.012244, 0.022290, 0.033752, 0.046049, 0.058433, 0.070094, 0.080216, 0.088051, 0.093009, 0.094705, 0.093009, 0.088051, 0.080216, 0.070094, 0.058433, 0.046049, 0.033752, 0.022290, 0.012244, 0.004031, -0.002146, -0.006271, -0.008511, -0.009164, -0.008615, -0.007277, -0.005539, -0.006052 };
    static Complex[] prev_val = new Complex[coefs.Length];
    static Complex[] LPF1(Complex[] data) // FIR filter (36 order, 40kHz passband edge, 150kHz stopband edge, 40dB stopband att)
    {
        for (int i = 0; i < data.Length; i++)
        {
            Array.Copy(prev_val, 1, prev_val, 0, prev_val.Length - 1);
            prev_val[prev_val.Length - 1] = data[i];
            Complex output = Complex.Zero;
            for (int j = 0; j < prev_val.Length; j++)
            {
                output += prev_val[j] * coefs[j];
            }
            data[i] = output;
        }
        return data;
    }


    static double[] coefs2 = { 0.0, 0.001741, -0.015102, 0.007302, 0.041937, -0.057247, -0.070744, 0.298657, 0.583333, 0.298657, -0.070744, -0.057247, 0.041937, 0.007302, -0.015102, 0.001741 };
    static double[] prev_val2 = new double[16];
    static double[] LPF2(double[] data) // FIR filter (16 order, 10kHz passband edge, 17kHz stopband edge, 50dB stopband att)
    {
        for (int i = 0; i < data.Length; i++)
        {
            Array.Copy(prev_val2, 1, prev_val2, 0, 15);
            prev_val2[15] = data[i];
            double output = 0;
            for (int j = 0; j < 16; j++)
            {
                output += prev_val2[j] * coefs2[j];
            }
            data[i] = output;
        }
        return data;
    }

    static void AddWaveHeader(Stream s, Int32 size) // wave header writer. Described at http://soundfile.sapp.org/doc/WaveFormat/
    {
        byte[] header = {   0x52, 0x49, 0x46, 0x46, 0x24, 0x00, 0x00, 0x80, 0x57, 0x41, 0x56, 0x45, // RIFF chunk descriptor
                                0x66, 0x6d, 0x74, 0x20, 0x10, 0x00, 0x00, 0x00, 0x01, 0x00, 0x01, 0x00, // fmt sub-chunk PCM mono 
                                0x80, 0xBB, 0x00, 0x00, 0x00, 0x77, 0x01, 0x00, 0x02, 0x00, 0x10, 0x00, // 48000 sample rate, 96000 byte rate, 16 bits per sample
                                0x64, 0x61, 0x74, 0x61, 0x00, 0x00, 0x00, 0x80}; // DATA
        Array.Copy(BitConverter.GetBytes(size + 36), 0, header, 4, 4);
        Array.Copy(BitConverter.GetBytes(size), 0, header, 40, 4);
        s.Write(header, 0, 44);
    }

    static double prev_I = 0;
    static double prev_Q = 0;
    static double[] DemodFMsamples(Complex[] IQ, double Fs)
    {
        double[] ret = new double[IQ.Length];
        double I;
        double Q;
        for (int i = 0; i < ret.Length; i++)
        {
            I = IQ[i].Real;
            Q = IQ[i].Imaginary;
            double dI = (I - prev_I) * Fs;
            double dQ = (Q - prev_Q) * Fs;

            ret[i] = (dI * Q - dQ * I) / (I * I + Q * Q); // https://ru.dsplib.org/content/signal_fm_demod/img_html/fmdemod_html_46eb685f.gif

            prev_I = I;
            prev_Q = Q;
        }
        return ret;
    }

    static double[] Decimate(double[] data, int K)
    {
        double[] ret = new double[data.Length / K];
        for (int i = 0; i < ret.Length; i++)
        {
            ret[i] = data[K * i];
        }
        return ret;
    }

    static byte[] ConvertToBuffer(double[] data)
    {
        byte[] buffer = new byte[data.Length * 2];
        for (int i = 0; i < data.Length; i++)
        {
            Int16 sample = (short)(data[i]);
            var bytes = BitConverter.GetBytes(sample);
            buffer[i * 2] = bytes[0];
            buffer[i * 2 + 1] = bytes[1];
        }
        return buffer;
    }

    // https://github.com/greatscottgadgets/hackrf/blob/4b8dbfc308bc2b7abfec614b8bd988ad33136110/host/libhackrf/src/hackrf.c#L1014
    private int ComputeBasebandFilterBwRoundDownLt(int bandwidthHz)
    {
        int i = 0;
        // Find the first entry >= requested bandwidth
        while (i < max2837_ft.Length && max2837_ft[i] != 0)
        {
            if (max2837_ft[i] >= bandwidthHz)
                break;
            i++;
        }

        // Round down (if not the first entry)
        if (i != 0)
            i--;
        return max2837_ft[i];
    }

}