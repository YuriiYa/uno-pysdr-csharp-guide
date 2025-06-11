using System.Numerics;
using System.Runtime.Serialization;
using nethackrf;
using ScottPlot;

namespace uno_palyground.Presentation;

public sealed partial class MainPage : Page
{

    private double SampleFrequency = 1024 * 1024 * 10;
    private int FFTSize = 1024 * 1024 * 8; // FFT size (number of complex samples)
    private int CenterFrequency = 99000000;    // 105.5 MHz
    HackRFInteraction _hackrfInteraction;
    public MainPage()
    {
        this.InitializeComponent();

        var freqDomainLab = new FrequencyDomain(WinUIPlot1.Plot);
        //freqDomainLab.PopulateFrequency();
        // freqDomainLab.DrawSpectogram();
        freqDomainLab.FFTSimulation();
        WinUIPlot1.Refresh();

        //Image1.Source = path;

        /* WinUIPlot1.Plot.YLabel("Power");
                    WinUIPlot1.Plot.XLabel("Frequency (Hz)");
                    //WinUIPlot1.Plot.Axes.GetXAxes().First().Range.Set(95e6, 200);
                    Task.Factory.StartNew(() =>
                     {
                         _hackrfInteraction = new HackRFInteraction();
                         var isStarted = _hackrfInteraction.Setup(
                              deviceSerialNumber: "0000000000000000930464dc242ea317",
                              centerFrequency: CenterFrequency,
                              //sampleFrequency: 48000 * 42, // 48kHz * 42
                              sampleFrequency: SampleFrequency,
                              fftSize: FFTSize,
                              lnaGain: 20, // dB
                              vgaGain: 40, // dB
                              filterBandwidth: 0.2 // MHz
                          );

                         if (!isStarted || !_hackrfInteraction.Start())
                         {
                             System.Console.WriteLine("Failed to start HackRF interaction.");
                             return;
                         }

                         _hackrfInteraction.ReadFromHackRF(OnDataReceived);

                     }, TaskCreationOptions.LongRunning);
            */

    }

    private void OnDataReceived(double[] values)
    {
        DispatcherQueue.TryEnqueue(() =>
        {
            double[] freq = FftSharp.FFT.FrequencyScale(values.Length, SampleFrequency, positiveOnly: true);
            WinUIPlot1.Plot.Clear();
            double CenterFreqinMHz = CenterFrequency;
            double[] absFreqs = freq.Select(f => f + CenterFreqinMHz).ToArray();

            var lengthValues = values.Length;
            var lengthAbsFreqs = absFreqs.Length;
            WinUIPlot1.Plot.Add.SignalXY(absFreqs, values);
            WinUIPlot1.Plot.Axes.AutoScale();
            WinUIPlot1.Plot.YLabel("Power");
            WinUIPlot1.Plot.XLabel("Frequency (Hz)");
            WinUIPlot1.Refresh();
        });
    }


}

public class HackRFInteraction
{
    private NetHackrf? _device;
    private Stream? _dataStream;
    public int FftSize { get; private set; }
    public int SampleFrequency { get; private set; }
    public bool Setup(
        string deviceSerialNumber,
        int centerFrequency,
        double sampleFrequency,
        int fftSize,
        int lnaGain,
        int vgaGain,
        double filterBandwidth
    )
    {
        var devices = NetHackrf.HackrfDeviceList(); // get list of all connected hackrf transceivers
        if (devices.Length == 0) // if no hackrfs discovered
        {
            System.Console.WriteLine("No hackrf devices were found");
            return false;
        }

        var foundDeviceByName = devices.FirstOrDefault(device => device.serial_number == deviceSerialNumber);
        if (foundDeviceByName == null)
        {
            System.Console.WriteLine($"No hackrf devices were found with the specified serial number {deviceSerialNumber} will be using the first one in the list");
        }

        _device = (foundDeviceByName ?? devices[0]).OpenDevice(); // connecting to the first transceiver in the list

        _device.CarrierFrequencyMHz = centerFrequency / 1000000;//  in  MHz
        FftSize = fftSize;
        SampleFrequency = (int)sampleFrequency;
        _device.SampleFrequencyMHz = sampleFrequency / 1000000;//  in  MHz
        _device.LNAGainDb = lnaGain;
        _device.VGAGainDb = vgaGain;
        _device.FilterBandwidthMHz = filterBandwidth;
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

    public void ReadFromHackRF(Action<double[]> onDataRecieved)
    {

        byte[] buffer;
        buffer = new byte[FftSize];
        //buffer = new byte[100000];

        try
        {
            for (int i = 0; i < 100; i++)
            {
                if (_dataStream == null || !_dataStream.CanRead) continue;
                _dataStream.Read(buffer, 0, buffer.Length); // reading interpolated IQ data from stream
                System.Console.WriteLine("demodulating...");
                var IQ_data = ConvertToIQ(buffer); // converting interpolated IQ data to complex array
                IQ_data = LPF1(IQ_data); // low pass filter to cutoff other frequencies
                var values = DemodFMsamples(IQ_data, 4200); // FM demodulator
                values = LPF2(values); // low pass filter to cutoff pilot tone, stereo fm subcarrier and RDS
                values = Decimate(values, 42); // changing sample frequency to 48kHz
                                               //var bytes = ConvertToBuffer(demod); // converting array of double to array of bytes
                                               //file_length += bytes.Length;
                                               //double[] values = IQ_data.Select(x => x.Real).ToArray(); // calculating power spectrum of the signal
                                               //FftSharp.FFT.Forward(values);
                                               //double[] values = FftSharp.FFT.Magnitude(IQ_data, positiveOnly: true); // calculating power spectrum of the signal
                                               //var window = new FftSharp.Windows.Welch();
                                               //window.ApplyInPlace(values);
                                               //FftSharp.Filter.LowPass(values, sampleRate: SampleFrequency, maxFrequency: 2000);
                                               //FftSharp.FFT.FftShift(values);

                onDataRecieved(values);

            }

        }
        catch (Exception ex)
        {
            System.Console.WriteLine($"Error while reading from HackRF: {ex.Message}");
            _device?.Reset(); // reset hackRF if something goes wrong
        }
    }

    static Complex[] ConvertToIQ(byte[] buffer)
    {
        Complex[] ret = new Complex[buffer.Length / 2];
        for (int i = 0; i < ret.Length; i++)
        {
            ret[i] = new Complex((sbyte)buffer[2 * i], (sbyte)buffer[2 * i + 1]);
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

}