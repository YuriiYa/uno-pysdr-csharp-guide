using System.Numerics;
using System.Runtime.Serialization;
using nethackrf;
using ScottPlot;
using HackRF.Namespace;

namespace uno_palyground.Presentation;

public static class TimeLapseHelper
{
    public static void PrintTime(Action action)
    {
        _ = PrintTime(() =>
         {
             action();
             return (false, false);
         });
    }

    public static (bool breakResult, bool continueResult) PrintTime(Func<(bool breakResult, bool continueResult)> action)
    {
        var watch = System.Diagnostics.Stopwatch.StartNew();
        var result = action();
        watch.Stop();

        var elapsedMs = watch.ElapsedMilliseconds;
        TimeSpan t = TimeSpan.FromMilliseconds(elapsedMs);
        string answer = string.Format("{0:D2}h:{1:D2}m:{2:D2}s:{3:D3}ms - time of {4}",
                        t.Hours,
                        t.Minutes,
                        t.Seconds,
                        t.Milliseconds,
                        action.Method.Name);

        Console.WriteLine(answer);
        return result;
    }
}
public sealed partial class MainPage : Page
{
    public MainPage()
    {
        this.InitializeComponent();

        // https://github.com/ScottPlot/ScottPlot/blob/main/src/ScottPlot5/ScottPlot5%20Demos/ScottPlot5%20WinForms%20Demo/Demos/SharedAxes.cs
        // WinUIPlot1.Plot.Axes.Link(WinUIPlotModulation.Plot, x: true, y: false);
        // WinUIPlotModulation.Plot.Axes.Link(WinUIPlot1.Plot, x: true, y: false);

        //var freqDomainLab = new FrequencyDomain(WinUIPlot1.Plot);
        //freqDomainLab.PopulateFrequency();
        // freqDomainLab.DrawSpectogram();
        // freqDomainLab.FFTSimulation();
        //var iqSampling = new IQSampling(WinUIPlot1.Plot);
        //iqSampling.CalculatePowerDpectralDensity();

        /* var d =new DigitalModulation(WinUIPlot1.Plot);
         d.QPSK();
         WinUIPlot1.Refresh();
 */
        //var noise = new Noise(WinUIPlot1.Plot, WinUIPlotModulation.Plot);
        // noise.GaussianNoise();
        //noise.ComplexNoise();

        /*var videoFileExample = new VideoFromFile(WinUIPlot1.Plot, WinUIPlotModulation.Plot);
        videoFileExample.DoShowVideoFile();*/

        // DecodePALVideoAsync(WinUIPlot1.Plot);


        var recording_time = 1; // seconds
        var center_freq = 1000 * 1000 * 100;  // 100 MHz
        var sample_rate = 1000 * 1000 * 10; //10 MHz
        var numSamples = (int)recording_time * sample_rate;//(int)recording_time * sample_rate; // 10 million samples for 1 second at 10 MHz
        var baseband_filter = 7.5e6;
        var lna_gain = 30; // 0 to 40 dB in 8 dB steps
        var vga_gain = 50; // 0 to 62 dB in 2 dB steps

        var configuration = new HackRFConfiguration(
                  recordingTimeInSec: recording_time,
                  deviceSerialNumber: "0000000000000000a09867dc386f51a3",
                  centerFrequency: center_freq,
                  sampleFrequency: sample_rate,
                  numSamples: numSamples,
                  lnaGain: lna_gain, // dB
                  vgaGain: vga_gain, // dB
                  filterBandwidth: baseband_filter, // MHz
                  fftSize: 2048 * 16
              );
        var d = new HackRF.Namespace.HackRF(WinUIPlot1.Plot, DispatcherQueue, WinUIPlotModulation.Plot, configuration, new HackRFInteraction());
        
        d.Start();
        WinUIPlot1.Refresh();
    }

    private Task DecodePALVideoAsync(Plot plot)
    {
        return Task.Factory.StartNew(() =>
           {
               try
               {
                   DecodePALVideo(plot);
               }
               catch (Exception ex)
               {
                   Console.WriteLine($"Decode error: {ex}");
               }
           }, TaskCreationOptions.LongRunning);
    }

    public void DecodePALVideo(Plot plot)
    {
        var filename = @"C:\projects\sdr\RFData\SDRSharp_20170122_171736Z_179100000Hz_IQ.wav";
        // for debian
        // var filename = @"/home/pi5-dos/SDRProj/video/SDRSharp_20170122_171736Z_179100000Hz_IQ.wav";

        if (!File.Exists(filename))
            throw new FileNotFoundException("IQ WAV file not found", filename);

        using var fs = new FileStream(filename, FileMode.Open, FileAccess.Read, FileShare.Read);

        var header = WavHelper.ReadWavHeader(fs);
        WavHelper.PrintWavHeader(header);

        // Initialize PAL decoder
        var palDecoder = new PALDecoder(plot, DispatcherQueue);

        // Decode PAL signal (assuming 10 MHz sample rate)
        // Ring buffer size: 2 frames worth (approx). Reserve after reading header so we know sample rate.
        int samplesPerLine = (int)(PALDecoder.PAL_LINE_DURATION * 10000000);
        int samplesPerFrame = samplesPerLine * PALDecoder.PAL_LINES_PER_FRAME;
        uno_palyground.PySDRGuide.IQWavReader.ConfigureRingBuffer(samplesPerFrame * 2);

        TimeLapseHelper.PrintTime(() =>
        {
            palDecoder.DecodePALSignal(sampleRate: 10000000, fs);
        });
    }

}