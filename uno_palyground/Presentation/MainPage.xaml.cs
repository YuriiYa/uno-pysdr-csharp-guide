using System.Numerics;
using System.Runtime.Serialization;
using nethackrf;
using ScottPlot;
using HackRF.Namespace;

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
        /*var d = new HackRF.Namespace.HackRF(WinUIPlot1.Plot, DispatcherQueue, WinUIPlotModulation.Plot);
        d.GetAndVisualizeHackRFData(new HackRFInteraction());
        WinUIPlot1.Refresh();*/
        //var noise = new Noise(WinUIPlot1.Plot, WinUIPlotModulation.Plot);
        // noise.GaussianNoise();
        //noise.ComplexNoise();

        /*var videoFileExample = new VideoFromFile(WinUIPlot1.Plot, WinUIPlotModulation.Plot);
        videoFileExample.DoShowVideoFile();*/

        DecodePALVideoAsync(WinUIPlot1.Plot);
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
        palDecoder.DecodePALSignal( sampleRate: 10000000, fs);
    }

}