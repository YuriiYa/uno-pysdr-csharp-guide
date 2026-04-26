using System.Numerics;
using uno_palyground.PySDRGuide;

namespace uno_palyground.PySDRGuide.Pipeline;

public class IqDemodulator
{
    public bool TryRead(Stream source, int samplesPerFrame, double[] frameBuffer)
    {
        int filled = 0;
        double sum = 0;
        while (filled < samplesPerFrame)
        {
            int need = samplesPerFrame - filled;
            if (source is FileStream fsrc)
            {
                var seg = IQWavReader.ReadIQIntoRingOptimized(fsrc, need);
                if (seg.IsEmpty) break;
                var first = seg.First;
                for (int i = 0; i < first.Length; i++)
                {
                    double mag = first[i].Magnitude;
                    frameBuffer[filled++] = mag;
                    sum += mag;
                }
                var second = seg.Second;
                for (int i = 0; i < second.Length && filled < samplesPerFrame; i++)
                {
                    double mag = second[i].Magnitude;
                    frameBuffer[filled++] = mag;
                    sum += mag;
                }
            }
            else
            {
                const int bytesPerSample = 2;
                int bytesNeeded = need * bytesPerSample;
                byte[] tmp = new byte[Math.Min(bytesNeeded, 8192)];
                int read = source.Read(tmp, 0, tmp.Length);
                if (read <= 0) break;
                int samplesRead = read / bytesPerSample;
                for (int s = 0; s < samplesRead && filled < samplesPerFrame; s++)
                {
                    double iVal = (sbyte)tmp[2 * s] / 128.0;
                    double qVal = (sbyte)tmp[2 * s + 1] / 128.0;
                    double mag = Math.Sqrt(iVal * iVal + qVal * qVal);
                    frameBuffer[filled++] = mag;
                    sum += mag;
                }
            }
        }
        if (filled < samplesPerFrame) return false;
        double dc = sum / samplesPerFrame;
        for (int i = 0; i < samplesPerFrame; i++) frameBuffer[i] -= dc;
        return true;
    }
}
