# python3 -m venv sdrenv
# ./sdrenv/Scripts/activate
# pip install matplotlib
# pip install scipy
# py .\wavFile.py

# Convert IQ samples in wav file to array of complex numbers and plot spectrum as function of time
# https://gist.github.com/asgaut/0b31455686ba87ffdac3d97028bec863

# Plot spectrum of wav file saved by SDR# as function of time (3D plot)
# Bandwidth = 32 kHz when input sample rate = 2.048 MHz (decimated by 64)
# Asgaut Eng/2016-10-02

import numpy as np
import scipy.io.wavfile as wav
import scipy.signal as signal
from scipy.fftpack import fft, fftfreq, fftshift
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fs, data = wav.read("C:\\projects\\sdr\\RFData\\SDRSharp_20170122_171736Z_179100000Hz_IQ.wav", mmap=True)
print("1 First few data:", data[0:10])

samples = data # np.fromfile('C:/projects/sdr/RFData/SDRSharp_20170122_171736Z_179100000Hz_IQ.wav', np.int16)

print("1 First few samples:", samples[0:10])
print("shape:", samples.shape)

print("np.iinfo(np.int16).max:", np.iinfo(np.int16).max)
samples = samples/np.iinfo(np.int16).max # convert to -1 to +1 (optional)
print("First few samples:", samples[0:10])
# Convert to complex IQ samples
iq_samples = samples[:, 0] + 1j * samples[:, 1]

print(iq_samples)

# Take one batch of samples equivalent to 8 MHz bandwidth
# For proper frequency resolution, use a power of 2 that covers the bandwidth well
batch_size = 8192  # Good FFT size for 8 MHz bandwidth
batch_samples = iq_samples[:batch_size]

# Sample rate from your WAV file
sample_rate = fs  # This should be 10 MHz from your wav.read()

# Apply window function to reduce spectral leakage
windowed_samples = batch_samples * np.hanning(len(batch_samples))
center_freq = 179100000  # Hz


# Apply FFT to convert to frequency domain
fft_result = np.fft.fft(windowed_samples)
fft_shifted = np.fft.fftshift(fft_result)  # Center the spectrum

# Calculate frequency bins relative to center frequency (179.1 MHz)

freqs = np.fft.fftshift(np.fft.fftfreq(batch_size, 1/sample_rate))
actual_freqs = center_freq + freqs  # Actual RF frequencies

# Convert to power spectral density (magnitude squared, in dB)
psd = 20 * np.log10(np.abs(fft_shifted) + 1e-12)  # Add small value to avoid log(0)

# Plot the frequency domain
plt.figure(figsize=(14, 8))
plt.plot(actual_freqs / 1e6, psd)  # Convert Hz to MHz for x-axis
plt.xlabel('Frequency (MHz)')
plt.ylabel('Power (dB)')
plt.title(f'PAL Broadcast Signal - Frequency Domain\nCenter: {center_freq/1e6} MHz, Bandwidth: 8 MHz')
plt.grid(True)

# Set x-axis limits to show the 8 MHz bandwidth
plt.xlim([(center_freq - 4e6)/1e6, (center_freq + 4e6)/1e6])

plt.show()

print(f"Batch size: {batch_size} samples")
print(f"Frequency resolution: {sample_rate/batch_size/1000:.1f} kHz")
print(f"Center frequency: {center_freq/1e6} MHz")
print(f"Frequency range: {(center_freq-sample_rate/2)/1e6:.1f} to {(center_freq+sample_rate/2)/1e6:.1f} MHz")

num_rows = len(windowed_samples) // batch_size # // is an integer division which rounds down
spectrogram = np.zeros((num_rows, batch_size))
for i in range(num_rows):
    spectrogram[i,:] = 10*np.log10(np.abs(np.fft.fftshift(np.fft.fft(windowed_samples[i*batch_size:(i+1)*batch_size])))**2)

plt.imshow(spectrogram, aspect='auto', extent = [sample_rate/-2/1e6, sample_rate/2/1e6, len(windowed_samples)/sample_rate, 0])
plt.xlabel("Frequency [MHz]")
plt.ylabel("Time [s]")
plt.show()