# python3 -m venv sdrenv
# ./sdrenv/Scripts/activate
# pip install scipy
# py .\wavFile.py

from scipy.io import wavfile
import numpy as np

# Replace 'audio.wav' with the path to your WAV file
try:
    sample_rate, data = wavfile.read('C:\\projects\\sdr\\RFData\\SDRSharp_20170122_171736Z_179100000Hz_IQ.wav')
    print(f"Sample Rate: {sample_rate} Hz")
    print(f"Data shape: {data.shape}")

    # The 'data' variable now holds the audio data as a NumPy array.
    #  - If the WAV file is mono (single channel), data will be a 1D array.
    #  - If the WAV file is stereo (2 channels), data will be a 2D array, where each row represents a sample,
    #    and each column represents a channel (left and right).

    # Example: Accessing the first sample of the first channel:
    if len(data.shape) == 1:
        first_sample = data[0]
    else:
        first_sample = data[0, 0]
    print(f"First sample: {first_sample}")

except FileNotFoundError:
    print("Error: The specified WAV file was not found.")
except Exception as e:
    print(f"An error occurred: {e}")