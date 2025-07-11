import numpy as np
from scipy.signal import find_peaks
from os.path import exists

from numpy import genfromtxt, max, abs

### inputs
S_thresh = 1e-5  # threshold for the populations std
spec_thres = 1e-3  # threshold for the spectrum

### checking files
# reading calculated data
input_file = 'autocorr_function.dat'
if not exists(input_file):
    print(f'ERROR: {input_file:s} not found!')
    exit(1)
S = genfromtxt(input_file).T

# reading reference data
input_file = 'autocorr_function_ref.dat'
if not exists(input_file):
    print(f'ERROR: {input_file:s} not found!')
    exit(1)
S_ref = genfromtxt(input_file).T

# calculating standard deviation
S_std = max(abs(S[1:] - S_ref[1:]))

### spectrum
t = S[0]
S = S[1] + 1j * S[2]

t = np.concatenate([-t[1:][::-1], t])
S = np.concatenate([np.conjugate(S[1:][::-1]), S])

# damping
S *= np.exp(-t ** 2 / (2 * 3500 ** 2))

# fft
spectrum = np.fft.ifft(S)
freq = 2 * np.pi * np.fft.fftfreq(len(t), d=t[1] - t[0])

# reoreder spectrum and frequency
freq = np.fft.fftshift(freq)
spectrum = np.fft.fftshift(spectrum)

# squared values
spectrum = np.abs(spectrum) ** 2

# locate all maxima in the spectrum
peaks, _ = find_peaks(np.abs(spectrum), height=0.00000001)

# reference energies got from imaginary time dynamics
ref = [-0.16074056, -0.14462691, -0.12923272, -0.11456857, -0.10064745, -0.08748475,
       -0.07509838, -0.06350902, -0.05274048, -0.0428203, -0.03378056, -0.02565916,
       -0.01850162, -0.01236403, -0.00731806, -0.00346031, -0.0009334]

spec_diff = np.max([abs(freq[peak]-ref[i]) for i, peak in enumerate(peaks)])

### deciding result
if (S_std < S_thresh) and (spec_diff < spec_thres):
    result = True
else:
    result = False

# printing result: True/False
print(result)


exit()

# get autocorrelation functions



sigma = 3500

data = np.genfromtxt('autocorr_function.dat').T
t = data[0]
S = data[1] + 1j * data[2]

t = np.concatenate([-t[1:][::-1], t])
S = np.concatenate([np.conjugate(S[1:][::-1]), S])

# damping
S *= np.exp(-t ** 2 / (2 * sigma ** 2))

spectrum = np.fft.ifft(S)
freq = 2 * np.pi * np.fft.fftfreq(len(t), d=t[1] - t[0])

# reoreder spectrum and frequency
freq = np.fft.fftshift(freq)
spectrum = np.fft.fftshift(spectrum)

# squared values
spectrum = np.abs(spectrum) ** 2

# locate all maxima in the spectrum
peaks, _ = find_peaks(np.abs(spectrum), height=0.00000001)
print("Peaks found at frequencies:", freq[peaks])
