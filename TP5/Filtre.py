import numpy as np
import matplotlib.pyplot as plt

#from reconstruction import laminogram
sinogram = np.loadtxt("data/sinogram-patient.txt")
filtered_sinogram = []
for line in sinogram:
    fftline = np.fft.fft(line)
    fftshift = np.fft.fftshift(fftline)
    freq =  np.fft.fftfreq(len(line), 1)
    freqshift = np.fft.fftshift(freq)
    fftline = np.conj(np.abs(freqshift)*fftshift)
    lineshift = np.fft.ifftshift(fftline)
    newline = np.fft.fft(lineshift)

    filtered_sinogram.append(newline.real)

#print(filtered_sinogram)
plt.imshow(filtered_sinogram)
plt.show()
np.savetxt("data/sino_fftfiltre.txt",filtered_sinogram)
