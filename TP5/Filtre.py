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
'''
def filtre(Sinogram):

    sino_fft = np.zeros((720, 336))

    for i in range(0, 720):
        y = data_sino_patient[i]
        yfft_i = np.fft.fft(y)
        yfft_real = yfft_i.real
        yfft_shift_real = np.fft.fftshift(yfft_real)
        yfft_im = yfft_i.imag
        yfft_shift_im = np.fft.fftshift(yfft_im)
        f = np.fft.fftfreq(len(y), 1)
        f_shift = np.fft.fftshift(f)
        f_filtered_real = np.zeros(336)
        f_filtered_im = np.zeros(336)
        for j in range(0, len(f_filtered_real)):
            f_filtered_real[j] = yfft_shift_real[j] * np.abs(f_shift[j])
            f_filtered_im[j] = yfft_shift_im[j] * np.abs(f_shift[j])

        y_real = np.fft.ifftshift(f_filtered_real)
        y_im = np.fft.ifftshift(f_filtered_im)
        y_final = np.fft.fft(y_real + y_im * 1j)
        sino_fft[i] = np.array(y_final.real)

    return np.flip(sino_fft, )'''