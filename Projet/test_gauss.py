import numpy as np
import matplotlib.pyplot as plt

def cond_init(x, x0, sigma=1e-9, k=5e10):
    return 1/np.sqrt(2*np.pi**2*sigma**2)*np.exp(-(x-x0)**2/(2*sigma**2))*np.exp(1j*k*x)

def cond_init_2D(x, y, x0, y0, sigma=5e-13, k=5e10):
    d = ((x-x0)**2 + (y-y0)**2)
    return np.exp(-d**2/(2*sigma**2))*np.exp(1j*k*x)

def fentes(x, y):
    temp = 0*(x**2 + y**2)
    temp[0:350, 500:550] = np.nan
    temp[450:550, 500:550] = np.nan
    temp[-351:-1, 500:550] = np.nan
    return  temp


L = 1e-5
N = 1000

x, y = np.meshgrid(np.linspace(0, L, N+1), np.linspace(0, L, N+1))

psi_2D = cond_init_2D(x, y, L/6, L/2)
psi_2D = psi_2D/np.max(psi_2D)

# Conditions frontières
f = fentes(x, y)


# Forçage de la fonction d'onde aux fentes
psi_2D[0:350, 500:550] = 0
psi_2D[450:550, 500:550] = 0
psi_2D[-351:-1, 500:550] = 0


plt.imshow(abs(psi_2D + f), extent=(0, L, 0, L), cmap='gray')
plt.ticklabel_format(style='sci', scilimits=(0,0))
plt.colorbar()
plt.show()