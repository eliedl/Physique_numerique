import numpy as np
import matplotlib.pyplot as plt
from scipy import constants


def tridiag(a, b, c, k1=-1, k2=0, k3=1):
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

def cond_init(x, x0, sigma=1e-10, k=5e10):
    return np.exp(-(x-x0)**2/(2*sigma**2))*np.exp(1j*k*x)


N = 1000
hbar = constants.value(u'Planck constant over 2 pi')
L = 1e-8
a = L/N
m_e = 9.109e-31
h = 1e-18


c = h*hbar/(2*m_e*a**2)

a1 = 1 + 1j*c*np.ones(N+1)
a2 = -0.5j*c*np.ones(N)

b1 = 1 + 1j*c*np.ones(N+1)
b2 = 0.5j*c*np.ones(N)


A = tridiag(a2, a1, a2)
B = tridiag(b2, b1, b2)
A_I = np.linalg.inv(A)


x = np.linspace(0, L, N+1)
psi = cond_init(x, L/2)
psi_p = np.zeros(N+1)


plt.plot(x[450:550], 1.6+psi.real[450:550], c='C3', label='0', alpha=0.2)

t = 0
t1 = 1e-17
t2 = 2e-16
t3 = 5e-16
t4 = 1e-15
t5 = 1.6e-15

epsilon = h/1000


while t<1e-14:
    psi_p = np.dot(A_I, psi)

    t+=h
    psi, psi_p = psi_p, psi

    #if abs(t - t1) < epsilon:
    #    plt.plot(x, psi.real)
    if abs(t - t2) < epsilon:
        plt.plot(x[505:610], 0.8+psi.real[505:610], c='C3', label='$2$', alpha=0.4)
    if abs(t - t3) < epsilon:
        plt.plot(x[550:700], psi.real[550:700], c='C3', label='$5$', alpha=0.6)
    if abs(t - t4) < epsilon:
        plt.plot(x[600:850], psi.real[600:850]-0.8, c='C3', label='$10$', alpha=0.8)
    if abs(t - t5) < epsilon:
        plt.plot(x[700:1000], psi.real[700:1000]-1.6, c='C3', label='$16$', alpha=1)

#plt.legend(frameon=False)
plt.axis('off')
plt.show()



#print(np.linalg.inv(A))





