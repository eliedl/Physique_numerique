import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
from FTCS import speed


def tridiag(a, b, c, k1=-1, k2=0, k3=1):
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

def cond_init(x, x0, sigma=1e-10, k=5e10):
    return np.exp(-(x-x0)**2/(2*sigma**2))*np.exp(1j*k*x)

def wave(L, N, tmax = 0.1, Nt = 1000, fps = 50, plot=False):

    # Parametres d'intégration
    #-------------------------
    v = 100           # Vitesse
    a = L/N           # Spacing
    h = 1e-8        # Time step
    epsilon = h/1000

    c = h*v**2/(2*a**2)


    # Initialisation variables d'intégration
    x = np.linspace(0, L, N+1)  # Ligne
    phi = np.zeros(N+1)          # Position initiale
    psi = speed(x)            # Vitesse

    phi_p = np.zeros(N + 1)
    psi_p = speed(x)

    Phi = np.vstack((phi.reshape((N+1, 1)), psi.reshape((N+1, 1))))
    Phi_p = np.vstack((phi_p.reshape((N+1, 1)), psi_p.reshape((N+1, 1))))

    a1 = -2*np.ones(N+1)
    a2 = 1*np.ones(N)

    A_m = c*tridiag(a2, a1, a2)
    I = np.identity(N+1)
    I_m = 0.5*h*I

    mat1 = np.hstack((I, -I_m))
    mat2 = np.hstack((-A_m, I))
    M_m = np.vstack((mat1, mat2))

    mat3 = np.hstack((I, I_m))
    mat4 = np.hstack((A_m, I))
    M_p = np.vstack((mat3, mat4))


    M_m_I = np.linalg.inv(M_m)
    M = np.dot(M_m_I, M_p)
    print(M)


    t1 = 0.001
    t2 = 0.05
    t3 = 0.1

    tend = t3 + epsilon
    t = 0

    res = np.zeros((N+1, 1))
    count = 0
    while t<tmax:
        print(t, abs(t-t1))
        Phi_p = np.dot(M, Phi)

        t+= h

        Phi, Phi_p = Phi_p, Phi

        if abs(t - t1) < epsilon:
            plt.plot(Phi)
            plt.show()
        if abs(t - t2) < epsilon:
            plt.plot(Phi)
            plt.show()
        if abs(t - t3) < epsilon:
            plt.plot(Phi)
            plt.show()
    return res

def schro(N=1000, L=1e-8, h=1e-18):

    hbar = constants.value(u'Planck constant over 2 pi')
    m_e = 9.109e-31

    a = L / N
    c = h*hbar/(2*m_e*a**2)

    a1 = 1 + 1j*c*np.ones(N+1)
    a2 = -0.5j*c*np.ones(N)

    b1 = 1 - 1j*c*np.ones(N+1)
    b2 = 0.5j*c*np.ones(N)

    A = tridiag(a2, a1, a2)
    B = tridiag(b2, b1, b2)
    A_I = np.linalg.inv(A)

    M = np.dot(A_I, B)

    x = np.linspace(0, L, N+1)
    psi = cond_init(x, L/2)
    psi_p = np.zeros(N+1)


    #plt.plot(x[450:550], 1.6+psi.real[450:550], c='C3', label='0', alpha=0.2)

    t = 0
    t1 = 1e-17
    t2 = 2e-16
    t3 = 5e-16
    t4 = 1e-15
    t5 = 1.6e-15

    epsilon = h/1000

    res = np.zeros((N+1, 1))
    count = 0
    while t<1e-14:
        #print(t/1e-14)
        psi_p = np.dot(M, psi)

        if count%10 == 0:
            res = np.hstack((res, psi.reshape(N + 1, 1)))
            print("\r", t)


        t+=h
        count += 1
        psi, psi_p = psi_p, psi

        #if abs(t - t1) < epsilon:
        #    plt.plot(x, psi.real)
        #if abs(t - t2) < epsilon:
        #    plt.plot(x[505:610], 0.8+psi.real[505:610], c='C3', label='$2$', alpha=0.4)
        #if abs(t - t3) < epsilon:
        #    plt.plot(x[550:700], psi.real[550:700], c='C3', label='$5$', alpha=0.6)
        #if abs(t - t4) < epsilon:
        #    plt.plot(x[600:850], psi.real[600:850]-0.8, c='C3', label='$10$', alpha=0.8)
        #if abs(t - t5) < epsilon:
        #    plt.plot(x[700:1000], psi.real[700:1000]-1.6, c='C3', label='$16$', alpha=1)

    #plt.legend(frameon=False)
    #plt.axis('off')
    #plt.show()

    return res


if __name__ == '__main__':
    #N = 1000
    #x = np.linspace(0, 1, N)
    #tmax = 0.090
    #wave(1, N, tmax = tmax, Nt = 5e6, fps = 500//tmax)

    s = schro()
    np.save('Schro_alt', s)




