import numpy as np
import matplotlib.pyplot as plt

def speed(x, C=1, L=1, d=0.1, sigma=0.3):
    return (C*x*(L-x)/L**2)*np.exp(-1*(x-d)**2/(2*sigma**2))

def propagation(L, N, tmax = 0.1, Nt = 1000, fps = 50, plot=False):

    # Parametres d'intégration
    #-------------------------
    v = 100           # Vitesse
    a = L/N           # Spacing
    h = tmax/Nt          # Time step
    epsilon = h/1000

    c = h*v**2/a**2


    # Initialisation variables d'intégration
    x = np.linspace(0, L, N+1)  # Ligne
    phi = np.zeros(N+1)          # Position initiale
    psi = speed(x)            # Vitesse

    phi_p = np.zeros(N + 1)
    psi_p = speed(x)

    t1 = 0.01
    t2 = 0.05
    t3 = 0.1

    tend = t3 + epsilon
    t = 0

    res = np.zeros((N+1, 1))
    count = 0
    while t<tmax:
        phi_p = phi + h*psi
        psi_p[1:N] = psi[1:N] + c*(phi[0:N-1] + phi[2:N+1]-2*phi[1:N])

        phi, phi_p = phi_p, phi
        psi, psi_p = psi_p, psi
        if t > count/fps:
            count += 1
            res = np.hstack((res, np.reshape(phi, (N+1, 1))))
            print("\r", t)

        t += h

        if plot:
            if abs(t - t1) < epsilon:
                plt.plot(phi, c='k')
                plt.axis('off')
                plt.show()

            if abs(t - t2) < epsilon:
                plt.plot(phi, c='k')
                plt.axis('off')
                plt.show()

            if abs(t-t3) < epsilon:
                plt.plot(phi, c='k')
                plt.axis('off')
                plt.show()

    return res

if __name__ == '__main__':

    N = 1000
    x = np.linspace(0, 1, N)
    phi = np.zeros(N+1)
    psi = speed(x)
    tmax = 0.090
    x = propagation(1, N, tmax = tmax, Nt = 5e6, fps = 500//tmax)
    np.save("carlos",x)


