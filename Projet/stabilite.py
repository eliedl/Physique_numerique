import numpy as np
import matplotlib.pyplot as plt

N = 1000
L = 1

x = np.linspace(0, L, N+1)
y = np.load('./wave_implicit.npy')

print(y.shape)

plt.plot(x, y[:, 11], c='C1', alpha=0.4, label='20 ms')
plt.plot(x, y[:, 22], c='C1', alpha=0.6, label='40 ms')
plt.plot(x, y[:, 33], c='C1', alpha=0.8, label='60 ms')
plt.plot(x, y[:, 44], c='C1', alpha=1, label='80 ms')
plt.ylim(-0.0005, 0.001)
plt.legend(frameon=False, fontsize=14, loc='upper center', ncol=2)
plt.axis('off')
plt.show()

plt.plot(x, y[:, 450], c='C1', alpha=0.4, label='600 ms')
plt.plot(x, y[:, 465], c='C1', alpha=0.6, label='620 ms')
plt.plot(x, y[:, 475], c='C1', alpha=0.8, label='640 ms')
plt.plot(x, y[:, 500], c='C1', alpha=1, label='660 ms')
plt.ylim(-0.0005, 0.001)
plt.legend(loc='upper center',frameon=False, fontsize=14, ncol=2)
plt.axis('off')
plt.show()

