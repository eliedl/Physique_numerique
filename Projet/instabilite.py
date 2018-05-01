import numpy as np
import matplotlib.pyplot as plt

N = 1000
L = 1

x = np.linspace(0, L, N+1)
y = np.load('./carlos.npy')
print(y.shape)

plt.plot(x, y[:, 4], c='C0', alpha=0.4, label='20 ms')
plt.plot(x, y[:, 9], c='C0', alpha=0.6, label='40 ms')
plt.plot(x, y[:, 14], c='C0', alpha=0.8, label='60 ms')
plt.plot(x, y[:, 19], c='C0', alpha=1, label='80 ms')
plt.ylim(-0.0005, 0.001)
plt.legend(frameon=False, fontsize=14, loc='upper center', ncol=2)
plt.axis('off')
plt.show()

plt.plot(x, y[:, 140], c='C0', alpha=0.4, label='600 ms')
plt.plot(x, y[:, 145], c='C0', alpha=0.6, label='620 ms')
plt.plot(x, y[:, 150], c='C0', alpha=0.8, label='640 ms')
plt.plot(x, y[:, 155], c='C0', alpha=1, label='660 ms')
plt.ylim(-0.0005, 0.001)
plt.legend(loc='upper center',frameon=False, fontsize=14, ncol=2)
plt.axis('off')
plt.show()

