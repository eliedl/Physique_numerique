import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

N = 1000
L = 1

x = np.linspace(0, L, N+1)
y = np.load('./carlos.npy')

#y = y[:, ::10]

#plt.plot(x, y[:, 0])
#plt.show()

print(y.shape)
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(0, 1), ylim=(-0.003, 0.003))
line, = ax.plot([], [], lw=2)

ax.set_xlabel('$x$')
ax.set_ylabel('Amplitude')

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    line.set_data(x, y[:, i+1])
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(y[0])-1, interval=20, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html

#plt.rcParams['animation.ffmpeg_path'] ='C:\\ffmpeg\\bin\\ffmpeg.exe'
#FFwriter = animation.FFMpegWriter(fps=30, extra_args=['-vcodec', 'libx264'])
#anim.save('animation_carlos_2.mp4', writer = FFwriter)

plt.show()