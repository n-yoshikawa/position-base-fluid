import matplotlib.animation as anm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import pbf


def W_poly6(r_vect, h):
    r = np.linalg.norm(r_vect)
    if r <= h:
        return 315.0/(64.0 * np.pi * h**9.0) * (h**2.0 - r**2.0)**3.0
    else:
        return 0

def dW_spiky(r_vect, h):
    r = np.linalg.norm(r_vect)
    if r <= h:
        return 45.0/(np.pi * h**6.0) * (h - r)**2.0 * (r_vect / r)
    else:
        return 0

N = 100  # Number of particle
t = 0    # time
dt = 0.1 # time step
# Initial positions and velocities
x = np.asarray([[0.5*r*np.cos(2.0*np.pi/10.0 * i), 0.5*r*np.sin(2.0*np.pi/10.0 * i), k + 10.0] for i in range(10) for r in range(1, 3) for k in range(5) ])
v = np.asarray([[0., 0., 0.] for i in range(N)])
gravity = np.asarray([[0., 0., -10.] for i in range(N)])

# Calculate initial density
# find neighboring particles
# Counting 8 neighbors is standard acoording to the prof
neighbor = []
for i in range(N):
    distance = np.asarray([np.linalg.norm(x[i] - x[j]) for j in range(N)])
    neighbor.append(list(distance.argsort()[:8]))


h = 100  # Parameter for kernel. I am not sure about the optimal value
# Calculate the density of the each particle at initial state
rho_init = [0.0 for _ in range(N)]  
for i in range(N):
    rho = 0.0
    for j in neighbor[i]:  
        rho_init[i] += W_poly6(x[i] - x[j], h)


def update(i, fig_title, A):
    global t, x, v, gravity
    if i != 0:
        ax.cla()
    ax.set_xlim3d(-10, 10)
    ax.set_ylim3d(-10, 10)
    ax.set_zlim3d(0, 20)
    t += dt

    print("t = {:.2}".format(t))

    x, v = pbf.step(x, v, gravity, rho_init, dt)

    ax.scatter(x[:, 0], x[:, 1], x[:, 2], color='red')

fig = plt.figure()
ax = Axes3D(fig)
ax.set_xlim3d(-10, 10)
ax.set_ylim3d(-10, 10)
ax.set_zlim3d(0, 20)

ani = anm.FuncAnimation(fig, update, fargs = ('Initial Animation! ', 2.0), interval = 100)
plt.show()
