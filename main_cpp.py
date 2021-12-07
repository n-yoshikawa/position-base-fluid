import matplotlib.animation as anm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import pbf


N = 1728 # Number of particle
#N=125
#N = 64
t = 0    # time
dt = 0.005 # time step

# Initial positions and velocities

dim = int(np.cbrt(N))
x = np.asarray([(i/dim/2)-0.25 for i in range(dim)])
y = np.asarray([(i/dim/2)-0.25 for i in range(dim)])
z = np.asarray([(i/dim/4)+0.5 for i in range(dim)])
X = []
for x_i in x:
    for y_i in x:
        for z_i in z:
            X.append([x_i,y_i,z_i])

x = np.array(X)
v = np.asarray([[0., 0., 0.] for i in range(N)])
gravity = np.asarray([[0., 0., -10.] for i in range(N)])

# Calculate initial density
boundary=0.5
epsilon=5000
# Calculate the density of the each particle at initial state
rho_init = [8500 for _ in range(N)]
hit_bottom = False
min_height_stack = np.min(z)
def update(i, fig_title, A):
    global t, x, v, gravity, hit_bottom
    if i != 0:
        ax.cla()
    ax.set_xlim3d(-0.5, 0.5)
    ax.set_ylim3d(-0.5, 0.5)
    ax.set_zlim3d(0, 1)
    t += dt

    print("t = {:.2}".format(t))
    # change value of h for visual purposes:
    if hit_bottom or (not hit_bottom and np.min(x[:,2]) < 0.00001):
        h = 0.1
        hit_bottom = True
    elif not hit_bottom:
        h = 100000
    x, v = pbf.step(x, v, gravity, rho_init, dt, h, boundary, epsilon)
    ax.scatter(x[:, 0], x[:, 1], x[:, 2], color='red')

fig = plt.figure()
ax = Axes3D(fig)
ax.set_xlim3d(-10, 10)
ax.set_ylim3d(-10, 10)
ax.set_zlim3d(0, 20)

ani = anm.FuncAnimation(fig, update, fargs = ('Initial Animation! ', 2.0), interval = 100)
plt.show()
