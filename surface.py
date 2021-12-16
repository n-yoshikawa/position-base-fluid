import copy
import numpy as np
import pbf
from mayavi import mlab

t = 0    # time
dt = 0.01 # time step

# Initial positions and velocities
dim = 10
x = np.linspace(-0.25, 0.25, num=dim)
y = np.linspace(-0.25, 0.25, num=dim)
z = np.linspace(0.0, 0.2, num=4)
X = []
for x_i in x:
    for y_i in x:
        for z_i in z:
            X.append([x_i,y_i,z_i])
N = len(X)

x = np.array(X)
v = np.asarray([[0., 0., 0.] for i in range(N)])
gravity = np.asarray([[0., 0., -10.] for i in range(N)])

# Calculate initial density
boundary=0.5
epsilon=5000
# Calculate the density of the each particle at initial state
rho_init = [5000 for _ in range(N)]
hit_bottom = False
min_height_stack = np.min(z)
def update():
    global t, x, v, gravity, hit_bottom
    t += dt
    print("t = {:.2}".format(t))
    # change value of h for visual purposes:
    if hit_bottom or (not hit_bottom and np.min(x[:,2]) < 0.00001):
        h = 0.1
        hit_bottom = True
    elif not hit_bottom:
        h = 100000
    boundary = np.asarray([-0.25, 0.25, -0.25, 0.25])
    boundary[0] += 0.1 * np.sin(2.0 * np.pi * t)
    x, v = pbf.step(x, v, gravity, rho_init, dt, h, boundary, epsilon)

def kernel(s):
    return max(0.0, (1-s**2)**3)

def W_poly6(r_vect, h):
    r = np.linalg.norm(r_vect)
    if r <= h:
        return 315.0/(64.0 * np.pi * h**9.0) * (h**2.0 - r**2.0)**3.0
    else:
        return 0

if __name__ == "__main__": 
    positions = []
    X = np.linspace(-0.25, 0.25, num=10)
    Y = np.linspace(-0.25, 0.25, num=10)
    Z = np.linspace(0.0, 0.5, num=10)
    R = 0.1
    surface = np.zeros((10, 10, 10))
    for frame in range(500):
        update()
        positions.append(copy.deepcopy(x))
        for i, r_x in enumerate(X):
            for j, r_y, in enumerate(Y):
                for k, r_z, in enumerate(Z):
                    r = np.asarray([r_x, r_y, r_z])
                    density = np.sum([W_poly6(r - p, 0.1) for p in x])
                    surface[i, j, k] = density
        mlab.contour3d(surface, contours=3)
        mlab.savefig(f'surface{frame}.png')
        #mlab.show()
        mlab.clf()
    np.save(f, np.asarray(positions))
