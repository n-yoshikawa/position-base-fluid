import matplotlib.animation as anm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Kernel functions
# reference: https://matthias-research.github.io/pages/publications/sca03.pdf
# slide by the PBF author: http://mmacklin.com/pbf_slides.pdf
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

fig = plt.figure()
ax = Axes3D(fig)
ax.set_xlim3d(-10, 10)
ax.set_ylim3d(-10, 10)
ax.set_zlim3d(0, 20)

N = 125  # Number of particle
t = 0    # time
dt = 0.1 # time step
# Initial positions and velocities
x = np.asarray([[i-2.0, j-2.0, k + 10.0] for i in range(5) for j in range(5) for k in range(5)])
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
for _ in range(5):
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
    # The only force is gravity
    force = gravity

    # apply forces
    v += dt * force

    # predict position
    x_pred = x + dt * v

    # find neighboring particles (should be optimized later)
    neighbor = []
    for i in range(N):
        distance = np.asarray([np.linalg.norm(x_pred[i] - x_pred[j]) for j in range(N)])
        neighbor.append(list(distance.argsort()[:8]))

    # the simulation loop (Algorithm 1)
    for _ in range(5):
        # calculate lambda_i
        lam = [0 for _ in range(N)]
        for i in range(N):
            rho = 0.0
            rho0 = rho_init[i]
            # the particle itself is included to calculate density
            for j in neighbor[i]:
                rho += W_poly6(x_pred[i] - x_pred[j], h)
            # Formula (1)
            C_i = rho / rho0 - 1.0
            Z = 0.0  # denominator for formula (11)
            # Formula (8)
            for k in neighbor[i]:
                if k == i: # if k = i
                    grad = np.zeros_like(x_pred[0])
                    for j in neighbor[i]:
                        if i != j: # Exclude i==j to avoid divide by zero (confirmed to the prof)
                            grad += dW_spiky(x_pred[i] - x_pred[j], h) / rho0
                    Z += np.linalg.norm(grad) ** 2.0
                else: # if k = j
                    grad = -dW_spiky(x_pred[i] - x_pred[k], h) / rho0
                    Z += np.linalg.norm(grad)**2.0
            # formula (11)
            epsilon = 10.0  # too big but this makes simulation stable
            lam[i] = -(C_i / (Z + epsilon))

        # calculate delta p_i
        delta_p = [0 for _ in range(N)]
        for i in range(N):
            rho0 = rho_init[i]
            # formula (12)
            for j in neighbor[i]:
                if i != j:  # Exclude i==j
                    delta_p[i] += (lam[i]+lam[j]) * dW_spiky(x_pred[i] - x_pred[j], h) / rho0
            # Collision detection and respose (14. of algorithm)
            surface_normal = np.asarray([0, 0, 1])  # floor normal
            # Implementation of 3.4 of Position Based Dynamics paper
            # https://www.cs.toronto.edu/~jacobson/seminar/mueller-et-al-2007.pdf
            if x_pred[i][2] < 0: # collision occurs
                if x[i][2] > 0:  # if ray x->p enters the floor
                    q = (x[i] * x[i][2] - x_pred[i] * x_pred[i][2]) / (x[i][2] - x_pred[i][2])
                else:  # if ray x->p lies completely inside the floor
                    q = x_pred[i] - surface_normal.dot(x_pred[i]) * surface_normal
                # Formula (7) of PBD paper where C(p) = (p-q) \dot n
                s = (x_pred[i] - q) * surface_normal
                # Formula (6) of PBD paper
                delta_p[i] -= s * surface_normal
        # update position (17. of the algorithm)
        for i in range(N):
            x_pred[i] += delta_p[i]

    for i in range(N):
        v[i] = (x_pred[i] - x[i]) / dt  # update velocity (21. of the algorithm)
        x[i] = x_pred[i]                # update position (23. of the algorithm)

    ax.scatter(x[:, 0], x[:, 1], x[:, 2], color='red')

ani = anm.FuncAnimation(fig, update, fargs = ('Initial Animation! ', 2.0), interval = 100)
plt.show()
