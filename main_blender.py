import copy
import numpy as np
import pbf

t = 0    # time
dt = 0.01 # time step

# Initial positions and velocities
dim = 10
x = np.linspace(-0.25, 0.25, num=dim)
y = np.linspace(-0.25, 0.25, num=dim)
z = np.linspace(0.4, 0.9, num=dim)
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
    x, v = pbf.step(x, v, gravity, rho_init, dt, h, boundary, epsilon)


if __name__ == "__main__": 
    positions = []
    for _ in range(100):
        update()
        positions.append(copy.deepcopy(x))

    with open('positions.npy', 'wb') as f:
        np.save(f, np.asarray(positions))
