import matplotlib.animation as anm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import time
import pbf
import argparse
import copy

t = 0    # time
dt = 0.01 # time step

args = None
v = None
x = None
gravity = None
# Calculate initial density
epsilon=5000
hit_bottom = False
ax = None
def init_params():
    global x, v, gravity, rho_init, args
    # Initial positions and velocities
    if args.particles == 1000:
        dim = 10
        x = np.linspace(-0.25, 0.25, num=dim)
        y = np.linspace(-0.25, 0.25, num=dim)
        z = np.linspace(0.0, 0.5, num=dim)
    elif args.particles == 8000:
        dim = 20
        x = np.linspace(-0.25*2, 0.25*2, num=dim)
        y = np.linspace(-0.25*2, 0.25*2, num=dim)
        z = np.linspace(0.05*2, 0.55*2, num=dim)

    X = []
    for x_i in x:
        for y_i in x:
            for z_i in z:
                X.append([x_i,y_i,z_i])
    N = len(X)
    
    x = np.array(X)
    v = np.asarray([[0., 0., 0.] for i in range(N)])
    gravity = np.asarray([[0., 0., -10.] for i in range(N)])
    # Calculate the density of the each particle at initial state
    rho_init = [5000 for _ in range(N)]

def update(i=None, fig_title=None, A=None):
    global t, x, v, gravity, hit_bottom, rho_init, args, ax
    if i != 0 and not args.save_pos:
        ax.cla()
    t += dt

    print("t = {:.4}".format(t))
    # change value of h for visual purposes:
    if hit_bottom or (not hit_bottom and np.min(x[:,2]) < 0.00001):
        h = 0.1
        hit_bottom = True
    elif not hit_bottom:
        h = 100000
    
    # set boundary constraints
    if args.particles == 1000:
        boundary = np.asarray([-0.5, 0.5, -0.5, 0.5, 1.25])
        if not args.save_pos:
            ax.set_xlim3d(-0.5, 0.5)
            ax.set_ylim3d(-0.5, 0.5)
            ax.set_zlim3d(0, 1)
    elif args.particles == 8000:
        boundary = np.asarray([-0.75, 0.75, -0.75, 0.75, 1.25])
        if not args.save_pos:
            ax.set_xlim3d(-0.75, 0.75)
            ax.set_ylim3d(-0.75, 0.75)
            ax.set_zlim3d(0, 1.25)
    if args.wave:
        if args.particles == 1000 and not args.save_pos:
            ax.set_xlim3d(-0.75, 0.75)
        elif args.particles == 8000 and not args.save_pos:
            ax.set_xlim3d(-1.0, 1.0)
        boundary[0] += 0.1 * np.sin(2.0 * np.pi * t)
    
    x, v = pbf.step(x, v, gravity, rho_init, dt, h, boundary, epsilon, args.pressure, args.vorticity, args.viscosity)
    if not args.save_pos:
        ax.scatter(x[:, 0], x[:, 1], x[:, 2], color='red')

def main():
    parser = argparse.ArgumentParser(description='Run particle fluid simulation.')
    parser.add_argument('--pressure', action='store_true', help='turn on artificial pressure kernel')
    parser.add_argument('--vorticity', action='store_true', help='turn on vorticity')
    parser.add_argument('--viscosity', action='store_true', help='turn on viscosity')
    parser.add_argument('--wave', action='store_true', help='generate waves (moving wall) instead of dropping particle cube')
    parser.add_argument('--side_view', action='store_true', help="change perspective of graph from bird's-eye view to side view")
    parser.add_argument('--particles', type=int, help='number of particles to simulate (1000 or 8000) --> default=1000', default=1000)
    parser.add_argument('--save_pos', action='store_true', help='save positions to numpy array called `positions.npy` --> can be used later for Blender renderings (note: will not show graph in order to run faster)')
    global args, ax
    args = parser.parse_args()
    assert args.particles in [1000, 8000], "Our system is only optimized for 1000 or 8000 particles"
    init_params()
    
    if args.save_pos:
        positions = []
        # run for 1000 time steps
        for _ in range(1000):
            update()
            positions.append(copy.deepcopy(x))

        with open('positions.npy', 'wb') as f:
            np.save(f, np.asarray(positions))
    else:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.set_xlim3d(-10, 10)
        ax.set_ylim3d(-10, 10)
        
        ax.set_zlim3d(0, 20)
        if args.side_view:
            ax.view_init(0,90)
        ani = anm.FuncAnimation(fig, update, fargs = ('Initial Animation! ', 2.0), interval = 100)
        plt.show()
if __name__ == "__main__":
    main()
