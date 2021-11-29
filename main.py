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

N = 100 # Number of particle
t = 0    # time
dt = 0.1 # time step


# Initial positions and velocities
x = np.asarray([[0.5*r*np.cos(2.0*np.pi/10.0 * i), 0.5*r*np.sin(2.0*np.pi/10.0 * i), k + 10.0] for i in range(10) for k in range(5) for r in range(1, 3)])
#x = [ [2/5 * i, 2/5 * j, 1/20 * k + 10] for i in range(5) for j in range(5) for k in range(20)]
#print(np.cbrt(N))
x = np.array(x)
v = np.asarray([[0., 0., 0.] for i in range(N)])
gravity = np.asarray([[0., 0., -10.] for i in range(N)])
vorticity = np.asarray([[0., 0., 0.] for i in range(N)])
# Calculate initial density
# find neighboring particles
# Counting 8 neighbors is standard acoording to the prof
NUM_NEIGHBOR=8

## vectorized implementation of finding nearest neighbours
x_3D = x[:, np.newaxis]
distance = np.sum((x_3D - x)**2, axis=2)
neighbor = np.argsort(distance, axis=-1)[:,:NUM_NEIGHBOR]

h = 100  # Parameter for kernel. I am not sure about the optimal value
# Calculate the density of the each particle at initial state
# ***POTENTIALLY REVISE *** I think this should just be density of water - Marta
rho_init = [0.0 for _ in range(N)]  
for i in range(N):
    rho = 0.0
    for j in neighbor[i]:  
        rho_init[i] += W_poly6(x[i] - x[j], h)

def update(i, fig_title, A):
    global t, x, v, gravity, vorticity, NUM_NEIGHBOR
    if i != 0:
        ax.cla()
    ax.set_xlim3d(-2, 2)
    ax.set_ylim3d(-2, 2)
    ax.set_zlim3d(0, 20)
    t += dt

    print("t = {:.2}".format(t))
    # The only forces are gravity and vorticity
    force = gravity + vorticity

    # (2. of Algorithm 1 )apply forces
    v += dt * force

    # (3. of Algorithm 1) predict position
    x_pred = x + dt * v

    # (6. of Algorithm 1) find neighboring particles (should be optimized later)
    # updated with vectorized implementation
    x_pred_3D = x_pred[:, np.newaxis]
    distance = np.sum((x_pred_3D - x_pred)**2, axis=2)
    neighbor = np.argsort(distance, axis=-1)[:,:NUM_NEIGHBOR]
    
    # the simulation loop (Algorithm 1)
    for _ in range(5):
        # calculate lambda_i
        # (10. in Algorithm 1)
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
                else: # if k = j
                    grad = -1 * dW_spiky(x_pred[i] - x_pred[k], h) / rho0
                
                Z += np.linalg.norm(grad)**2.0
            # formula (11)
            epsilon = 1e-5
            lam[i] = C_i / (Z + epsilon) # the sign is flipped, but this yields better result

        # calculate delta p_i
        # (13. of Algorithm 1)
        delta_p = [0 for _ in range(N)]
        for i in range(N):
    
            rho0 = rho_init[i]
            # formula (12)
            for j in neighbor[i]:
                if i != j:  # Exclude i==j
                    # artificial pressure term:
                    s_corr = -0.1 * (W_poly6(x_pred[i]-x_pred[j],h)/W_poly6(0.001*h,h))**4
                    s_corr = 0
                    delta_p[i] += (lam[i]+lam[j] + s_corr) * dW_spiky(x_pred[i] - x_pred[j], h) / rho0
            # Collision detection between particles:
            '''
            for j in neighbor[i]:
                if i != j:
                    dist = np.linalg.norm(x_pred[i] - x_pred[j])
                    if dist <= 0.01: # collision
                        normal = (x_pred[i] - x_pred[j]) / dist
                        alpha_i = DAMPING*np.dot(normal.T, v[i])
                        alpha_j = DAMPING*np.dot(normal.T, v[j])
                        v[i] -= alpha_i*normal
                        v[j] -= alpha_j*normal
            '''
            # Collision detection and respose (14. of algorithm)
            surface_normal = np.asarray([0, 0, 1])  # floor normal
            wall1_normal = np.asarray([1,0,0])
            wall3_normal = np.asarray([-1,0,0])
            wall2_normal = np.asarray([0,1,0])
            wall4_normal = np.asarray([0,-1,0])
            ceiling_normal = np.asarray([0,0,-1])
            # Implementation of 3.4 of Position Based Dynamics paper
            # https://www.cs.toronto.edu/~jacobson/seminar/mueller-et-al-2007.pdf
            def compute_collision(dim, normal, x, x_pred,boundary):
                if (x[dim] > boundary and boundary <= 0) or (x[dim] < boundary and boundary > 0): # if ray x->p enters wall
                #if (x[dim] > boundary): # if ray x->p enters wall
                    #q = x_pred
                    #q[dim] = 0
                    q = (x * x[dim] - x_pred * x_pred[dim]) / (x[dim] - x_pred[dim]) # computes average in dimension of normal
                else: # if ray x->p lies completely inside wall
                    q = x_pred - normal.dot(x_pred) * normal # this will always be 0 in the dimensions that is being corrected ... normal.dot(x_pred) * normal isolates dimension being corrected; x_pred - normal.dot(x_pred) * normal sets dimension being corrected to 0
                #normal = ()
                # Formula (7) of PBD paper where C(p) = (p-q) \dot n
                s = (x_pred - q).dot(normal) 
                return s
            DAMPING=0.01
            if x_pred[i][2] < 0: # collision occurs on floor
                #s = compute_collision(2, surface_normal, x[i], x_pred[i],boundary=0)
                x_pred[i][2] = 0 # prevents it from getting stuck in corners
                #alpha = np.dot(v[i], surface_normal)
                #v[i] -= DAMPING*alpha * surface_normal
                # Formula (6) of PBD paper
                #delta_p[i] -= s * surface_normal
            
            if x_pred[i][0] < -1: # collision hits wall 1
                #s = compute_collision(0, wall1_normal, x[i], x_pred[i],boundary=-1)
                x_pred[i][0] = -1 
                #alpha = min([0, np.dot(v[i], wall1_normal)])
                #alpha=np.dot(v[i], wall1_normal)
                #v[i] -= DAMPING*alpha * wall1_normal
                #if s > 0.5:
                 #   print("left wall", s, x[i], x_pred[i])
                #delta_p[i] -= DAMPING * s * wall1_normal
            elif x_pred[i][0] > 1: # collision hits wall 3
                #s = compute_collision(0, wall3_normal, x[i], x_pred[i],boundary=1)
                x_pred[i][0] = 1 
                #alpha = min([0, np.dot(v[i], wall3_normal)])
                #alpha = np.dot(v[i], wall3_normal)
                #v[i] -= DAMPING*alpha * wall3_normal
                #if s > 0.5:
                 #   print("right wall", s, x[i], x_pred[i])
                #delta_p[i] -= DAMPING*s * wall3_normal
            if x_pred[i][1] < -1: # collision hits wall 2
                #s = compute_collision(1, wall2_normal, x[i], x_pred[i],boundary=-2)
                x_pred[i][1] = -1 
                #alpha = min([0, np.dot(v[i], wall2_normal)])
                #alpha = np.dot(v[i], wall2_normal)
                #v[i] -= DAMPING*alpha * wall2_normal
                #delta_p[i] -= s * wall2_normal
            if x_pred[i][1] > 1: # collision hits wall 4
                #s = compute_collision(1, wall4_normal, x[i], x_pred[i],boundary=2)
                x_pred[i][1] = 1 
                #alpha = min([0, np.dot(v[i], wall4_normal)])
                #alpha = np.dot(v[i], wall4_normal)
                #v[i] -= DAMPING*alpha * wall4_normal
                #delta_p[i] -= s * wall4_normal
            

        # update position (17. of the algorithm)
        for i in range(N):
            x_pred[i] += delta_p[i]

    for i in range(N):
        v[i] = (x_pred[i] - x[i]) / dt  # update velocity (21. of the algorithm)
    for i in range(N):
   ## VORTICITY (22. of the algorithm)
        omega_i = 0
        
        for j in neighbor[i]:
            if i != j:
                omega_i += (v[j] - v[i]) * dW_spiky(x_pred[i] - x_pred[j], h)
        eta=0
        for j in neighbor[i]:
            if i != j:
                eta += dW_spiky(x_pred[i] - x_pred[j], h) * np.linalg.norm(omega_i)
        
        # added small constant in case denominator = 0 (source: https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf section 5.1)
        N_vort = eta / (np.linalg.norm(eta)+10e-20)
        
        epsilon_vort = 1e-2
        #epsilon_vort = 0
        vorticity[i] = epsilon*(N_vort * omega_i)

        ## VISCOSITY (22. of the algorithm)
        c = 0.01
        viscosity = 0
        for j in neighbor[i]:
            viscosity += (v[j] - v[i]) * W_poly6(x_pred[i] - x_pred[j], h)
        v[i] += c * viscosity
    for i in range(N):
        x[i] = x_pred[i]                # update position (23. of the algorithm)

    ax.scatter(x[:, 0], x[:, 1], x[:, 2], color='blue')

fig = plt.figure()
ax = Axes3D(fig)
ax.set_xlim3d(-5, 5)
ax.set_ylim3d(-5, 5)
ax.set_zlim3d(0, 20)

ani = anm.FuncAnimation(fig, update, fargs = ('Initial Animation! ', 2.0), interval = 100)
plt.show()
