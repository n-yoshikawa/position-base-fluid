import numpy as np
import sys
import mcubes
import os
from tqdm import tqdm 

binDim = 15
dirname = f"marching_cubes_bin{binDim}"
os.makedirs(name=dirname, exist_ok=True)
min_bound = -0.5
max_bound = 0.5
bins = np.arange(binDim)/binDim * (max_bound - min_bound) + min_bound
print(bins)
data = np.load("positions.npy")
for timestep in tqdm(range(1000)):
    particles = data[timestep]
    inds = np.digitize(particles, bins)
    
    mesh = np.full((binDim, binDim, binDim), False, dtype=bool)
    for i in inds:
        x, y, z = i -1
        mesh[x][y][z] = True
    mesh = np.swapaxes(mesh, 1, 2)
    vertices, triangles = mcubes.marching_cubes(mesh, 0.5)
    mcubes.export_obj(vertices, triangles, f"{dirname}/bin{binDim}_step{timestep}.obj")
