import numpy as np
import sys
import mcubes
import os
from tqdm import tqdm 
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bins", type=int, default=30, help="number of bins/cells along each dimension")
    parser.add_argument("--path", type=str, default="positions_8000p_1000ts.npy", help="path to numpy array of particle positions at each time step (dimension of array should be num_timesteps x num_particles x (x,y,z))")
    args = parser.parse_args()

    binDim = args.bins
    data = np.load(args.path)
    num_particles = data.shape[1]
    dirname = f"marching_cubes_bin{binDim}_{num_particles}p_padded"
    os.makedirs(name=dirname, exist_ok=True)

    assert num_particles in [1000,8000], "System tested only for 1000 or 8000 particles."

    if num_particles == 8000:
        min_bound = -0.75
        max_bound = 1.25
    elif num_particles == 1000:
        min_bound = -0.5
        max_bound = 1
    
    # each unit has width (max_bound - min_bound)/binDim
    unit_width = (max_bound - min_bound)/binDim
    bins = np.arange(binDim)/binDim * (max_bound - min_bound) + min_bound
    bins = np.append([min_bound - unit_width],bins) 
    bins = np.append(bins, [max_bound, max_bound + unit_width]) 
    for timestep in tqdm(range(data.shape[0])):
        particles = data[timestep]
        inds = np.digitize(particles, bins)
        
        mesh = np.full((binDim+3, binDim+3, binDim+3), False, dtype=bool)
        for i in inds:
            x, y, z = i
            mesh[x][y][z] = True
        mesh = np.swapaxes(mesh, 1, 2)
        vertices, triangles = mcubes.marching_cubes(mesh, 0.5)
        mcubes.export_obj(vertices, triangles, f"{dirname}/bin{binDim}_step{timestep}.obj")

if __name__ == "__main__":
    main()
