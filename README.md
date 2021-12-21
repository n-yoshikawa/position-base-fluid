# Position Based Fluid
Final project for CSC2549: Physics-Based Animation\
University of Toronto, December 2021

## Requirements

- numpy
- matplotlib
- pybind11 (can be installed using `pip3 install pybind11`)

### Compile
```
c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) -I/usr/include/eigen3/ pbf.cpp -o pbf$(python3-config --extension-suffix)
```

Note that on MacOSX we had some issues with the compilation. This command worked for us:

```
c++ -O3 -Wall -shared -std=c++11 -fPIC -Wl,-undefined,dynamic_lookup $(python3 -m pybind11 --includes) -I/usr/include/Eigen/ pbf.cpp -o pbf$(python3-config --extension-suffix)
```
If you are on MacOSX, you might also see this error: 
`/usr/local/lib/python3.8/site-packages/pybind11/include/pybind11/eigen.h:30:10: fatal error: 'Eigen/Core' file not found #include <Eigen/Core>`\
We solved it by replacing `# include <Eigen/Core>` with `# include <eigen3/Eigen/Core>` in the `eigen.h` file.


### Run Fluid Simulation

The main script for running the simulation is:
```
python3 main.py
```
You can run it as is, but there are also a number of flags you can enter to change different parameters. These parameters are listed below and can also be listed by running `python3 main.py --help`

`--particles` : number of particles to run (our graphics are only scaled to 1000 or 8000 at this point; default is 1000 particles)\
`--pressure` : turns on articifial pressure kernel \  
`--vorticity` : turns on vorticity parameter\
`--viscosity` : turns on viscosity parameter \
`--wave` : move the wall in the `x`-dimension to generate waves\
`--side_view` : visualize the simulation from a side-view instead of bird's-eye view \
`--save_pos` : saves the particle positions for each timestep to a NumPy array called `positions.npy`, where the array has the following dimensions: `(timestep, num_particles, 3)`; this is useful for Blender renderings and isosurface generation\

For example, to run a simulation for 1000 particles with the pressure kernel and vorticity, you would run 
`python3 main.py --pressure --vorticity`

### Generate Surface Meshes using Marching Cubes

To generate an isosurface for the particles at each timestep, you can run the following script:

```
python3 marching_cubes.py
```
By default, it will take in a NumPy array called `positions.npy` and divide the cube into a 30x30x30 grid, but you can change this using the following command-line arguments: 

`--bins` : number of bins in each dimension (default is 30)\
`--path` : path to NumPy array \

This will save `.obj` files in a folder called `marching_cubes_bin{binDim}_{num_particles}p_padded` in the current working directory.

```
# Video
- [Rendered by matplotlib] (https://streamable.com/zm1kml)
- [1000 particles rendered by Blender] (https://streamable.com/b3sieg)
- [8000 particles rendered by Blender] (https://streamable.com/8q3pkf)
- [Mesh extracted using Marching Cubes] (https://streamable.com/gj5mix)
- [Mesh rendered by Blender Procedural Textures] (https://streamable.com/wu6xyc)
