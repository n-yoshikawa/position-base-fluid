# Position Based Fluid
Final project for Physics Based Animation

## Requirement

- numpy
- matplotlib

## Usage (Python ver)

```
python main.py
```

# C++ Implementation
Additional requirement: `pybind11`.

```
pip install pybind11
```

### Compile
```
c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) -I/usr/include/eigen3/ example.cpp -o example$(python3-config --extension-suffix)
```

### Run
```
python main_cpp.py
```

# Note for Blender
Blender rendering requires `positions.npy` which contains a numpy array whose shape is `(time_frame, num_particles, 3)`. See `main_blender.py` for example.

# Video
- [Rendered by matplotlib](https://streamable.com/zm1kml)
- [Rendered by Blender](https://streamable.com/4ni7qn)

<div style="width:100%;height:0px;position:relative;padding-bottom:56.250%;"><iframe src="https://streamable.com/e/4ni7qn" frameborder="0" width="100%" height="100%" allowfullscreen style="width:100%;height:100%;position:absolute;left:0px;top:0px;overflow:hidden;"></iframe></div>
