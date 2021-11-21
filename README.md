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

# Video
[link](https://streamable.com/zm1kml)
