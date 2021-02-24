# fmm

## Building with CMake

```
mkdir build
cd build
cmake ..
make
```

## Building python interface with f2py

First, build the library with CMake.
In the build folder, run
```
f2py  --f90exec=gfortran  -m octree -c --opt="-ffree-line-length-none -Ofast -march=native -mtune=native -Wall" ../src/fmm.F90 CMakeFiles/fmm.dir/src/tensors*.o only: field_fmm field_direct
```
