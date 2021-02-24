# fmm

## Building with CMake

```
mkdir build
cd build
cmake ..
make
```

## Building python interface with f2py

After building the library with CMake,
run f2py in the build folder
```
f2py  --f90exec=gfortran  -m octree -c --opt="-ffree-line-length-none -Ofast -march=native -mtune=native -Wall" ../src/fmm.F90 CMakeFiles/fmm.dir/src/tensors*.o only: field_fmm field_direct
```
