# fmm

## Building with CMake

```
mkdir build
cd build
cmake ..
make
```

## Building python interface with f2py

Enable the python interface with CMake,
and build 
```
mkdir build
cd build
cmake .. -DENABLE_PYTHON_INTERFACE=ON
make
```

## Getting fields with the Fortran interface
Import the module in your source
```
    use octree, only: field_direct, field_fmm
```

### Direct
For standard direct fields, call
```
    subroutine field_direct(comm, coordinates, multipoles, exclusions, field_order, field)
```
where the arguments are
```
        integer, intent(in) :: comm               ! MPI comm, if you don't use MPI put in any dummy integer
        real(8), intent(in) :: coordinates(:, :)  ! Coordinates of particles, shape is (N x 3)
        real(8), intent(in) :: multipoles(:, :)   ! Multipoles of particles, shape is (N x multipole_dim); multipoles are packed in "alphabetical" order: q, x, y, z, xx, xy, xz, yy, yz, zz, ...
        integer, intent(in) :: exclusions(:, :)   ! exclusions(i, :) are the sites from which interactions of particle i are excluded. Padded with zeros.
        integer, intent(in) :: field_order        ! fields are calculated up to and including this order, e.g., 0 -> -potential, 1-> field, 2-> field gradient, and so on ...
```
The result is placed in packed form in 
```
        real(8), intent(out) :: field(size(coordinates, 1), (field_order + 1)*(field_order + 2)*(field_order + 3)/6)
```

### FMM
To compute fields with FMM, call
```
    subroutine field_direct(comm, coordinates, multipoles, exclusions, theta, ncrit, expansion_order, field_order, field)
```
arguments are the same as for the direct field, with the need to additionally specify
```
        real(8), intent(in) :: theta           ! angle opening criterion; between 0. and 1.0. A value of 0.5 is probably OK. Smaller values give higher accuracy, but is slower. 
        integer, intent(in) :: ncrit           ! max number of particles in a leaf node. 64 is probably fine. 
        integer, intent(in) :: expansion_order ! max expansion order. Higher is more accuracte. 5 is probably fine
```
### Damping
Both direct and FMM fields support damping. Choose between
- `AMOEBA`: amoeba style damping
- `THOLE `: thole style damping
- `ERF   `: another kind of damping

To activate damping, supply the damping type and damping factors with the optional arguments to `field_direct` or `field_fmm`:
```
        character(len=6), intent(in), optional :: damp_type ! "AMOEBA", "THOLE " or "ERF   "
        real(8), intent(in), optional :: damping_factors(:) ! damping factors, factorized such that a_ij = damping_factors(i) * damping_factors(j), probably something like sqrt(a) / alpha_i**(1.0/6.0)
```
