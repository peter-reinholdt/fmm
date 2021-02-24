cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(fmm LANGUAGES Fortran VERSION 1.4.2)
include(GNUInstallDirs)

SET(SOURCES 
	src/fmm.F90 
	src/tensors.F90 
	src/tensors_recursive.F90
	src/tensors_damped_amoeba.F90 
	src/tensors_damped_erf.F90 
	src/tensors_damped_thole.F90 
	)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -Ofast -mtune=native")
add_library(fmm SHARED ${SOURCES})
