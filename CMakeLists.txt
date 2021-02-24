cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(fmm LANGUAGES Fortran VERSION 1.4.2)
include(GNUInstallDirs)

option(ENABLE_PYTHON_INTERFACE "Builds Python bindings using f2py" OFF)

SET(SOURCES 
	src/fmm.F90 
	src/tensors.F90 
	src/tensors_recursive.F90
	src/tensors_damped_amoeba.F90 
	src/tensors_damped_erf.F90 
	src/tensors_damped_thole.F90 
	)

set(Fortran_FLAGS "-ffree-line-length-none -Ofast -mtune=native")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Fortran_FLAGS}")
add_library(fmm SHARED ${SOURCES})

if (ENABLE_PYTHON_INTERFACE)
	set(pymodule_name "octree")
	set(generated_module_file ${CMAKE_CURRENT_BINARY_DIR}/${pymodule_name}${PYTHON_EXTENSION_MODULE_SUFFIX})

	add_custom_target(${pymodule_name} ALL
	DEPENDS ${generated_module_file} fmm
	)
	add_custom_command(
	OUTPUT ${generated_module_file}
	COMMAND "f2py"
		-m ${pymodule_name}
		-c
		--f90flags=${Fortran_FLAGS}
		${CMAKE_SOURCE_DIR}/src/fmm.F90
		-lfmm
		-L${CMAKE_CURRENT_BINARY_DIR}
		only: field_fmm field_direct
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	)
endif()