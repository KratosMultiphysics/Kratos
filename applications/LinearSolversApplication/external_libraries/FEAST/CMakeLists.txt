project(FEAST4 Fortran)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(OpenMP)
if (OPENMP_FOUND)
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif (OPENMP_FOUND)

# disable warnings
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")

if(${CMAKE_COMPILER_IS_GNUG77})
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w -O3 -ffree-line-length-none -ffixed-line-length-none -cpp -fPIC -DMKL -std=legacy")
  message("additional default options were set for gfortran")
  message("----------------------****************---------------------- CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")
else(${CMAKE_COMPILER_IS_GNUG77})
  message(FATAL_ERROR "Only the GNU Fortran compiler is supported.")
endif(${CMAKE_COMPILER_IS_GNUG77})

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/4.0/src/ )

set( CODE90
  ${CMAKE_CURRENT_SOURCE_DIR}/4.0/src/kernel/feast_tools.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/4.0/src/kernel/dzfeast.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/4.0/src/kernel/libnum.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/4.0/src/kernel/feast_aux.f90 )

set( CODE90_sparse
  ${CMAKE_CURRENT_SOURCE_DIR}/4.0/src/sparse/dzfeast_sparse.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/4.0/src/sparse/dzlsprim.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/4.0/src/sparse/sclsprim.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/4.0/src/sparse/dzfeast_pev_sparse.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/4.0/src/sparse/dzifeast_pev_sparse.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/4.0/src/sparse/dzifeast_sparse.f90 )

###############################################################
add_library( feast4 STATIC ${CODE90} ${CODE90_sparse} )
