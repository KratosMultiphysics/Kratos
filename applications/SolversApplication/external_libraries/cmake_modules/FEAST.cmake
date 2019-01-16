project(FEAST Fortran)
include_directories( ${FEAST_DIR} )
include_directories( ${FEAST_DIR}/src )
include_directories( ${FEAST_DIR}/include )

find_package(OpenMP)
if (OPENMP_FOUND)
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif (OPENMP_FOUND)

# disable warnings
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")

# options for gfortran
if(${CMAKE_COMPILER_IS_GNUG77})
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w -O3 -ffree-line-length-none -ffixed-line-length-none -cpp -fPIC")
  message(STATUS "gfortran : CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")
endif(${CMAKE_COMPILER_IS_GNUG77})

# enable_language (Fortran)
set( CODE90_functions_wrapper ${FEAST_DIR}/src/kernel/f90_functions_wrapper.f90 )

set( CODE90 ${FEAST_DIR}/src/kernel/feast_tools.f90
            ${FEAST_DIR}/src/kernel/dzfeast.f90
            ${FEAST_DIR}/src/kernel/scfeast.f90
            ${FEAST_DIR}/src/kernel/feast_aux.f90 )


# set( CODE90_banded ${FEAST_DIR}/src/banded/dzfeast_banded.f90
#                     ${FEAST_DIR}/src/banded/scfeast_banded.f90  )
#
# set( CODE90_dense ${FEAST_DIR}/src/dense/dzfeast_dense.f90
#                     ${FEAST_DIR}/src/dense/scfeast_dense.f90  )

set( CODE90_sparse ${FEAST_DIR}/src/sparse/dzfeast_sparse.f90
                    ${FEAST_DIR}/src/sparse/scfeast_sparse.f90  )


###############################################################

add_definitions(-DINCLUDE_FEAST)

add_library( feast3 STATIC ${CODE90_functions_wrapper} ${CODE90} )
target_link_libraries(feast3 ${BLAS_LIBRARIES} )

# add_library( feast_dense STATIC ${CODE90_dense} )
# target_link_libraries(feast_dense ${BLAS_LIBRARIES} )
#
# add_library( feast_banded STATIC ${CODE90_banded} )
# target_link_libraries(feast_banded ${BLAS_LIBRARIES} )

# add_library( feast3_sparse STATIC ${CODE90_sparse} )
# target_link_libraries(feast3_sparse ${BLAS_LIBRARIES} )

#set(FEAST_LIBRARIES feast3 ${LAPACK_LIBRARIES} gfortran m)
set(FEAST_LIBRARIES feast3 ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} gfortran m)
