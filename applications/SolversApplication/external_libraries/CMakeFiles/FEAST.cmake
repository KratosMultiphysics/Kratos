project(FEAST Fortran)
include_directories( ${FEAST_DIR} )
include_directories( ${FEAST_DIR}/src )
include_directories( ${FEAST_DIR}/include )

find_package(OpenMP)
if (OPENMP_FOUND)
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif (OPENMP_FOUND)

# disable warnings
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w -fPIC")
message(STATUS "CMAKE_Fortran_COMPILER : ${CMAKE_Fortran_COMPILER}")
message(STATUS "CMAKE_COMPILER_IS_GNUG77 : ${CMAKE_COMPILER_IS_GNUG77}")

# options for gfortran
if(${CMAKE_COMPILER_IS_GNUG77})
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -m64 -ffree-line-length-none -ffixed-line-length-none -cpp")
  message(STATUS "CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")
endif(${CMAKE_COMPILER_IS_GNUG77})

# enable_language (Fortran)
set( CODE90_functions_wrapper ${FEAST_DIR}/src/kernel/f90_functions_wrapper.f90 )

set( CODE90 ${FEAST_DIR}/src/kernel/feast_tools.f90
            ${FEAST_DIR}/src/kernel/dzfeast.f90
            ${FEAST_DIR}/src/kernel/scfeast.f90
            ${FEAST_DIR}/src/kernel/feast_aux.f90 )


# set( CODE90_banded ${FEAST_DIR}/src/banded/dzfeast_banded.f90
#                     ${FEAST_DIR}/src/banded/scfeast_banded.f90 )
# set( CODE90_dense ${FEAST_DIR}/src/dense/dzfeast_dense.f90
#                     ${FEAST_DIR}/src/dense/scfeast_dense.f90 )
#set( CODE90_sparse ${FEAST_DIR}/src/sparse/dzfeast_sparse.f90
#                    ${FEAST_DIR}/src/sparse/scfeast_sparse.f90 )


###############################################################

if(NOT BLAS_FOUND)
  find_package(BLAS)
endif(NOT BLAS_FOUND)

if(NOT LAPACK_FOUND)
  find_package(LAPACK)
endif(NOT LAPACK_FOUND)

message(STATUS "FEAST blas: ${BLAS_LIBRARIES}")
message(STATUS "FEAST lapack: ${LAPACK_LIBRARIES}")

add_library( external_feast STATIC ${CODE90_functions_wrapper} ${CODE90} )
target_link_libraries(external_feast ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
#target_link_libraries(external_feast ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} gfortran m)

# add_library( external_feast_dense STATIC ${CODE90_dense} )
# target_link_libraries(external_feast_dense ${BLAS_LIBRARIES} )
#
# add_library( external_feast_banded STATIC ${CODE90_banded} )
# target_link_libraries(external_feast_banded ${BLAS_LIBRARIES} )

# add_library( external_feast_sparse STATIC ${CODE90_sparse} )
# target_link_libraries(external_feast_sparse ${BLAS_LIBRARIES} )
