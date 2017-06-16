project(FEAST Fortran) 
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(OpenMP)
if (OPENMP_FOUND)
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif (OPENMP_FOUND)

# disable warnings
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")

if(${CMAKE_COMPILER_IS_GNUG77})
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w -O3 -ffree-line-length-none -ffixed-line-length-none -cpp -fPIC")
  message("additional default options were set for gfortran")
  message("----------------------****************---------------------- CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")
endif(${CMAKE_COMPILER_IS_GNUG77})
##TODO: compiler settings for other compilers!

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/3.0/src/ )


# enable_language (Fortran)

set( CODE90_functions_wrapper ${CMAKE_CURRENT_SOURCE_DIR}/3.0/src/kernel/f90_functions_wrapper.f90 )

set( CODE90 ${CMAKE_CURRENT_SOURCE_DIR}/3.0/src/kernel/feast_tools.f90 
            ${CMAKE_CURRENT_SOURCE_DIR}/3.0/src/kernel/dzfeast.f90 
            ${CMAKE_CURRENT_SOURCE_DIR}/3.0/src/kernel/scfeast.f90 
            ${CMAKE_CURRENT_SOURCE_DIR}/3.0/src/kernel/feast_aux.f90 )

                                              
# set( CODE90_banded ${CMAKE_CURRENT_SOURCE_DIR}/3.0/src/banded/dzfeast_banded.f90 
#                     ${CMAKE_CURRENT_SOURCE_DIR}/3.0/src/banded/scfeast_banded.f90  )
# 
# set( CODE90_dense ${CMAKE_CURRENT_SOURCE_DIR}/3.0/src/dense/dzfeast_dense.f90 
#                     ${CMAKE_CURRENT_SOURCE_DIR}/3.0/src/dense/scfeast_dense.f90  )

set( CODE90_sparse ${CMAKE_CURRENT_SOURCE_DIR}/3.0/src/sparse/dzfeast_sparse.f90 
                    ${CMAKE_CURRENT_SOURCE_DIR}/3.0/src/sparse/scfeast_sparse.f90  )

                    
###############################################################
add_library( feast STATIC ${CODE90_functions_wrapper} ${CODE90} )
target_link_libraries(feast ${BLAS_LIBRARIES} )

# add_library( feast_dense STATIC ${CODE90_dense} )
# target_link_libraries(feast_dense ${BLAS_LIBRARIES} )
# 
# add_library( feast_banded STATIC ${CODE90_banded} )
# target_link_libraries(feast_banded ${BLAS_LIBRARIES} )

add_library( feast_sparse STATIC ${CODE90_sparse} )
target_link_libraries(feast_sparse ${BLAS_LIBRARIES} )


