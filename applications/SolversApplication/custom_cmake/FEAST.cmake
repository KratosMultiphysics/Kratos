# FEAST automatic download
# Variables:  FEAST_DIR

find_path(FEAST_FOUND NAMES feast.h HINTS "${FEAST_DIR}/include")

if(FEAST_FOUND)
  message(STATUS "Feast found in ${FEAST_DIR}")
else(FEAST_FOUND)
  message(STATUS "FEAST_DIR not defined: ${FEAST_DIR}")
  # Default version
  set(feast_install_dir "${EXTERNAL_LIBRARIES_DIR}")

  # Default FEAST_DIR
  set(FEAST_DIR "${feast_install_dir}/feast")
  message(STATUS "set FEAST_DIR: ${FEAST_DIR}")

  set(feast_file_name "feast-kratos")
  set(feast_packed_file "${feast_file_name}.zip")
  set(feast_packed_dir "${feast_install_dir}/${feast_packed_file}")

  if(NOT EXISTS ${feast_packed_dir})
    message(STATUS "${feast_packed_file} not found in ${feast_install_dir}")
    message(STATUS "Downloading ${feast_packed_file} from https://github.com/PFEM/feast/archive/kratos.zip to ${feast_install_dir} ...")
    file(DOWNLOAD https://github.com/PFEM/feast/archive/kratos.zip ${feast_packed_dir} SHOW_PROGRESS EXPECTED_MD5 8335ec9c8a52dd28eeefeda76c427b78)
    message(STATUS "${feast_packed_file} downloaded")
  endif(NOT EXISTS ${feast_packed_dir})

  find_path(FEAST_FOUND NAMES feast.h HINTS "${FEAST_DIR}/include")

  if((NOT EXISTS ${FEAST_DIR}) OR (NOT FEAST_FOUND))
    message(STATUS "Unpacking ${feast_packed_file} in ${feast_install_dir} ...")
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${feast_packed_file}
      WORKING_DIRECTORY ${feast_install_dir})
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${feast_install_dir}/${feast_file_name} ${FEAST_DIR}
      WORKING_DIRECTORY ${feast_install_dir})
    execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory ${feast_install_dir}/${feast_file_name}
      WORKING_DIRECTORY ${feast_install_dir})
  endif((NOT EXISTS ${FEAST_DIR}) OR (NOT FEAST_FOUND))
endif(FEAST_FOUND)


# Compile FEAST library
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
