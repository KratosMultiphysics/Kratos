cmake_minimum_required(VERSION 2.8)
project(TEST_MORSE_CMAKE_MODULES_FIND Fortran C CXX)

# location of Morse modules to get some specific macros
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../")
include(MorseInit)

# main variable: control the list of libraries to find thanks to find_package
# pay attention that package names must be given with capital letters
set(PACKAGES "" CACHE STRING "List of packages to find, ex: BLAS;STARPU;PASTIX")

# specific components to look for with packages
set(LAPACKE_COMPONENTS "" CACHE STRING "List of specific dependencies to look for with LAPACKE, ex: TMG")
set(QUARK_COMPONENTS "" CACHE STRING "List of specific dependencies to look for with QUARK, ex: HWLOC")
set(PASTIX_COMPONENTS "" CACHE STRING "List of specific dependencies to look for with PASTIX, ex: MPI;SEQ;STARPU;STARPU_CUDA;STARPU_FXT;SCOTCH;PTSCOTCH;METIS")
set(PETSC_COMPONENTS "" CACHE STRING "List of specific dependencies to look for with PETSc, ex: CXX|C")
set(FFTW_COMPONENTS "" CACHE STRING "List of specific dependencies to look for with FFTW, ex: MKL;ESSL;THREADS;OMP;SIMPLE;LONG;QUAD")
set(MUMPS_COMPONENTS "" CACHE STRING "List of specific dependencies to look for with MUMPS, ex: MPI;SEQ;SCOTCH;PTSCOTCH;METIS;PARMETIS;OPENMP")
set(CHAMELEON_COMPONENTS "" CACHE STRING "List of specific dependencies to look for with CHAMELEON, ex: STARPU;QUARK;CUDA;MPI;FXT")
set(STARPU_COMPONENTS "" CACHE STRING "List of specific dependencies to look for with STARPU, ex: HWLOC;CUDA;MPI;BLAS;MAGMA;FXT;SIMGRID")
set(PARSEC_COMPONENTS "" CACHE STRING "List of specific dependencies to look for with PARSEC, ex: HWLOC;CUDA;MPI;AYUDAME")

# to enable ctests
option(ENABLE_CTEST "Enable Testing: will test all supported packages" OFF)

foreach(_library ${PACKAGES})

  if (${_library}_COMPONENTS)
    find_package(${_library} COMPONENTS ${${_library}_COMPONENTS})
  else()
    find_package(${_library})
  endif()
  if (${_library}_FOUND)
    if (${_library}_LIBRARIES)
      message(STATUS "${_library}_LIBRARIES found: ${${_library}_LIBRARIES}")
    else()
      message(WARNING "${_library}_LIBRARIES not found: ${${_library}_LIBRARIES}")
    endif()
    if (${_library}_STATIC_LIBRARIES)
      message(STATUS "${_library}_STATIC_LIBRARIES found: ${${_library}_STATIC_LIBRARIES}")
    endif()
    if (${_library}_LIBRARIES_DEP)
      message(STATUS "${_library}_LIBRARIES_DEP found: ${${_library}_LIBRARIES_DEP}")
    endif()
    if (${_library}_LIBRARY_DIRS)
      message(STATUS "${_library}_LIBRARY_DIRS found: ${${_library}_LIBRARY_DIRS}")
    else()
      message(WARNING "${_library}_LIBRARY_DIRS not found: ${${_library}_LIBRARY_DIRS}")
    endif()
    if (${_library}_LIBRARY_DIRS_DEP)
      message(STATUS "${_library}_LIBRARY_DIRS_DEP found: ${${_library}_LIBRARY_DIRS_DEP}")
    endif()
    if (${_library}_INCLUDE_DIRS)
      message(STATUS "${_library}_INCLUDE_DIRS found: ${${_library}_INCLUDE_DIRS}")
    else()
      message(WARNING "${_library}_INCLUDE_DIRS not found: ${${_library}_INCLUDE_DIRS}")
    endif()
    if (${_library}_INCLUDE_DIRS_DEP)
      message(STATUS "${_library}_INCLUDE_DIRS_DEP found: ${${_library}_INCLUDE_DIRS_DEP}")
    endif()
    if (${_library}_CFLAGS_OTHER)
      message(STATUS "${_library}_CFLAGS_OTHER found: ${${_library}_CFLAGS_OTHER}")
    endif()
    if (${_library}_LDFLAGS_OTHER)
      message(STATUS "${_library}_LDFLAGS_OTHER found: ${${_library}_LDFLAGS_OTHER}")
    endif()
  else()
    message(FATAL_ERROR "${_library} NOT FOUND !!")
  endif()

endforeach()

# Add CTest rules
if (ENABLE_CTEST)

  enable_testing()
  include(CTest)

  set(CTEST_PACKAGES_LIST
      BLAS
      BLASEXT
      CBLAS
      CHAMELEON
      CPPCHECK
      EZTRACE
      FABULOUS
      FFTW
      FXT
      GTG
      HQR
      HWLOC
      HYPRE
      LAPACK
      LAPACKE
      LAPACKEXT
      METIS
      MPIEXT
      MUMPS
      PAMPA
      PAPI
      PARMETIS
      PARSEC
      PASTIX
      PETSc
      PTSCOTCH
      QUARK
      SCALAPACK
      SCOTCH
      SIMGRID
      SPM
      STARPU
      SUITESPARSE
      TMG)

  foreach(_package ${CTEST_PACKAGES_LIST})
    add_test(FIND${_package} cmake ${CMAKE_SOURCE_DIR} -DPACKAGES=${_package})
  endforeach()

endif()

###
### END CMakeLists.txt
###
