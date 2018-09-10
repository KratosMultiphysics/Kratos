project (KratosMultiphysics)
cmake_minimum_required (VERSION 2.8.6)



# Setting some policies
if(POLICY CMP0054)
  cmake_policy(SET CMP0054 OLD)
endif(POLICY CMP0054)

# Uses INTERFACE_LINK_LIBRARIES instead of LINK_INTERFACE_LIBRARIES
if(POLICY CMP0022)
  cmake_policy(SET CMP0022 NEW)
endif(POLICY CMP0022)



# Set here the version number **** only update upon tagging a release!
set (KratosMultiphysics_MAJOR_VERSION 6)
set (KratosMultiphysics_MINOR_VERSION 0)
set (KratosMultiphysics_PATCH_VERSION 0)

# Define custom compiler build types
SET( CMAKE_CONFIGURATION_TYPES Release RelWithDebInfo Debug FullDebug Custom)
SET( BASIC_DEBUG_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
#SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" )
#SET( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" )
SET( CMAKE_CXX_FLAGS_DEBUG "${BASIC_DEBUG_FLAGS} -DNDEBUG -DKRATOS_DEBUG" )
SET( CMAKE_CXX_FLAGS_FULLDEBUG "${BASIC_DEBUG_FLAGS} -DKRATOS_DEBUG" )
#SET( CMAKE_CXX_FLAGS_CUSTOM "${CMAKE_CXX_FLAGS_CUSTOM}" )

# Define internal CMake flags needed
SET( CMAKE_C_FLAGS_FULLDEBUG "${CMAKE_C_FLAGS_DEBUG}" )
SET( CMAKE_SHARED_LINKER_FLAGS_FULLDEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}" )
SET( CMAKE_EXE_LINKER_FLAGS_FULLDEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG}" )
SET( CMAKE_C_FLAGS_CUSTOM "" )
SET( CMAKE_CXX_FLAGS_CUSTOM "" )
SET( CMAKE_SHARED_LINKER_FLAGS_CUSTOM "" )
SET( CMAKE_EXE_LINKER_FLAGS_CUSTOM "" )
SET( CMAKE_MODULE_LINKER_FLAGS_FULLDEBUG "${CMAKE_MODULE_LINKER_FLAGS}")
SET( CMAKE_MODULE_LINKER_FLAGS_CUSTOM "${CMAKE_MODULE_LINKER_FLAGS}")

if(USE_COTIRE MATCHES ON)
    set(CMAKE_SKIP_RPATH ON)
endif(USE_COTIRE MATCHES ON)

# If build mode is not defined, assume Release
if(NOT CMAKE_BUILD_TYPE)
  message("Undefined CMAKE_BUILD_TYPE, building in 'Release' mode.")
  set(CMAKE_BUILD_TYPE Release)
else(NOT CMAKE_BUILD_TYPE)
  list(FIND CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE} FOUND_MODE)
  if(FOUND_MODE EQUAL -1)
    message("Unknown CMAKE_BUILD_TYPE '${CMAKE_BUILD_TYPE}', using 'Release'.")
    set(CMAKE_BUILD_TYPE Release)
  endif(FOUND_MODE EQUAL -1)
endif(NOT CMAKE_BUILD_TYPE)

# Get subversion data. This is done automagically by the cmakes
include (GenerateExportHeader)

# Search the SHA1 Associated with the commit in the HEAD
find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    OUTPUT_VARIABLE KratosMultiphysics_SHA1_NUMBER
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
else(GIT_FOUND)
  message("Git was not found on your system. SHA1 number will be set to 0.")
  set (KratosMultiphysics_SHA1_NUMBER 0)
endif(GIT_FOUND)

set (KratosMultiphysics_BUILD_TYPE ${CMAKE_BUILD_TYPE})
if( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set (KratosMultiphysics_BUILD_SUFFIX "")
else( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set (KratosMultiphysics_BUILD_SUFFIX "-${CMAKE_BUILD_TYPE}" )
endif( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

configure_file(
  "${PROJECT_SOURCE_DIR}/kratos_version.h.in"
  "${PROJECT_SOURCE_DIR}/kratos/includes/kratos_version.h"
)

# Set compiler flags
if(${CMAKE_COMPILER_IS_GNUCXX})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops -Wall -std=c++11")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.0)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsuggest-override")
  endif()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -funroll-loops -Wall")
  set(CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")
  message("additional default options were set for gcc")
  message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
  message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
endif(${CMAKE_COMPILER_IS_GNUCXX})

if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  if(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC  -funroll-loops -Wall -Wno-unused-local-typedef -Wno-unknown-pragmas  ")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops -Wall -Wno-unknown-pragmas  ")
    set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
    message("additional default options were set for clang compiler")
    message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
    message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
  endif(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
  if(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0)
      # pybind requires min. version 17, but we need at least 18:
      message( "DEPRECATED: detected compiler as Intel " ${CMAKE_CXX_COMPILER_VERSION} )
      message( FATAL_ERROR "Please use Version 18 or greater")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC  -funroll-loops -lpthread -wd654 -wd10010  ")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops -lpthread -wd654 -wd10010  ")
    set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
    message("additional default options were set for intel compiler")
    message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
    message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
  endif(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
else(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  if(${CMAKE_C_COMPILER} MATCHES "icc.*$")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC  -funroll-loops  ") #-lpthread -wd654 -wd10010  ")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops  ") #-lpthread -wd654 -wd10010  ")
    set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
    message("additional default options were set for intel compiler")
    message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
    message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
  endif(${CMAKE_C_COMPILER} MATCHES "icc.*$")
endif(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")



# Old non-compatible versions of VS
if(${MSVC80})
  message( "DEPRECATED: detected compiler as MSVC80")
  message( FATAL_ERROR "Please use VisualStudio 2015 or greater")
endif(${MSVC80})

if(${MSVC90})
  message( "DEPRECATED: detected compiler as MSVC90")
  message( FATAL_ERROR "Please use VisualStudio 2015 or greater")
endif(${MSVC90})

if(${MSVC10})
  message( "DEPRECATED: detected compiler as MSVC10")
  message( FATAL_ERROR "Please use VisualStudio 2015 or greater")
endif(${MSVC10})

if(${MSVC11})
  message( "DEPRECATED: detected compiler as MSVC11")
  message( FATAL_ERROR "Please use VisualStudio 2015 or greater")
endif(${MSVC11})

if(${MSVC12})
  message( "DEPRECATED: detected compiler as MSVC12")
  message( FATAL_ERROR "Please use VisualStudio 2015 or greater")
endif(${MSVC12})

# Common flags for MSVC
if(${MSVC})
  message("Detected compiler as MSVC")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W1 /bigobj /EHsc -DBOOST_ALL_NO_LIB -D_SCL_SECURE_NO_WARNINGS")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W1 /bigobj /EHsc -DBOOST_ALL_NO_LIB -D_SCL_SECURE_NO_WARNINGS")
  string( REPLACE "/W3" "" CMAKE_C_FLAGS ${CMAKE_C_FLAGS} )
  string( REPLACE "/W3" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} )
endif(${MSVC})

# Specific flags for different versions of MSVC
if(${MSVC14})
  # message("Adding additional flags for MSVC14")
  # Nothing needed right now
endif(${MSVC14})

# If no test policy enable by default
if(NOT DEFINED KRATOS_BUILD_TESTING)
	message("KRATOS_BUILD_TESTING not defined. Setting to ON")
	SET (KRATOS_BUILD_TESTING ON)
endif(NOT DEFINED KRATOS_BUILD_TESTING)

# If install prefix is set, install all the files.
if(DEFINED KRATOS_INSTALL_PREFIX)
	SET(INSTALL_TESTING_FILES ON)
	SET(INSTALL_PYTHON_FILES ON)
endif(DEFINED KRATOS_INSTALL_PREFIX)

# Set kratos specific module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")

if(USE_COTIRE MATCHES ON)
    include(cotire)
endif(USE_COTIRE MATCHES ON)

################### PYBIND11
include(pybind11Tools)

# Reset pybind11 config and remove -LTO on FullDebug to speedup linking time,
if(CMAKE_BUILD_TYPE MATCHES FullDebug)
    set(PYBIND11_LTO_CXX_FLAGS "" CACHE INTERNAL "")
    set(PYBIND11_LTO_LINKER_FLAGS "" CACHE INTERNAL "")
endif(CMAKE_BUILD_TYPE MATCHES FullDebug)

######################################################################################
######################################################################################
######################################################################################

## Find all dependencies

##*****************************
#find and include OpenMP if possible
message(KRATOS_EXCLUDE_OPENMP ${KRATOS_EXCLUDE_OPENMP})
if(${KRATOS_EXCLUDE_OPENMP} MATCHES ON)
	message("**********************************************************")
	message("OpenMP calculation was excluced artificially. Remove the var")
	message("KRATOS_EXCLUDE_OPENMP from the configuration script or set it to off to include it")
    message("**********************************************************")
    if(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -diag-disable 3180 ") # Disabling warning abt unrecognized omp pragmas
    endif(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
else(${KRATOS_EXCLUDE_OPENMP} MATCHES ON)
	message("OpenMP included")
	find_package(OpenMP)
	if (OPENMP_FOUND)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif (OPENMP_FOUND)
endif(${KRATOS_EXCLUDE_OPENMP} MATCHES ON)

##*****************************
# Finding and including BOOST library (version should not matter anymore)
find_package(Boost)

set(Boost_USE_STATIC_LIBS   OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_REALPATH ON)

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

message("Boost Include: ${Boost_INCLUDE_DIRS}")
message("Boost Linkdir: ${Boost_LIBRARY_DIRS}")

##*****************************

# Follow the symbolic links for the python lib only if needed
if(${INSTALL_EMBEDDED_PYTHON} MATCHES ON )
  foreach( l ${PYTHON_LIBRARIES})
    GET_FILENAME_COMPONENT(aaa ${l} REALPATH)
    set(PYTHON_LIBRARIES_REALPATH ${PYTHON_LIBRARIES_REALPATH} ${aaa})
  endforeach(l)

  set(PYTHON_LIBRARIES ${PYTHON_LIBRARIES_REALPATH} CACHE STRING "exact file of the PYTHON library" FORCE)
endif(${INSTALL_EMBEDDED_PYTHON} MATCHES ON )

##*****************************
# Finding blas and lapack if needed
set(BLAS_INCLUDE_NEEDED OFF)
if(${EXTERNAL_SOLVERS_APPLICATION} MATCHES ON )
  set(BLAS_INCLUDE_NEEDED ON)

  # If prescribed set the fortran mangling
  if(DEFINED ${KRATOS_SUPERLU_FORTRAN_MANGLING})
    add_definitions(${KRATOS_SUPERLU_FORTRAN_MANGLING})
  else(DEFINED ${KRATOS_SUPERLU_FORTRAN_MANGLING})
    add_definitions(-DADD_)
  endif(DEFINED ${KRATOS_SUPERLU_FORTRAN_MANGLING})
endif(${EXTERNAL_SOLVERS_APPLICATION} MATCHES ON )

if(${TRILINOS_APPLICATION} MATCHES ON )
  set(BLAS_INCLUDE_NEEDED ON)
endif(${TRILINOS_APPLICATION} MATCHES ON )

if(${USE_TETGEN_NONFREE_TPL} MATCHES ON )
  find_package(Tetgen REQUIRED)
endif(${USE_TETGEN_NONFREE_TPL} MATCHES ON )

# BLAS and LAPACK
if(${BLAS_INCLUDE_NEEDED} MATCHES ON )

  # Linux needs Fortran
  if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
 	  enable_language( Fortran )
    include(FortranCInterface)
 	  FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

  # Linking to the absolute value of the blas
  find_package(BLAS REQUIRED)
  include_directories(${BLAS_INCLUDE_DIRS})

  foreach( l ${BLAS_LIBRARIES})
    GET_FILENAME_COMPONENT(aaa ${l} REALPATH)
		set(BLAS_LIBRARIES_REALPATH ${BLAS_LIBRARIES_REALPATH} ${aaa})
  endforeach(l)

  set(BLAS_LIBRARIES ${BLAS_LIBRARIES_REALPATH} CACHE STRING "exact file of the blas library" FORCE)
  link_directories(${BLAS_LIBRARIES})

  set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}" CACHE STRING "Linker flags for shared libs" FORCE)
  message("linker flags with blas =" ${BLAS_LINKER_FLAGS})
  message("BLA_VENDOR = " ${BLA_VENDOR})

	# Linking to the absolute value of the lapack
  find_package(LAPACK REQUIRED)
  message("LAPACK = " ${LAPACK_LIBRARIES})
  foreach( l ${LAPACK_LIBRARIES})
    GET_FILENAME_COMPONENT(aux ${l} REALPATH)
    set(LAPACK_LIBRARIES_REALPATH ${LAPACK_LIBRARIES_REALPATH} ${aux})
  endforeach(l)

  set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES_REALPATH} CACHE STRING "exact file of the lapack library" FORCE)
  message("LAPACK modified = " ${LAPACK_LIBRARIES})

  include_directories(${LAPACK_INCLUDE_DIRS})
  link_directories(${LAPACK_LIBRARIES})

endif(${BLAS_INCLUDE_NEEDED} MATCHES ON )

##*****************************
#finding Mpi
if(${METIS_APPLICATION} MATCHES ON )
  set(MPI_NEEDED ON)
endif(${METIS_APPLICATION} MATCHES ON )
if(${TRILINOS_APPLICATION} MATCHES ON )
  set(MPI_NEEDED ON)
endif(${TRILINOS_APPLICATION} MATCHES ON )
if(${MPI_SEARCH_APPLICATION} MATCHES ON )
  set(MPI_NEEDED ON)
endif(${MPI_SEARCH_APPLICATION} MATCHES ON )

if(${MPI_NEEDED} MATCHES ON )
  find_package(MPI REQUIRED)
  add_definitions( -DKRATOS_USING_MPI )
  set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
  set(CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
  set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

  # Note: that this will be included in the compilation of the kratos core, so that it is not needed to include it somewhere else
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  # find_package(Boost COMPONENTS mpi_python REQUIRED)
endif(${MPI_NEEDED} MATCHES ON )

# Find zlib (for gidpost)
find_package(ZLIB)

if( ZLIB_FOUND )
  include_directories( ${ZLIB_INCLUDE_DIRS} )
else( ZLIB_FOUND )
  message( "-- Preparing local ZLIB compilation.")
  include_directories( ${CMAKE_SOURCE_DIR}/external_libraries/zlib )
  # Compile our own
  add_subdirectory(external_libraries/zlib)
  # we have to add this for zconf
  include_directories( ${CMAKE_BINARY_DIR}/external_libraries/zlib )
  set( ZLIB_LIBRARIES zlib )
endif( ZLIB_FOUND )

# including the AMatrix
if( AMATRIX_DIR )
    message("AMATRIX_DIR : ${AMATRIX_DIR}")
    include_directories ( "${AMATRIX_DIR}/include" )
    add_definitions(-DKRATOS_USE_AMATRIX)
else()
    message( "AMatrix not found. Please define the AMATRIX_DIR pointing to the AMatrix root directory.")
endif()

##echo user options
message( " ")
message("CMAKE_BUILD_TYPE ........ ${CMAKE_BUILD_TYPE}")
message( " ")
######################################################################################
######################################################################################
######################################################################################
#include internal dependencies
include_directories( ${CMAKE_SOURCE_DIR}/external_libraries )

# defines needed
add_definitions( -DKRATOS_PYTHON )

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  add_definitions( -fPIC )
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

# include subdirectories
add_subdirectory(external_libraries/gidpost)
add_subdirectory(kratos)
add_subdirectory(applications)

#compiling the triangle library
add_subdirectory(${CMAKE_SOURCE_DIR}/external_libraries/triangle)


# Embedded python - compiles "krun and packages the python lib needed"
if(${INSTALL_EMBEDDED_PYTHON} MATCHES ON )
  add_subdirectory(embedded_python)
endif(${INSTALL_EMBEDDED_PYTHON} MATCHES ON )

# Mpi
if(${MPI_NEEDED} MATCHES ON )
  add_subdirectory(external_libraries/mpi_python)
endif(${MPI_NEEDED} MATCHES ON )

# Set installation directory
if(DEFINED KRATOS_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX ${KRATOS_INSTALL_PREFIX} )
  message("-------------------------------- user defined install dir ${CMAKE_INSTALL_PREFIX}")
else(DEFINED KRATOS_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} )
  message("--------------------------------  standard install dir ${CMAKE_INSTALL_PREFIX}")
endif(DEFINED KRATOS_INSTALL_PREFIX)


################################################################################
# install core files for the KratosMultiphysics python module
# note that this does not install KratosLoader.py, as it is different for installation and packaging runs (see next code block)
install(FILES "${CMAKE_SOURCE_DIR}/kratos/python_interface/__init__.py" DESTINATION KratosMultiphysics )
install(FILES "${CMAKE_SOURCE_DIR}/kratos/python_interface/kratos_utilities.py" DESTINATION KratosMultiphysics )
install(FILES "${CMAKE_SOURCE_DIR}/kratos/python_interface/kratos_globals.py" DESTINATION KratosMultiphysics )
install(FILES "${CMAKE_SOURCE_DIR}/kratos/python_interface/application_importer.py" DESTINATION KratosMultiphysics )
install(FILES "${CMAKE_SOURCE_DIR}/kratos/python_interface/kratos_unittest.py" DESTINATION KratosMultiphysics RENAME KratosUnittest.py )

# Remove the tags in the event of multiple versions of boost being found
# in the same directory
# list(REMOVE_ITEM Boost_LIBRARIES "debug" "optimized")

# Install the libraries in the libs folder
install(FILES ${Boost_LIBRARIES} DESTINATION libs)
install(FILES ${EXTRA_INSTALL_LIBS} DESTINATION libs)

# Install blas and lapack
if(${BLAS_INCLUDE_NEEDED} MATCHES ON )
	message("installed blas = " ${BLAS_LIBRARIES})
	install(FILES ${BLAS_LIBRARIES} DESTINATION libs)
	message("installed lapack = " ${LAPACK_LIBRARIES})
	install(FILES ${LAPACK_LIBRARIES} DESTINATION libs)
endif(${BLAS_INCLUDE_NEEDED} MATCHES ON )

# Kratos Testing. Install everything except sources
if(${INSTALL_TESTING_FILES} MATCHES ON )
  install(DIRECTORY ${CMAKE_SOURCE_DIR}/kratos/tests/ DESTINATION ${CMAKE_INSTALL_PREFIX}/kratos/tests
    PATTERN "*.git" EXCLUDE
    PATTERN "*.c" EXCLUDE
    PATTERN "*.h" EXCLUDE
    PATTERN "*.cpp" EXCLUDE
    PATTERN "*.hpp" EXCLUDE
  )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

################################################################################
if(${INSTALL_PYTHON_FILES} MATCHES ON)
    #take care: do not indent the lines below
file(WRITE ${CMAKE_SOURCE_DIR}/packaging_aux/KratosLoader.py "import sys
import os.path
kratos_libs=os.path.abspath(os.path.join(os.path.dirname(__file__),'../libs'))
kratos_applications=os.path.abspath(os.path.join(os.path.dirname(__file__),'../applications'))
kratos_scripts=os.path.abspath(os.path.join(os.path.dirname(__file__),'../kratos/python_scripts'))
kratos_tests=os.path.abspath(os.path.join(os.path.dirname(__file__),'../kratos/tests'))
sys.path.append(kratos_libs)
sys.path.append(kratos_applications)
sys.path.append(kratos_scripts)
sys.path.append(kratos_tests)
")

    file(WRITE ${CMAKE_SOURCE_DIR}/packaging_aux/kratos.conf "${CMAKE_INSTALL_PREFIX}/libs")

    install(FILES ${CMAKE_SOURCE_DIR}/packaging_aux/KratosLoader.py DESTINATION KratosMultiphysics )
    install(FILES ${CMAKE_SOURCE_DIR}/packaging_aux/kratos.conf DESTINATION packaging )
    install(FILES ${CMAKE_SOURCE_DIR}/packaging_aux/script.py DESTINATION . )

    SET(CPACK_GENERATOR "DEB")
    SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Kratos Team") #required
    set(CPACK_PACKAGE_FILE_NAME  "Kratos_${KratosMultiphysics_VERSION_MAJOR}_${KratosMultiphysics_VERSION_MINOR}")
    #SET(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
    SET(CPACK_DEBIAN_PACKAGE_DEPENDS ${CPACK_DEBIAN_PACKAGE_DEPENDS} "libc6, libgcc1, libgfortran3, libgomp1, libopenmpi1.3, libparmetis3.1, libstdc++6, openmpi-bin, ${KRATOS_EXTRA_PACKAGE_DEPENDS} " )
#    SET(CPACK_DEBIAN_PACKAGE_DEPENDS ${CPACK_DEBIAN_PACKAGE_DEPENDS} " openmpi-bin, libboost-mpi-python-dev, ${KRATOS_EXTRA_PACKAGE_DEPENDS}" )
    Set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
    set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "${CMAKE_SOURCE_DIR}/packaging_aux/postinst;${CMAKE_SOURCE_DIR}/packaging_aux/prerm;")

    INCLUDE(CPack)
else(${INSTALL_PYTHON_FILES} MATCHES ON)
    #take care: do not indent the lines below
    file(WRITE ${CMAKE_INSTALL_PREFIX}/KratosMultiphysics/KratosLoader.py "import sys
kratos_libs=\"${CMAKE_INSTALL_PREFIX}/libs\"
kratos_applications=\"${CMAKE_SOURCE_DIR}/applications\"
kratos_scripts=\"${CMAKE_INSTALL_PREFIX}/kratos/python_scripts\"
kratos_tests=\"${CMAKE_INSTALL_PREFIX}/kratos/tests\"
sys.path.append(kratos_libs)
sys.path.append(kratos_applications)
sys.path.append(kratos_scripts)
sys.path.append(kratos_tests) ")
endif(${INSTALL_PYTHON_FILES} MATCHES ON)
