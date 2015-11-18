cmake_minimum_required (VERSION 2.8.6)
project (KratosMultiphysics)

#set here the version number **** only update upon tagging a release!
set (KratosMultiphysics_VERSION_MAJOR 3)
set (KratosMultiphysics_VERSION_MINOR 3)

#get subversion data. This is done automagically by the cmakes
include (GenerateExportHeader)

find_package(Subversion QUIET)
if(Subversion_FOUND)
  get_filename_component(ABS_PROJECT_SOURCE_DIR ${PROJECT_SOURCE_DIR} REALPATH)
  include(FindSubversion)
  Subversion_WC_INFO(${ABS_PROJECT_SOURCE_DIR} KRATOS_SVN_REVISION)
  configure_file(
      "${PROJECT_SOURCE_DIR}/kratos_version.h.in"
      "${PROJECT_SOURCE_DIR}/kratos/includes/kratos_version.h"
      )
else(Subversion_FOUND)
  # if svn is not found use 0 as a default placeholder
  set (KRATOS_SVN_REVISION_WC_REVISION 0)
  configure_file(
    "${PROJECT_SOURCE_DIR}/kratos_version.h.in"
    "${PROJECT_SOURCE_DIR}/kratos/includes/kratos_version.h"
    )
endif(Subversion_FOUND) 

#set compiler flags
if(${CMAKE_COMPILER_IS_GNUCXX}) 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops -ffast-math -Wall -Wno-strict-aliasing ")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -funroll-loops -ffast-math -Wall -Wno-strict-aliasing")
  set(CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")
  message("additional default options were set for gcc")
  message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
  message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
endif(${CMAKE_COMPILER_IS_GNUCXX}) 
if(${CMAKE_C_COMPILER} MATCHES "icc.*$") 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC  -funroll-loops  -lpthread -wd654 -wd10010  ")
  set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops  -lpthread -wd654 -wd10010  ")
  set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
  message("additional default options were set for intel compiler")
  message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
  message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
endif(${CMAKE_C_COMPILER} MATCHES "icc.*$") 

if(${MSVC80})
  message("detected compiler as MSVC80")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES")
endif(${MSVC80}) 

if(${MSVC90})
  message("detected compiler as MSVC90")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES")
endif(${MSVC90}) 

if(${MSVC10})
  message("detected compiler as MSVC10")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES")
endif(${MSVC10}) 

if(${MSVC11})
  message("detected compiler as MSVC11")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS")
endif(${MSVC11}) 

if(${MSVC12})
  message("detected compiler as MSVC12")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS")
endif(${MSVC12}) 

if(${MSVC14})
  message("detected compiler as MSVC14")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS")
endif(${MSVC14}) 

#tell the linker to give an error if undefined functions are found
# set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")

#set kratos specific module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")

######################################################################################
######################################################################################
######################################################################################

# if(DEFINED KRATOS_EXTRA_LINK_LIBRARIES)
#     link_directories(${KRATOS_EXTRA_LINK_PATH})
#     find_library(EXTRA_LINK_LIBS ${KRATOS_EXTRA_LINK_LIBRARIES} PATH ${KRATOS_EXTRA_LINK_PATH} )
# endif(DEFINED KRATOS_EXTRA_LINK_LIBRARIES)



##find all dependencies

##*****************************
#find and include OpenMP if possible
message(KRATOS_EXCLUDE_OPENMP ${KRATOS_EXCLUDE_OPENMP})
if(${KRATOS_EXCLUDE_OPENMP} MATCHES ON)
	message("**********************************************************")
	message("OpenMP calculation was excluced artificially. Remove the var")
	message("KRATOS_EXCLUDE_OPENMP from the configuration script or set it to off to include it")
	message("**********************************************************")
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
#finding and including BOOST library
set(Boost_USE_STATIC_LIBS   OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_REALPATH ON)
find_package(Boost COMPONENTS python REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

##*****************************
#finding and including PYTHON
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
link_directories(${PYTHON_LIBRARIES})



#follow the symbolic links for the python lib only if needed
if(${INSTALL_EMBEDDED_PYTHON} MATCHES ON )
    foreach( l ${PYTHON_LIBRARIES})
    GET_FILENAME_COMPONENT(aaa ${l} REALPATH)
      set(PYTHON_LIBRARIES_REALPATH ${PYTHON_LIBRARIES_REALPATH} ${aaa})
    endforeach(l)
    set(PYTHON_LIBRARIES ${PYTHON_LIBRARIES_REALPATH} CACHE STRING "exact file of the PYTHON library" FORCE)
endif(${INSTALL_EMBEDDED_PYTHON} MATCHES ON )

##*****************************
#finding blas and lapack if needed
set(BLAS_INCLUDE_NEEDED OFF)
if(${EXTERNAL_SOLVERS_APPLICATION} MATCHES ON )
    set(BLAS_INCLUDE_NEEDED ON)

    ##if prescribed set the fortran mangling
    if(DEFINED ${KRATOS_SUPERLU_FORTRAN_MANGLING})
      add_definitions(${KRATOS_SUPERLU_FORTRAN_MANGLING})
    else(DEFINED ${KRATOS_SUPERLU_FORTRAN_MANGLING})
      add_definitions(-DADD_)
    endif(DEFINED ${KRATOS_SUPERLU_FORTRAN_MANGLING})
endif(${EXTERNAL_SOLVERS_APPLICATION} MATCHES ON )
if(${TRILINOS_APPLICATION} MATCHES ON )
    set(BLAS_INCLUDE_NEEDED ON)
endif(${TRILINOS_APPLICATION} MATCHES ON )

if(${BLAS_INCLUDE_NEEDED} MATCHES ON )
     if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
 	enable_language( Fortran )
        include(FortranCInterface) 
 	FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")
     endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
   
    #linking to the absolute value of the blas
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



	#linking to the absolute value of the lapack
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
set(MPI_NEEDED OFF)
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

    #note that this will be included in the compilation of the kratos core, so that it is not needed to include it somewhere else
    include_directories(${MPI_INCLUDE_PATH})
#    find_package(Boost COMPONENTS mpi_python REQUIRED)
endif(${MPI_NEEDED} MATCHES ON )


##echo user options
message( " ")
message("CMAKE_BUILD_TYPE ........ ${CMAKE_BUILD_TYPE}")
message( " ")
######################################################################################
######################################################################################
######################################################################################
#include internal dependencies
include_directories( ${CMAKE_SOURCE_DIR}/external_libraries )
include_directories( ${CMAKE_SOURCE_DIR}/external_libraries/omptl )
include_directories( ${CMAKE_SOURCE_DIR}/external_libraries/zlib )
# we have to add this for zconf
include_directories( ${CMAKE_BINARY_DIR}/external_libraries/zlib )
# include_directories( ${CMAKE_SOURCE_DIR}/external_libraries/gidpost )

#defines needed
add_definitions( -DKRATOS_PYTHON )
add_definitions( -DNDEBUG )
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    add_definitions( -fPIC )
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

#include subdirectories
add_subdirectory(external_libraries/zlib)
add_subdirectory(external_libraries/gidpost)
add_subdirectory(kratos)
add_subdirectory(applications)

#embedded python - compiles "krun and packages the python lib needed"
if(${INSTALL_EMBEDDED_PYTHON} MATCHES ON )
  add_subdirectory(embedded_python)
endif(${INSTALL_EMBEDDED_PYTHON} MATCHES ON )



if(${MPI_NEEDED} MATCHES ON )
   add_subdirectory(external_libraries/mpi_python)
endif(${MPI_NEEDED} MATCHES ON )

#set installation directory
if(DEFINED KRATOS_INSTALL_PREFIX)
   set(CMAKE_INSTALL_PREFIX ${KRATOS_INSTALL_PREFIX} )	
   message("-------------------------------- user defined install dir ${CMAKE_INSTALL_PREFIX}")
else(DEFINED KRATOS_INSTALL_PREFIX)
   set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} )
   message("--------------------------------  standard install dir ${CMAKE_INSTALL_PREFIX}")	
endif(DEFINED KRATOS_INSTALL_PREFIX)


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("inc_dirs = ${inc_dirs}")

################################################################################
# install core files for the KratosMultiphysics python module
# note that this does not install KratosLoader.py, as it is different for installation and packaging runs (see next code block)
install(FILES "${CMAKE_SOURCE_DIR}/kratos/python_interface/__init__.py" DESTINATION KratosMultiphysics )
install(FILES "${CMAKE_SOURCE_DIR}/kratos/python_interface/kratos_utilities.py" DESTINATION KratosMultiphysics )
install(FILES "${CMAKE_SOURCE_DIR}/kratos/python_interface/kratos_globals.py" DESTINATION KratosMultiphysics )
install(FILES "${CMAKE_SOURCE_DIR}/kratos/python_interface/application_importer.py" DESTINATION KratosMultiphysics )

message("boost python lib used = " ${Boost_PYTHON_LIBRARY_RELEASE})
#install(FILES ${Boost_PYTHON_LIBRARY_RELEASE} DESTINATION libs)
#install(FILES ${EXTRA_INSTALL_LIBS} DESTINATION libs)

#install blas and lapack
if(${BLAS_INCLUDE_NEEDED} MATCHES ON )
	message("installed blas = " ${BLAS_LIBRARIES})
#	install(FILES ${BLAS_LIBRARIES} DESTINATION libs)
	message("installed lapack = " ${LAPACK_LIBRARIES})
#	install(FILES ${LAPACK_LIBRARIES} DESTINATION libs)
endif(${BLAS_INCLUDE_NEEDED} MATCHES ON )

################################################################################
if(${INSTALL_PYTHON_FILES} MATCHES ON)
    #take care: do not indent the lines below
file(WRITE ${CMAKE_SOURCE_DIR}/packaging_aux/KratosLoader.py "import sys
# # kratos_libs=\"${CMAKE_INSTALL_PREFIX}/libs\"
# # kratos_applications=\"${CMAKE_INSTALL_PREFIX}/applications\"
# # kratos_scripts=\"${CMAKE_INSTALL_PREFIX}/kratos/python_scripts\"
# # sys.path.append(kratos_libs)
# # sys.path.append(kratos_applications) 
# # sys.path.append(kratos_scripts) 
import os.path
kratos_libs=os.path.join(os.path.dirname(__file__),'../libs')
kratos_applications=os.path.join(os.path.dirname(__file__),'../applications')
kratos_scripts=os.path.join(os.path.dirname(__file__),'../kratos/python_scripts')
sys.path.append(kratos_libs)
sys.path.append(kratos_applications) 
sys.path.append(kratos_scripts) 
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
sys.path.append(kratos_libs)
sys.path.append(kratos_applications) 
sys.path.append(kratos_scripts) ")
endif(${INSTALL_PYTHON_FILES} MATCHES ON)
