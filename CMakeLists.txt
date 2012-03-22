cmake_minimum_required (VERSION 2.8)
project (KratosBaseProject)

#set compiler flags
if(${CMAKE_COMPILER_IS_GNUCXX}) 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops -ffast-math -Wall -Wno-strict-aliasing -ansi")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -funroll-loops -ffast-math -Wall -Wno-strict-aliasing")
  message("additional default options were set for gcc")
  message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
  message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
endif(${CMAKE_COMPILER_IS_GNUCXX}) 
if(${CMAKE_C_COMPILER} MATCHES "icc.*$") 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -ansi -funroll-loops  -lpthread -wd654 -wd10010  ")
  set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops  -lpthread -wd654 -wd10010  ")
  message("additional default options were set for intel compiler")
  message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
  message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
endif(${CMAKE_C_COMPILER} MATCHES "icc.*$") 

if(${MSVC80})
  message("detected compiler as MSVC80")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS")
endif(${MSVC80}) 

if(${MSVC90})
  message("detected compiler as MSVC90")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS")
endif(${MSVC90}) 

if(${MSVC10})
  message("detected compiler as MSVC10")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc -D_SCL_SECURE_NO_WARNINGS")
endif(${MSVC10}) 

#tell the linker to give an error if undefined functions are found
set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")

#set kratos specific module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")

######################################################################################
######################################################################################
######################################################################################
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

    find_package(BLAS REQUIRED)
    include_directories(${BLAS_INCLUDE_DIRS})
    link_directories(${BLAS_LIBRARIES})
    set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}" CACHE STRING "Linker flags for shared libs" FORCE)
    message("linker flags with blas =" ${BLAS_LINKER_FLAGS})
    message("BLA_VENDOR = " ${BLA_VENDOR})

    find_package(LAPACK REQUIRED)
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

if(${MPI_NEEDED} MATCHES ON )
    find_package(MPI REQUIRED)
    add_definitions( -DKRATOS_USING_MPI )
    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

    #note that this will be included in the compilation of the kratos core, so that it is not needed to include it somewhere else
    include_directories(${MPI_INCLUDE_PATH})
    find_package(Boost COMPONENTS mpi serialization REQUIRED)
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
# include_directories( ${CMAKE_SOURCE_DIR}/external_libraries/gidpost )

#defines needed
add_definitions( -DKRATOS_PYTHON )
add_definitions( -DNDEBUG )
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    add_definitions( -fPIC )
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

#include subdirectories
add_subdirectory(external_libraries/gidpost)
add_subdirectory(kratos)
add_subdirectory(applications)

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
install(FILES "${CMAKE_SOURCE_DIR}/kratos/python_interface/kratos_globals.py" DESTINATION KratosMultiphysics )
install(FILES "${CMAKE_SOURCE_DIR}/kratos/python_interface/application_importer.py" DESTINATION KratosMultiphysics )


################################################################################
if(${INSTALL_PYTHON_FILES} MATCHES ON)
    #take care: do not indent the lines below
file(WRITE ${CMAKE_SOURCE_DIR}/packaging_aux/KratosLoader.py "import sys
kratos_libs=\"${CMAKE_INSTALL_PREFIX}/libs\"
kratos_applications=\"${CMAKE_INSTALL_PREFIX}/applications\"
sys.path.append(kratos_libs)
sys.path.append(kratos_applications) ")

    file(WRITE ${CMAKE_SOURCE_DIR}/packaging_aux/kratos.conf "${CMAKE_INSTALL_PREFIX}/libs")

    install(FILES ${CMAKE_SOURCE_DIR}/packaging_aux/KratosLoader.py DESTINATION KratosMultiphysics ) 
    install(FILES ${CMAKE_SOURCE_DIR}/packaging_aux/kratos.conf DESTINATION packaging )
    install(FILES ${CMAKE_SOURCE_DIR}/packaging_aux/script.py DESTINATION . )


    install(FILES ${Boost_PYTHON_LIBRARIES} DESTINATION libs)
    install(FILES ${EXTRA_INSTALL_LIBS} DESTINATION libs)

#    GET_FILENAME_COMPONENT(PYTHON_LIBRARIES_REALPATH ${PYTHON_LIBRARIES} REALPATH)
#    install(FILES ${PYTHON_LIBRARIES_REALPATH} DESTINATION libs)
#    install(FILES ${PYTHON_LIBRARIES} DESTINATION libs)



    SET(CPACK_GENERATOR "DEB")
    SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Kratos Team") #required
    #SET(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
    SET(CPACK_DEBIAN_PACKAGE_DEPENDS ${CPACK_DEBIAN_PACKAGE_DEPENDS} "libblas3gf | libblas.so.3gf | libatlas3gf-base, libc6 (>= 2.7), libgcc1 (>= 1:4.1.1), libgfortran3 (>= 4.6), libgomp1 (>= 4.4), liblapack3gf | liblapack.so.3gf | libatlas3gf-base, libopenmpi1.3, libparmetis3.1, libstdc++6 (>= 4.4.0), ${KRATOS_EXTRA_PACKAGE_DEPENDS} " )
#    SET(CPACK_DEBIAN_PACKAGE_DEPENDS ${CPACK_DEBIAN_PACKAGE_DEPENDS} " openmpi-bin, libboost-mpi-python-dev, ${KRATOS_EXTRA_PACKAGE_DEPENDS}" )
    Set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
    set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "${CMAKE_SOURCE_DIR}/packaging_aux/postinst;${CMAKE_SOURCE_DIR}/packaging_aux/prerm;")

    INCLUDE(CPack)
else(${INSTALL_PYTHON_FILES} MATCHES ON)
    #take care: do not indent the lines below
    file(WRITE ${CMAKE_INSTALL_PREFIX}/KratosMultiphysics/KratosLoader.py "import sys
kratos_libs=\"${CMAKE_INSTALL_PREFIX}/libs\"
kratos_applications=\"${CMAKE_SOURCE_DIR}/applications\"
sys.path.append(kratos_libs)
sys.path.append(kratos_applications) ")
endif(${INSTALL_PYTHON_FILES} MATCHES ON)
