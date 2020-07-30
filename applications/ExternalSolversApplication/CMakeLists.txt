set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** Configuring KratosExternalSolversApplication ****")

include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/SuperLU_4.3 )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries )

################### PYBIND11
include(pybind11Tools)

########### Add Libblas and lapack to compilation
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

## Add subdirectories for compilation of dependent libraries
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/SuperLU_4.3)

if(${INCLUDE_PASTIX} MATCHES ON)
    add_definitions(-DINCLUDE_PASTIX)
    include_directories( ${PASTIX_INSTALL_DIR} )
    FIND_LIBRARY(PASTIX_LINK_LIB pastix  ${PASTIX_INSTALL_DIR})
    FIND_LIBRARY(SCOTCH_ERR_LIB scotcherr ${SCOTCH_INSTALL_DIR})
    FIND_LIBRARY(SCOTCH_LIB scotch  ${SCOTCH_INSTALL_DIR})
    message("Pastix library found =" ${PASTIX_LINK_LIB})
    add_definitions( -DFORCE_NOMPI )
    add_library(kratos_pastix STATIC custom_external_libraries/solvePASTIX.c  )
    target_link_libraries(kratos_pastix  ${PASTIX_LINK_LIB} ${SCOTCH_LIB} ${SCOTCH_ERR_LIB} rt  ${BLAS_LIBRARIES} )

    install(FILES ${SCOTCH_LIB} ${SCOTCH_ERR_LIB} ${BLAS_LIBRARIES} DESTINATION libs )
else(${INCLUDE_PASTIX} MATCHES ON)
    set(PASTIX_LINK_LIB "")
endif(${INCLUDE_PASTIX} MATCHES ON)

if(${INCLUDE_FEAST} MATCHES ON)
  include(CMakeAddFortranSubdirectory)
  cmake_add_fortran_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/FEAST NO_EXTERNAL_INSTALL )
  if (NOT DEFINED FEAST_ROOT)
    set(FEAST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/FEAST/3.0 )
  endif(NOT DEFINED FEAST_ROOT)
  include_directories( ${FEAST_ROOT}/include )
#   link_directories( ${FEAST_ROOT}/lib/x64 )
  add_definitions(-DINCLUDE_FEAST)
endif(${INCLUDE_FEAST} MATCHES ON)

## Generate variables with the sources
set( KRATOS_EXTERNAL_SOLVERS_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/externalsolvers_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_linear_solvers_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/externalsolvers_python_application.cpp
)

###############################################################
## Define library Kratos which defines the basic python interface
pybind11_add_module(KratosExternalSolversApplication MODULE THIN_LTO ${KRATOS_EXTERNAL_SOLVERS_APPLICATION_SOURCES})

if(${INCLUDE_PASTIX} MATCHES ON)
        target_link_libraries(KratosExternalSolversApplication PRIVATE KratosCore superlu  kratos_pastix rt)
else(${INCLUDE_PASTIX} MATCHES ON)
        target_link_libraries(KratosExternalSolversApplication PRIVATE KratosCore superlu  )
endif(${INCLUDE_PASTIX} MATCHES ON)

if(${INCLUDE_FEAST} MATCHES ON)
  set(FEAST_LINK_LIBS feast ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} gfortran m)
  target_link_libraries(KratosExternalSolversApplication PRIVATE ${FEAST_LINK_LIBS})
endif(${INCLUDE_FEAST} MATCHES ON)

set_target_properties(KratosExternalSolversApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
# if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#     set_target_properties(KratosExternalSolversApplication PROPERTIES SUFFIX .pyd)
# endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#
# # changing the .dylib suffix to .so
# if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#     set_target_properties(KratosExternalSolversApplication PROPERTIES SUFFIX .so)
# endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosExternalSolversApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ExternalSolversApplication.py KratosMultiphysics/ExternalSolversApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
        PATTERN "*.git" EXCLUDE
        PATTERN "*.c" EXCLUDE
        PATTERN "*.h" EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.hpp" EXCLUDE
  )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosExternalSolversApplication DESTINATION libs )

# Install blas and lapack
install(FILES ${BLAS_LIBRARIES} DESTINATION libs)
install(FILES ${LAPACK_LIBRARIES} DESTINATION libs)

