#
#   KRATOS _______
#         / ____(_)___ ____  ____
#        / __/ / / __ `/ _ \/ __ \
#       / /___/ / /_/ /  __/ / / /
#      /_____/_/\__, /\___/_/ /_/ SolversApplication
#              /____/
#
#   Author: Thomas Oberbichler
#

set( CMAKE_INCLUDE_CURRENT_DIR ON )

message( "**** configuring KratosEigenSolversApplication ****" )

################### PYBIND11
include(pybind11Tools)

if( USE_EIGEN_MKL MATCHES ON AND NOT DEFINED ENV{MKLROOT} )
    message( FATAL_ERROR "MKLROOT not defined" )
endif()

add_definitions( -DEIGEN_DEFAULT_TO_ROW_MAJOR -DEIGEN_MPL2_ONLY )

include_directories(
    ${KRATOS_SOURCE_DIR}/kratos
        SYSTEM external_libraries/eigen3 # Including as System to suppress compile-warnings from it
)

if( USE_EIGEN_MKL MATCHES ON )
    message( "MKLROOT = $ENV{MKLROOT}" )

    if( NOT MSVC )
        set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -L$ENV{MKLROOT}/lib/intel64 -Wl,--no-as-needed -lpthread -lm -ldl" )
    endif()

    include_directories( $ENV{MKLROOT}/include )
    add_definitions( -DUSE_EIGEN_MKL -DEIGEN_USE_MKL_ALL )
endif()

## generate variables with the sources
set( KRATOS_EIGENSOLVERS_APPLICATION_CORE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/eigen_solvers_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_factories/dense_linear_solver_factory.cpp
)

# Sources for the Python module
file(
    GLOB_RECURSE
    KRATOS_EIGEN_PYTHON_INTERFACE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

###############################################################
## EigenSolversApplication core library (C++ parts)
add_library( KratosEigenSolversCore SHARED ${KRATOS_EIGENSOLVERS_APPLICATION_CORE_SOURCES} )
target_link_libraries( KratosEigenSolversCore PUBLIC KratosCore )
set_target_properties( KratosEigenSolversCore PROPERTIES COMPILE_DEFINITIONS "EIGENSOLVERS_APPLICATION=EXPORT,API")

## EigenSolversApplication python module
pybind11_add_module( KratosEigenSolversApplication MODULE THIN_LTO ${KRATOS_EIGEN_PYTHON_INTERFACE_SOURCES} )
target_link_libraries( KratosEigenSolversApplication PRIVATE KratosEigenSolversCore )
set_target_properties( KratosEigenSolversApplication PROPERTIES PREFIX "")

if( USE_EIGEN_MKL MATCHES ON )
    if( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
        find_library(MKL_RT_LIB mkl_rt)
        message( "mkl_rt.lib found at: ${MKL_RT_LIB}")
        target_link_libraries( KratosEigenSolversCore PUBLIC ${MKL_RT_LIB} )
    elseif( ${CMAKE_CXX_COMPILER_ID} MATCHES Clang )
        message( FATAL_ERROR "Clang does not yet support MKL" )
    else( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
        target_link_libraries( KratosEigenSolversCore PUBLIC mkl_rt )
    endif( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
endif()

# changing the .dll suffix to .pyd
if( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
    set_target_properties( KratosEigenSolversCore PROPERTIES SUFFIX .pyd )
    set_target_properties( KratosEigenSolversApplication PROPERTIES SUFFIX .pyd )
endif( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosEigenSolversCore)
    cotire(KratosEigenSolversApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/EigenSolversApplication.py KratosMultiphysics/EigenSolversApplication/__init__.py )

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
install( TARGETS KratosEigenSolversCore DESTINATION libs )
install( TARGETS KratosEigenSolversApplication DESTINATION libs )
