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

add_library(KratosLinearSolversCore SHARED)

set( CMAKE_INCLUDE_CURRENT_DIR ON )

message( "**** configuring KratosLinearSolversApplication ****" )

################### PYBIND11
include(pybind11Tools)

add_definitions( -DEIGEN_MPL2_ONLY)

include_directories(${KRATOS_SOURCE_DIR}/kratos)

set(KRATOS_LINEAR_SOLVERS_APPLICATION_LIB_INC "external_libraries/eigen3" "external_libraries/spectra1/include")

if( USE_EIGEN_MKL MATCHES ON )
    if( DEFINED ENV{MKLROOT} )
        message( "-- MKLROOT = $ENV{MKLROOT}" )

        set(KRATOS_LINEAR_SOLVERS_APPLICATION_LIB_INC ${KRATOS_LINEAR_SOLVERS_APPLICATION_LIB_INC} $ENV{MKLROOT}/include)

        link_directories("$ENV{MKLROOT}/lib")

        if( NOT MSVC )
            set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -L$ENV{MKLROOT}/lib/intel64 -Wl,--no-as-needed -lpthread -lm -ldl" )
        endif()
    elseif( DEFINED ENV{CONDA_PREFIX} )
        message("-- Found Conda environment: $ENV{CONDA_PREFIX}")

        if( WIN32 )
            set(KRATOS_LINEAR_SOLVERS_APPLICATION_LIB_INC ${KRATOS_LINEAR_SOLVERS_APPLICATION_LIB_INC} $ENV{CONDA_PREFIX}/Library/include)
            link_directories("$ENV{CONDA_PREFIX}/Library/lib")
        else()
            set(KRATOS_LINEAR_SOLVERS_APPLICATION_LIB_INC ${KRATOS_LINEAR_SOLVERS_APPLICATION_LIB_INC} $ENV{CONDA_PREFIX}/include)
            link_directories("$ENV{CONDA_PREFIX}/lib")
        endif()
    else()
        message( FATAL_ERROR "MKLROOT not defined" )
    endif()

    add_definitions( -DUSE_EIGEN_MKL -DEIGEN_USE_MKL_ALL )
    target_sources(KratosLinearSolversCore PRIVATE
                   "${CMAKE_CURRENT_SOURCE_DIR}/custom_solvers/mkl_solver_base.cpp"
                   "${CMAKE_CURRENT_SOURCE_DIR}/custom_solvers/mkl_ilu.cpp")
endif()

if( USE_EIGEN_FEAST MATCHES ON )
    if( WIN32 )
        message( FATAL_ERROR "FEAST4 is not yet supported under Windows")
    endif()
    if( USE_EIGEN_MKL MATCHES OFF OR NOT DEFINED USE_EIGEN_MKL )
        message( FATAL_ERROR "FEAST4 can only be used with MKL" )
    endif()
    include(CMakeAddFortranSubdirectory)
    cmake_add_fortran_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external_libraries/FEAST NO_EXTERNAL_INSTALL )
    set(FEAST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/external_libraries/FEAST/4.0 )
    set(KRATOS_LINEAR_SOLVERS_APPLICATION_LIB_INC ${KRATOS_LINEAR_SOLVERS_APPLICATION_LIB_INC} ${CMAKE_CURRENT_SOURCE_DIR}/external_libraries/FEAST/4.0/include)
    add_definitions(-DUSE_EIGEN_FEAST)
    # MKL implements some of FEAST's function. This definition excludes the definitions of MKL's FEAST
    # functions in mkl/include/mkl_solvers_ee.h. Otherwise, some functions are defined multiple times
    # defining both as MKL uses either, depending on the version
    add_definitions(-D_MKL_SOLVERS_EE_H)
    add_definitions(-D__MKL_SOLVERS_EE_H)
endif()

## generate variables with the sources
target_sources(KratosLinearSolversCore PRIVATE
               ${CMAKE_CURRENT_SOURCE_DIR}/linear_solvers_application.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/custom_factories/dense_linear_solver_factory.cpp)

# Sources for the Python module
file(
    GLOB_RECURSE
    KRATOS_LINEARSOLVERS_APPLICATION_PYTHON_INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

###############################################################
## LinearSolversApplication core library (C++ parts)
target_link_libraries( KratosLinearSolversCore PUBLIC KratosCore )
set_target_properties( KratosLinearSolversCore PROPERTIES COMPILE_DEFINITIONS "LINEARSOLVERS_APPLICATION=EXPORT,API")

target_include_directories(KratosLinearSolversCore SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${KRATOS_LINEAR_SOLVERS_APPLICATION_LIB_INC})

## LinearSolversApplication python module
pybind11_add_module( KratosLinearSolversApplication MODULE THIN_LTO ${KRATOS_LINEARSOLVERS_APPLICATION_PYTHON_INTERFACE} )
target_link_libraries( KratosLinearSolversApplication PRIVATE KratosLinearSolversCore )
set_target_properties( KratosLinearSolversApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosLinearSolversCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosLinearSolversApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

if( USE_EIGEN_MKL MATCHES ON )
    if( MSVC )
        find_library(MKL_RT_LIB mkl_rt)
        if(NOT MKL_RT_LIB)
            message( FATAL_ERROR "mkl_rt.lib not found")
        else(NOT MKL_RT_LIB)
            message( "mkl_rt.lib found at: ${MKL_RT_LIB}")
            target_link_libraries( KratosLinearSolversCore PUBLIC ${MKL_RT_LIB} )
        endif(NOT MKL_RT_LIB)
    else( MSVC )
        target_link_libraries( KratosLinearSolversCore PUBLIC mkl_rt )
        if( USE_EIGEN_FEAST MATCHES ON )
            target_link_libraries(KratosLinearSolversCore PUBLIC feast4 gfortran m)
        endif()
    endif( MSVC )
endif()

# changing the .dll suffix to .pyd
if( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
    set_target_properties( KratosLinearSolversCore PROPERTIES SUFFIX .pyd )
    set_target_properties( KratosLinearSolversApplication PROPERTIES SUFFIX .pyd )
endif( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/LinearSolversApplication.py KratosMultiphysics/LinearSolversApplication/__init__.py )

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
install( TARGETS KratosLinearSolversCore DESTINATION libs )
install( TARGETS KratosLinearSolversApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosLinearSolversCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosLinearSolversApplication" PARENT_SCOPE)