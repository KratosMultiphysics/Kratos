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

if( NOT DEFINED EIGEN_ROOT )
	message( FATAL_ERROR "EIGEN_ROOT not defined" )
endif()

if( USE_EIGEN_MKL MATCHES ON AND NOT DEFINED ENV{MKLROOT} )
	message( FATAL_ERROR "MKLROOT not defined" )
endif()

add_definitions( -DEIGEN_DEFAULT_TO_ROW_MAJOR -DEIGEN_MPL2_ONLY )

include_directories(
	${CMAKE_SOURCE_DIR}/kratos
        SYSTEM ${EIGEN_ROOT} # Including as System to suppress compile-warnings from it
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
set( KRATOS_EIGENSOLVERS_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/eigen_solvers_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/eigen_solvers_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_solvers_to_python.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module( KratosEigenSolversApplication MODULE THIN_LTO ${KRATOS_EIGENSOLVERS_APPLICATION_SOURCES} )
target_link_libraries( KratosEigenSolversApplication PRIVATE KratosCore )

if( USE_EIGEN_MKL MATCHES ON )
    if( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
        find_library(MKL_RT_LIB mkl_rt)
        message( "mkl_rt.lib found at: ${MKL_RT_LIB}")
        target_link_libraries( KratosEigenSolversApplication PRIVATE ${MKL_RT_LIB} )
    else( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
        target_link_libraries( KratosEigenSolversApplication PRIVATE mkl_rt )
    endif( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
endif()

set_target_properties( KratosEigenSolversApplication PROPERTIES PREFIX "" )
install( TARGETS KratosEigenSolversApplication DESTINATION libs )

# changing the .dll suffix to .pyd
if( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
	set_target_properties( KratosEigenSolversApplication PROPERTIES SUFFIX .pyd )
endif( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )

if( ${INSTALL_PYTHON_FILES} MATCHES ON )
	get_filename_component ( CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME )
	install( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE )
endif( ${INSTALL_PYTHON_FILES} MATCHES ON )

if(USE_COTIRE MATCHES ON)
    cotire(KratosEigenSolversApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
install( FILES "${CMAKE_CURRENT_SOURCE_DIR}/EigenSolversApplication.py" DESTINATION KratosMultiphysics )
