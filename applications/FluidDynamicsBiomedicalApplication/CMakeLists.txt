set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFluidDynamicsBiomedicalApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/kratos/applications/FluidDynamicsApplication)

## define application core sources (pure C++ parts)
set( KRATOS_FLUID_DYNAMICS_BIOMEDICAL_APPLICATION_CORE_SOURCES

  # application class and variables
  ${CMAKE_CURRENT_SOURCE_DIR}/fluid_dynamics_biomedical_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fluid_dynamics_biomedical_application_variables.cpp

  # utilities
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/wss_statistics_utilities.cpp

  # elements (from FluidDynamicsApplication)
  # conditions (from FluidDynamicsApplication)
  # fluid constitutive laws (from FluidDynamicsApplication)
  # processes (from FluidDynamicsApplication)
)

## define python interface sources
set( KRATOS_FLUID_DYNAMICS_BIOMEDICAL_PYTHON_INTERFACE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/fluid_dynamics_biomedical_python_application.cpp
)

## define fluid dynamics test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
 	file(GLOB_RECURSE KRATOS_FLUID_DYNAMICS_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################################
## FluidDynamicsBiomedicalApplication core library (C++ parts)
add_library( KratosFluidDynamicsBiomedicalCore SHARED ${KRATOS_FLUID_DYNAMICS_BIOMEDICAL_APPLICATION_CORE_SOURCES} ${KRATOS_FLUID_DYNAMICS_TESTING_SOURCES} )
target_link_libraries( KratosFluidDynamicsBiomedicalCore PUBLIC KratosCore KratosFluidDynamicsCore)
set_target_properties( KratosFluidDynamicsBiomedicalCore PROPERTIES COMPILE_DEFINITIONS "FLUID_DYNAMICS_BIOMEDICAL_APPLICATION=EXPORT,API")

## FluidDynamicsBiomedicalApplication python module
pybind11_add_module( KratosFluidDynamicsBiomedicalApplication MODULE THIN_LTO ${KRATOS_FLUID_DYNAMICS_BIOMEDICAL_PYTHON_INTERFACE_SOURCES} )
target_link_libraries( KratosFluidDynamicsBiomedicalApplication PRIVATE KratosFluidDynamicsBiomedicalCore)
set_target_properties( KratosFluidDynamicsBiomedicalApplication PROPERTIES PREFIX "")

###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosFluidDynamicsBiomedicalApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosFluidDynamicsBiomedicalApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

###############################################################################
if(USE_COTIRE MATCHES ON)
    cotire(KratosFluidDynamicsBiomedicalCore)
    cotire(KratosFluidDynamicsBiomedicalApplication)
endif(USE_COTIRE MATCHES ON)

###############################################################################
## installing the resulting libraries

install(TARGETS KratosFluidDynamicsBiomedicalCore DESTINATION libs )
install(TARGETS KratosFluidDynamicsBiomedicalApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/FluidDynamicsBiomedicalApplication.py" DESTINATION "KratosMultiphysics/FluidDynamicsBiomedicalApplication" RENAME "__init__.py")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

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

###############################################################################
## If trilinos is avaliable, switch on the Fluid trilinos extension

if((${MPI_NEEDED} MATCHES ON) AND (${TRILINOS_APPLICATION} MATCHES ON))
  message("Enabling trilinos extension for KratosFluidDynamicsBiomedicalApplication")
  add_subdirectory(trilinos_extension)
endif((${MPI_NEEDED} MATCHES ON) AND (${TRILINOS_APPLICATION} MATCHES ON) )
