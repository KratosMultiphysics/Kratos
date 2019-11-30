set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosRANSApplication ****")

# check for fluid dynamics application
if(NOT (${FLUID_DYNAMICS_APPLICATION} MATCHES ON))
  message(FATAL_ERROR "RANSApplication requires FluidDynamicsApplication, Please re-compile with -DFLUID_DYNAMICS_APPLICATION=ON flag.")
endif(NOT (${FLUID_DYNAMICS_APPLICATION} MATCHES ON))

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication )

if((${MPI_NEEDED} MATCHES ON) AND (${TRILINOS_APPLICATION} MATCHES ON))
  message("Enabling trilinos extension for KratosRANSApplication")
  add_subdirectory(trilinos_extension)
endif((${MPI_NEEDED} MATCHES ON) AND (${TRILINOS_APPLICATION} MATCHES ON) )

## RANSApplication Core sources
file(GLOB_RECURSE KRATOS_RANS_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/rans_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/rans_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

file(GLOB_RECURSE KRATOS_RANS_APPLICATION_PYTHON_INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

if(${KRATOS_BUILD_TESTING} MATCHES ON)
  file(GLOB_RECURSE
    KRATOS_RANS_APPLICATION_TESTING
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp/*.cpp
  )
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_library(KratosRANSCore SHARED ${KRATOS_RANS_APPLICATION_CORE} ${KRATOS_RANS_APPLICATION_TESTING})
target_link_libraries(KratosRANSCore PUBLIC KratosCore KratosFluidDynamicsCore)
set_target_properties(KratosRANSCore PROPERTIES COMPILE_DEFINITIONS "RANS_APPLICATION=EXPORT,API")

pybind11_add_module(KratosRANSApplication MODULE THIN_LTO ${KRATOS_RANS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosRANSApplication PRIVATE KratosRANSCore)
set_target_properties(KratosRANSApplication PROPERTIES PREFIX "")

if(USE_COTIRE MATCHES ON)
    cotire(KratosRANSCore)
    cotire(KratosRANSApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosRANSCore DESTINATION libs )
install(TARGETS KratosRANSApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosRANSApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosRANSApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/RANSApplication.py" DESTINATION "KratosMultiphysics/RANSApplication" RENAME "__init__.py")
