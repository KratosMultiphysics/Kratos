set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMeshMovingApplication ****")

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## MeshMoving Core sources
file(GLOB_RECURSE KRATOS_MESH_MOVING_APPLICATION_CORE
	${CMAKE_CURRENT_SOURCE_DIR}/mesh_moving_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
)

## MeshMoving python interface sources
file(GLOB_RECURSE KRATOS_MESH_MOVING_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

## define mesh moving test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
 	file(GLOB_RECURSE KRATOS_MESH_MOVING_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
  kratos_add_gtests(TARGET KratosMeshMovingCore SOURCES "${KRATOS_MESH_MOVING_TESTING_SOURCES}" USE_CUSTOM_MAIN)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_library(KratosMeshMovingCore SHARED ${KRATOS_MESH_MOVING_APPLICATION_CORE})
target_link_libraries(KratosMeshMovingCore PUBLIC KratosCore)
set_target_properties(KratosMeshMovingCore PROPERTIES COMPILE_DEFINITIONS "MESH_MOVING_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosMeshMovingApplication MODULE THIN_LTO ${KRATOS_MESH_MOVING_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosMeshMovingApplication PRIVATE KratosMeshMovingCore)
set_target_properties(KratosMeshMovingApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosMeshMovingCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosMeshMovingApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosMeshMovingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosMeshMovingApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/MeshMovingApplication.py KratosMultiphysics/MeshMovingApplication/__init__.py )

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
install(TARGETS KratosMeshMovingCore DESTINATION libs )
install(TARGETS KratosMeshMovingApplication DESTINATION libs )

## If trilinos is avaliable, switch on the Fluid trilinos extension
if((${USE_MPI} MATCHES ON) AND (${TRILINOS_FOUND}))
  message("Enabling trilinos extension for KratosMeshMovingApplication")
  add_subdirectory(trilinos_extension)
endif((${USE_MPI} MATCHES ON) AND (${TRILINOS_FOUND}))

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosMeshMovingCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosMeshMovingApplication" PARENT_SCOPE)