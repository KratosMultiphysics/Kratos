set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosOptimizationApplication ****")

# find_package(OpenSubdiv REQUIRED)

################### PYBIND11
include(pybind11Tools)

# include_directories( ${OPENSUBDIV_INCLUDE_DIR} )
include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## Optimization core sources
file( GLOB_RECURSE KRATOS_OPTIMIZATION_APPLICATION_CORE
	${CMAKE_CURRENT_SOURCE_DIR}/optimization_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/optimization_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_controls/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_responsess/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
    # could add here option to include subd if needed to avoid failing build
)

## Optimization python interface sources
file( GLOB_RECURSE KRATOS_OPTIMIZATION_APPLICATION_PYTHON_INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

# Optimization testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_OPTIMIZATION_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
    kratos_add_gtests(TARGET KratosOptimizationCore SOURCES "${KRATOS_OPTIMIZATION_APPLICATION_TESTING_SOURCES}" USE_CUSTOM_MAIN)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_library(KratosOptimizationCore SHARED ${KRATOS_OPTIMIZATION_APPLICATION_CORE})
target_link_libraries(KratosOptimizationCore PUBLIC KratosCore PRIVATE ${OPENSUBDIV_LIBRARIES})
set_target_properties(KratosOptimizationCore PROPERTIES COMPILE_DEFINITIONS "OPTIMIZATION_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosOptimizationApplication MODULE THIN_LTO ${KRATOS_OPTIMIZATION_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosOptimizationApplication PRIVATE KratosOptimizationCore)
set_target_properties(KratosOptimizationApplication PROPERTIES PREFIX "")

# Compiling the fcpw library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/fcpw)
# add_subdirectory(custom_external_libraries/fcpw)
# target_link_libraries(KratosOptimizationCore fcpw)
target_include_directories(KratosOptimizationCore PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/fcpw/include)
target_include_directories(KratosOptimizationCore PRIVATE ${FCPW_EIGEN_INCLUDES})
target_include_directories(KratosOptimizationCore PRIVATE ${FCPW_ENOKI_INCLUDES})

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosOptimizationCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosOptimizationApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosOptimizationApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosOptimizationApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/OptimizationApplication.py KratosMultiphysics/OptimizationApplication/__init__.py )

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
install(TARGETS KratosOptimizationCore DESTINATION libs )
install(TARGETS KratosOptimizationApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosOptimizationCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosOptimizationApplication" PARENT_SCOPE)