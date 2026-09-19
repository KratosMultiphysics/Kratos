set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosRailwayApplication ****")

# PYBIND11
include(pybind11Tools)

add_compile_definitions(_USE_MATH_DEFINES)

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/LinearSolversApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/GeoMechanicsApplication)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/LinearSolversApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/GeoMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/RailwayApplication )


## Railway Core sources
file(GLOB_RECURSE KRATOS_RAILWAY_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/railway_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp_railway_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)


## Railway python interface sources
file(GLOB_RECURSE KRATOS_RAILWAY_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosRailwayCore SHARED ${KRATOS_RAILWAY_APPLICATION_CORE})
target_link_libraries(KratosRailwayCore PUBLIC KratosCore KratosStructuralMechanicsCore KratosLinearSolversCore KratosGeoMechanicsCore ${CMAKE_DL_LIBS})
set_target_properties(KratosRailwayCore PROPERTIES COMPILE_DEFINITIONS "RAILWAY_APPLICATION=EXPORT,API")


###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosRailwayApplication MODULE THIN_LTO ${KRATOS_RAILWAY_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosRailwayApplication PRIVATE KratosRailwayCore)
set_target_properties(KratosRailwayApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosRailwayCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosRailwayApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# Add version info to DLL's
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.rc.in ${CMAKE_CURRENT_BINARY_DIR}/version.rc @ONLY)
    target_sources(KratosRailwayCore PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/version.rc)
    target_sources(KratosRailwayApplication PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/version.rc)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosRailwayApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosRailwayApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/RailwayApplication.py KratosMultiphysics/RailwayApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )


# Install targets
install(TARGETS KratosRailwayCore DESTINATION libs )
install(TARGETS KratosRailwayApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosRailwayCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosRailwayApplication" PARENT_SCOPE)