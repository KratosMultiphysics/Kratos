set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosGeoMechanicsApplication ****")

# PYBIND11
include(pybind11Tools)

add_compile_definitions(_USE_MATH_DEFINES)

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/LinearSolversApplication)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/LinearSolversApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/GeoMechanicsApplication )

## GeoMechanics Core sources
file(GLOB_RECURSE KRATOS_GEO_MECHANICS_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/geo_mechanics_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp_geomechanics_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/geo_mechanics_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_retention/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_operations/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_workflows/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/external_bindings/*.cpp
)


# Shared setup functionality between unit test and microbenchmarking
file(GLOB_RECURSE TestSetupSources
        ${CMAKE_CURRENT_SOURCE_DIR}/test_setup_utilities/*.cpp
)

add_library(KratosGeoMechanicsTestSetupUtilities SHARED ${TestSetupSources})
target_link_libraries(KratosGeoMechanicsTestSetupUtilities KratosGeoMechanicsCore)
set_target_properties(KratosGeoMechanicsTestSetupUtilities PROPERTIES COMPILE_DEFINITIONS "KRATOS_GEO_TEST_UTILS=IMPORT,API" CXX_STANDARD 20)

## GeoMechanics testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_GEO_MECHANICS_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
    kratos_add_gtests(TARGET KratosGeoMechanicsCore SOURCES "${KRATOS_GEO_MECHANICS_APPLICATION_TESTING_SOURCES}")
    target_link_libraries(KratosGeoMechanicsCoreTest KratosGeoMechanicsTestSetupUtilities)
    set_target_properties("KratosGeoMechanicsCoreTest" PROPERTIES CXX_STANDARD 20)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## GeoMechanics python interface sources
file(GLOB_RECURSE KRATOS_GEO_MECHANICS_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosGeoMechanicsCore SHARED ${KRATOS_GEO_MECHANICS_APPLICATION_CORE})
target_link_libraries(KratosGeoMechanicsCore PUBLIC KratosCore KratosStructuralMechanicsCore KratosLinearSolversCore ${CMAKE_DL_LIBS})
set_target_properties(KratosGeoMechanicsCore PROPERTIES COMPILE_DEFINITIONS "GEO_MECHANICS_APPLICATION=EXPORT,API" CXX_STANDARD 20)


###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosGeoMechanicsApplication MODULE THIN_LTO ${KRATOS_GEO_MECHANICS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosGeoMechanicsApplication PRIVATE KratosGeoMechanicsCore)
set_target_properties(KratosGeoMechanicsApplication PROPERTIES PREFIX "" CXX_STANDARD 20)

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosGeoMechanicsCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosGeoMechanicsApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# Add version info to DLL's
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.rc.in ${CMAKE_CURRENT_BINARY_DIR}/version.rc @ONLY)
    target_sources(KratosGeoMechanicsCore PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/version.rc)
    target_sources(KratosGeoMechanicsApplication PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/version.rc)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosGeoMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosGeoMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/GeoMechanicsApplication.py KratosMultiphysics/GeoMechanicsApplication/__init__.py )
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_registry_lists.py KratosMultiphysics/GeoMechanicsApplication/python_registry_lists.py )

if(CMAKE_UNITY_BUILD MATCHES ON)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_small_strain_element.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_small_strain_FIC_element.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
endif(CMAKE_UNITY_BUILD MATCHES ON)

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )


# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    kratos_python_install_tests(
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tests"
        DESTINATION "applications/${CURRENT_DIR_NAME}" 
        INSTALL_PYTHON_USING_LINKS "${INSTALL_PYTHON_USING_LINKS}"
        EXCLUSION_PATTERNS "\.git$" "\.c$" "\.h$" "\.cpp$" "\.hpp$"
    )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosGeoMechanicsCore DESTINATION libs )
install(TARGETS KratosGeoMechanicsApplication DESTINATION libs )
install(TARGETS KratosGeoMechanicsTestSetupUtilities DESTINATION libs)

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosGeoMechanicsCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosGeoMechanicsApplication" PARENT_SCOPE)

if(${KRATOS_BUILD_BENCHMARK} MATCHES ON)
    file(GLOB_RECURSE KRATOS_BENCHMARK_SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/*.cpp
    )

    foreach(file ${KRATOS_BENCHMARK_SOURCES})
        get_filename_component(filename ${file} NAME_WE)
        add_executable(${filename} ${file})
        target_link_libraries(${filename} PUBLIC KratosGeoMechanicsCore KratosGeoMechanicsTestSetupUtilities benchmark::benchmark)
        set_target_properties(${filename} PROPERTIES COMPILE_DEFINITIONS "KRATOS_GEO_BENCHMARK=IMPORT,API")
        set_target_properties(${filename} PROPERTIES CXX_STANDARD 20)

        install(TARGETS ${filename} DESTINATION benchmark)
    endforeach(file ${KRATOS_BENCHMARK_SOURCES})
endif(${KRATOS_BUILD_BENCHMARK} MATCHES ON)