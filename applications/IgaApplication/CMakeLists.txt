set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosIgaApplication ****")

################### PYBIND11
include(pybind11Tools)

if(NOT DEFINED ANURBS_ROOT)
    message(FATAL_ERROR
        "ANURBS_ROOT not defined!\n"
        "  - Download ANurbs from https://github.com/oberbichler/ANurbs\n"
        "  - Add '-DANURBS_ROOT=<path to ANurbs>' to your build script")
else()
    message(STATUS "ANURBS_ROOT = ${ANURBS_ROOT}")
endif()

include_directories(
    ${CMAKE_SOURCE_DIR}/kratos
    ${CMAKE_SOURCE_DIR}/external_libraries/triangle
    ${ANURBS_ROOT}/include
    ${ANURBS_ROOT}/external_libraries
)


# generate variables with the sources
set(KRATOS_IGA_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/iga_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/iga_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/iga_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp

    # utilities
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/iga_flags.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/brep_json_io.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/nurbs_brep_modeler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/embedded_iga_modeler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/embedded_iga/embedded_iga_tessellation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/embedded_iga/embedded_iga_triangulation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/embedded_iga/embedded_iga_error_estimation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/brep_topology/brep_edge.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/brep_topology/brep_face.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/brep_topology/brep_model.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/brep_topology/brep_vertex.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/brep_topology/brep_boundary_loop.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/brep_topology/brep_trimming_curve.cpp

    ${ANURBS_ROOT}/external_libraries/ClipperLib/clipper.cpp
)

# Define python interface sources
file(GLOB_RECURSE KRATOS_IGA_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

# Generate variables with the testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
  file(GLOB_RECURSE KRATOS_IGA_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosIgaApplication MODULE THIN_LTO ${KRATOS_IGA_APPLICATION_SOURCES} ${KRATOS_IGA_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosIgaApplication PRIVATE KratosCore triangle)
set_target_properties(KratosIgaApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosIgaApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosIgaApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosIgaCore)
    cotire(KratosIgaApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/IgaApplication.py KratosMultiphysics/IgaApplication/__init__.py )

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
install(TARGETS KratosIgaCore DESTINATION libs )
install(TARGETS KratosIgaApplication DESTINATION libs )
