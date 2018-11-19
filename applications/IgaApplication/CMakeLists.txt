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
    ${ANURBS_ROOT}/include
)

# generate variables with the sources
set(KRATOS_IGA_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/iga_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/iga_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/iga_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/iga_truss_element.cpp
    # utilities
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/iga_flags.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/brep_json_io.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/nurbs_brep_modeler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/brep_topology/brep_edge.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/brep_topology/brep_face.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/brep_topology/brep_model.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/brep_topology/brep_vertex.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/brep_topology/brep_boundary_loop.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/brep_topology/brep_trimming_curve.cpp
    # elements
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/base_discrete_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/surface_base_discrete_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_kl_discrete_element.cpp
)


# define library Kratos which defines the basic python interface
pybind11_add_module(KratosIgaApplication MODULE THIN_LTO ${KRATOS_IGA_APPLICATION_SOURCES})
target_link_libraries(KratosIgaApplication PRIVATE KratosCore)
set_target_properties(KratosIgaApplication PROPERTIES PREFIX "")

if(USE_COTIRE MATCHES ON)
    cotire(KratosIgaApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosIgaApplication DESTINATION libs)

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosIgaApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosIgaApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/IgaApplication.py" DESTINATION KratosMultiphysics)
