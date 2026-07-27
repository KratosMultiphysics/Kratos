set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMeshioPlusPlusApplication ****")

################### PYBIND11
include(pybind11Tools)

# meshio++ (https://github.com/loumalouomega/meshioplusplus) is consumed as a normal,
# independently built external dependency - nothing is vendored into Kratos.
#
# Build and install it with the C++ API and the KRATOS mesh backend, e.g.:
#
#   cmake -S <meshioplusplus> -B build -DMESHIOPLUSPLUS_INSTALL_CPP=ON \
#         -DMESHIOPLUSPLUS_INSTALL_CPP_BACKENDS="KRATOS" \
#         -DMESHIOPLUSPLUS_MESH_BACKEND=KRATOS -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF \
#         -DBUILD_SHARED_LIBS=ON
#   cmake --build build && cmake --install build --prefix <prefix>
#
# then point this configure at it with -Dmeshioplusplus_DIR=<prefix>/lib/cmake/meshioplusplus
#
# EXACT is deliberate, not stylistic: the meshio++ C++ API makes no ABI promise (Mesh,
# ModelPart and GeometricalEntity are header-defined types whose layout changes with the
# headers), so the application must be rebuilt whenever meshio++ is. See doc/cpp_api.md.
# All three version components are required: under SameMajorVersion, EXACT is a full string
# comparison, so "9.2" would not match an installed "9.2.0".
find_package(meshioplusplus 9.2.0 EXACT CONFIG REQUIRED COMPONENTS CXX)

if(NOT TARGET meshioplusplus::core_kratos)
    message(FATAL_ERROR
        "MeshioPlusPlusApplication requires the meshio++ KRATOS mesh backend, but the install at "
        "'${meshioplusplus_DIR}' does not carry it (it shipped: ${MESHIOPLUSPLUS_MESH_BACKENDS}). "
        "Reconfigure meshio++ with -DMESHIOPLUSPLUS_INSTALL_CPP_BACKENDS=\"KRATOS\" (or a list "
        "including it).")
endif()

message(STATUS "meshio++ ${meshioplusplus_VERSION}: backends '${MESHIOPLUSPLUS_MESH_BACKENDS}', "
               "parallel '${MESHIOPLUSPLUS_PARALLEL_BACKEND}', HDF5 '${MESHIOPLUSPLUS_WITH_HDF5}', "
               "netCDF '${MESHIOPLUSPLUS_WITH_NETCDF}'")

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## MeshioPlusPlus Core sources
file(GLOB_RECURSE KRATOS_MESHIOPLUSPLUS_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/meshioplusplus_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_modelers/*.cpp
)

## MeshioPlusPlus testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_MESHIOPLUSPLUS_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
    kratos_add_gtests(TARGET KratosMeshioPlusPlusCore SOURCES "${KRATOS_MESHIOPLUSPLUS_APPLICATION_TESTING_SOURCES}" USE_CUSTOM_MAIN)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## MeshioPlusPlus python interface sources
file(GLOB_RECURSE KRATOS_MESHIOPLUSPLUS_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosMeshioPlusPlusCore SHARED ${KRATOS_MESHIOPLUSPLUS_APPLICATION_CORE})
target_link_libraries(KratosMeshioPlusPlusCore PUBLIC KratosCore meshioplusplus::core_kratos)
set_target_properties(KratosMeshioPlusPlusCore PROPERTIES COMPILE_DEFINITIONS "KRATOS_MESHIOPLUSPLUS_APPLICATION=EXPORT,API")

# meshio++ requires C++20; KratosCore may be built with an older standard, so require it here.
target_compile_features(KratosMeshioPlusPlusCore PUBLIC cxx_std_20)

# The library version is not available through the C++ API (mio_version() is C-API only),
# so it is passed through from the CMake package for the start-up log message.
target_compile_definitions(KratosMeshioPlusPlusCore PRIVATE KRATOS_MESHIOPLUSPLUS_VERSION="${meshioplusplus_VERSION}")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosMeshioPlusPlusApplication MODULE THIN_LTO ${KRATOS_MESHIOPLUSPLUS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosMeshioPlusPlusApplication PRIVATE KratosMeshioPlusPlusCore)
set_target_properties(KratosMeshioPlusPlusApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosMeshioPlusPlusApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosMeshioPlusPlusApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/MeshioPlusPlusApplication.py KratosMultiphysics/MeshioPlusPlusApplication/__init__.py )
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_registry_lists.py KratosMultiphysics/MeshioPlusPlusApplication/python_registry_lists.py )

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
install(TARGETS KratosMeshioPlusPlusCore DESTINATION libs)
install(TARGETS KratosMeshioPlusPlusApplication DESTINATION libs)

# Where the meshio++ shared library actually lives, needed both for the RPATH below and for
# the optional TPL copy further down.
get_target_property(MESHIOPLUSPLUS_CORE_LOCATION meshioplusplus::core_kratos IMPORTED_LOCATION_RELEASE)
if(NOT MESHIOPLUSPLUS_CORE_LOCATION)
    get_target_property(MESHIOPLUSPLUS_CORE_LOCATION meshioplusplus::core_kratos IMPORTED_LOCATION)
endif()
if(MESHIOPLUSPLUS_CORE_LOCATION)
    get_filename_component(MESHIOPLUSPLUS_LIB_DIR "${MESHIOPLUSPLUS_CORE_LOCATION}" DIRECTORY)
endif()

# Kratos installs with an RPATH of just "$ORIGIN", which is enough for its own libraries but
# not for meshio++: it is an external dependency that normally lives in its own prefix
# (~/.local/meshioplusplus, /opt/..., a Conan cache), none of which the dynamic linker
# searches by default. Without its directory recorded here, importing the Python module
# fails with "libmeshioplusplus_core_kratos.so.0: cannot open shared object file" unless the
# user sets LD_LIBRARY_PATH by hand. $ORIGIN is kept first so an INSTALL_TPL_LIBRARIES copy
# in libs/ still wins.
if(MESHIOPLUSPLUS_LIB_DIR AND NOT WIN32)
    set_target_properties(KratosMeshioPlusPlusCore KratosMeshioPlusPlusApplication PROPERTIES
        INSTALL_RPATH "$ORIGIN;${MESHIOPLUSPLUS_LIB_DIR}")
endif()

# Install dependent libraries (the meshio++ shared library and, when it was built against
# them, its own HDF5/netCDF dependencies are resolved by the dynamic linker at load time).
if(${INSTALL_TPL_LIBRARIES} MATCHES ON)
    if(MESHIOPLUSPLUS_CORE_LOCATION)
        if(WIN32)
            install(DIRECTORY "${MESHIOPLUSPLUS_LIB_DIR}/" DESTINATION libs FILES_MATCHING PATTERN "*meshioplusplus*.dll")
        elseif(APPLE)
            install(DIRECTORY "${MESHIOPLUSPLUS_LIB_DIR}/" DESTINATION libs FILES_MATCHING PATTERN "*meshioplusplus*.dylib")
        else()
            file(GLOB MESHIOPLUSPLUS_LIB_FILES "${MESHIOPLUSPLUS_LIB_DIR}/*meshioplusplus*.so*")
            install(FILES ${MESHIOPLUSPLUS_LIB_FILES} DESTINATION libs)
        endif()
    endif()
endif()

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosMeshioPlusPlusCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosMeshioPlusPlusApplication" PARENT_SCOPE)
