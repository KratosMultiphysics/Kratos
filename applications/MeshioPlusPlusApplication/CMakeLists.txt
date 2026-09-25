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
#         -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release
#   cmake --build build && cmake --install build --prefix <prefix>
#
# Upgrading from a meshio++ older than v16.0.0 requires rebuilding this application once: the
# C++ variants' SOVERSION is the ABI version, so the library this application needs is
# libmeshioplusplus_core_kratos.so.16 rather than an earlier .so.N.
#
# then point this configure at it with -Dmeshioplusplus_DIR=<prefix>/lib/cmake/meshioplusplus
#
# The pin is on the ABI version, not the release version. meshio++ carries
# MESHIOPLUSPLUS_ABI_VERSION, a binary-compatibility counter that moves only when the installed
# headers stop being compatible with an already-compiled consumer - so a release that provably
# cannot affect this application (v9.3.0 added one inline constexpr string) needs no rebuild,
# where a "X.Y.Z EXACT" pin would have forced one: v12.0.0 is a pure version-number bump over
# v11.6.0 and holds the ABI at 13, as v10.20.0 did over v10.19.0 at ABI 11.
#
# ABI 6 -> 11 is, unlike 3 -> 6, a run this application genuinely feels. Where ABI 5 and 6 came
# from format side-channel structs (MedInfo, OpenFoamInfo) it never names, 7 through 11 are all
# option/metadata structs it now passes by value: RefineOptions gained mRecordHierarchy (7),
# RemeshOptions grew twice (8, 9), SmoothMethod gained an explicit : std::uint8_t underlying
# type (10) and MeshMetadata gained the provenance fields (11, sizeof 256 -> 288).
#
# ABI 11 -> 13 touches this application far less: 12 (v10.35.0) is a Tier B ODR bump with no
# struct this application names changing size - NDArray::Size()'s inline body used to report 0
# for a rank-0 (scalar) array, so a 0-d field_data value was silently destroyed by every
# operation that clones a mesh; a plain relink against the new .so fixes it, no source change
# needed. 13 (v11.4.0) grows OpenFoamInfo (96 -> 128 bytes) with a new mRegion field for
# multi-region case selection, a struct this application does not name at ABI 11's level of
# integration (it reads OpenFOAM through the registry map, not the direct read_openfoam(...,
# OpenFoamInfo&) overload), so the layout change is inert until that overload is adopted.
#
# ABI 13 -> 16 is a rebuild, not a relink. 14 (v14.0.0) grows ReadOptions 56 -> 72 bytes with
# mPiece/mPieceSet, the piece selector of partitioned files; 15 (v15.0.0) puts mGhosts (a
# GhostPolicy) into its tail padding, so the size holds but a consumer built against v14 headers
# leaves that byte indeterminate. This application builds ReadOptions by value for every read,
# so both reach it. 16 (v16.0.0) appends CellType::Triangle7 before CellType::Custom, moving
# Custom 76 -> 77 and growing the inline cell-type tables; nothing here names either value, but
# the inlined tables are compiled into this application. Kratos has no seven-node triangle, so
# a triangle7 block (a Code_Aster TRIA7, say) still cannot become a ModelPart entity.
#
# v16.1.0 -> v16.13.0 holds ABI 16: every installed-header change is an addition (the release
# rows of meshio++'s doc/abi_reviews.md), so nothing already compiled here is affected.
#
# A single ABI rather than a window, deliberately. Supporting several meant a
# MESHIOPLUSPLUS_VERSION_AT_LEAST guard per feature, and the operations this application now
# exposes (remesh v10.10.0, optimize_volume v10.14.0, estimate_error v10.2.0, the gid format
# v10.18-10.20, compute_curvature v10.37.0, repair/shrinkwrap/sobolev_deform v10.38.0,
# vtkhdf v14.0.0, pvtu/pvtp/pvd v15.0.0, pcd/xyz v15.1.0, lsdyna v15.2.0, frd v15.3.0,
# compute_normals/gltf v15.4.0, tensor_invariants/binary OpenFOAM v15.5.0, nastran_h5 v15.7.0,
# code_aster v16.0.0, mphbin v16.1.0) would have needed one each. ABI 16 is v16.0.0 onwards, so
# everything compiles unconditionally.
#
# There is one release floor on top of the ABI pin, because the format table and a few direct
# overloads this application calls only exist from later ABI-16 releases: elmer, febio, femap,
# libmesh, mfem, patran, z88 and the results formats (abaqus_fil, ansys_rst, lsdyna_d3plot, ...)
# arrived between v16.2.0 and v16.13.0, as did read_mfem(path, grid functions, ReadOptions),
# write_mfem(path, mesh, bool), read_patran(path, results) and read_z88/write_z88's flags. An
# older ABI-16 install is binary compatible but does not *have* those symbols, so it is refused
# here by release rather than failing at link time.
#
# Getting it wrong is not silent: the C++ variants' SOVERSION is the ABI version, and every
# consumer translation unit carries a link-time sentinel naming it, so a header/library skew is
# a link error rather than memory corruption. See meshio++'s doc/abi.md for the criterion.
find_package(meshioplusplus CONFIG REQUIRED COMPONENTS CXX)

if(NOT MESHIOPLUSPLUS_ABI_VERSION EQUAL 16)
    message(FATAL_ERROR
        "MeshioPlusPlusApplication needs meshio++ ABI 16, but the install at "
        "'${meshioplusplus_DIR}' reports ABI ${MESHIOPLUSPLUS_ABI_VERSION} "
        "(release ${meshioplusplus_VERSION}). ABI 16 is meshio++ v16.0.0 onwards, "
        "which includes v16.13.0.")
endif()

if(meshioplusplus_VERSION VERSION_LESS 16.13.0)
    message(FATAL_ERROR
        "MeshioPlusPlusApplication needs meshio++ v16.13.0 or newer (for its format table and "
        "the elmer/mfem/patran/z88 entry points), but the install at '${meshioplusplus_DIR}' "
        "is release ${meshioplusplus_VERSION}.")
endif()

if(NOT TARGET meshioplusplus::core_kratos)
    message(FATAL_ERROR
        "MeshioPlusPlusApplication requires the meshio++ KRATOS mesh backend, but the install at "
        "'${meshioplusplus_DIR}' does not carry it (it shipped: ${MESHIOPLUSPLUS_MESH_BACKENDS}). "
        "Reconfigure meshio++ with -DMESHIOPLUSPLUS_INSTALL_CPP_BACKENDS=\"KRATOS\" (or a list "
        "including it).")
endif()

message(STATUS "meshio++ ${meshioplusplus_VERSION} (ABI ${MESHIOPLUSPLUS_ABI_VERSION}): "
               "backends '${MESHIOPLUSPLUS_MESH_BACKENDS}', "
               "parallel '${MESHIOPLUSPLUS_PARALLEL_BACKEND}', HDF5 '${MESHIOPLUSPLUS_WITH_HDF5}', "
               "netCDF '${MESHIOPLUSPLUS_WITH_NETCDF}', ADIOS2 '${MESHIOPLUSPLUS_WITH_ADIOS2}', "
               "TecIO '${MESHIOPLUSPLUS_WITH_TECIO}', bzip2 '${MESHIOPLUSPLUS_WITH_BZIP2}'")

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
