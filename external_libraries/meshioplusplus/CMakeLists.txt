# meshio++ C++ core, vendored for Kratos Multiphysics. See README.md for
# provenance and resync instructions. Built as a STATIC library that is
# linked into KratosCore (see kratos/CMakeLists.txt, KRATOS_EXTERNAL_LIBRARIES).
#
# The compile definitions below are PUBLIC on purpose: the consumers of this
# target (KratosCore, in particular kratos/input_output/meshioplusplus_io.cpp)
# include meshio++ headers and must see the same mesh-backend and MESHIOPLUSPLUS_HAS_*
# macros the library itself was compiled with.
set(CMAKE_INCLUDE_CURRENT_DIR ON)

file(GLOB_RECURSE MESHIOPLUSPLUS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/cpp/src/*.cpp)

add_library(meshioplusplus_core STATIC
    ${MESHIOPLUSPLUS_SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/third_party/pugixml/pugixml.cpp
)

set_target_properties(meshioplusplus_core PROPERTIES
    POSITION_INDEPENDENT_CODE ON        # linked into the shared KratosCore
    CXX_VISIBILITY_PRESET hidden        # keep meshio++/pugixml symbols out of the DSO surface
    VISIBILITY_INLINES_HIDDEN ON
)

target_compile_features(meshioplusplus_core PUBLIC cxx_std_20)

target_include_directories(meshioplusplus_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/third_party/pugixml
)

# Mesh backend: fixed to KRATOS for the Kratos integration (Mesh == KratosMesh).
# The pybind bindings of meshio++ (which would force the MESHIO backend) are not built.
target_compile_definitions(meshioplusplus_core PUBLIC MESHIOPLUSPLUS_MESH_BACKEND_KRATOS)

# --- zlib: the Kratos root CMakeLists guarantees ZLIB_LIBRARIES (system or
# --- external_libraries/zlib) and its include dirs before this subdirectory is added.
target_compile_definitions(meshioplusplus_core PUBLIC MESHIOPLUSPLUS_HAS_ZLIB)
if(ZLIB_INCLUDE_DIRS)
    target_include_directories(meshioplusplus_core PUBLIC ${ZLIB_INCLUDE_DIRS})
endif()
target_link_libraries(meshioplusplus_core PUBLIC ${ZLIB_LIBRARIES})
message(STATUS "meshio++ : zlib enabled (${ZLIB_LIBRARIES})")

# --- HDF5 (optional): enables the med, cgns, h5m, hmf formats and the HDF
# --- data path of xdmf. NOTE: find_package must run HERE - the root scope
# --- sets HDF5_FOUND to FALSE for gidpost before adding subdirectories.
option(KRATOS_MESHIOPLUSPLUS_HDF5 "Enable the HDF5-backed meshio++ formats (med, cgns, h5m, hmf, xdmf-HDF)" ON)
set(MESHIOPLUSPLUS_HDF5_ENABLED FALSE)
if(KRATOS_MESHIOPLUSPLUS_HDF5)
    find_package(HDF5 QUIET COMPONENTS C)
    if(NOT HDF5_FOUND)
        set(HDF5_PREFER_PARALLEL ON)
        find_package(HDF5 QUIET COMPONENTS C)
    endif()
    if(HDF5_FOUND AND HDF5_IS_PARALLEL)
        find_package(MPI QUIET COMPONENTS C)
        if(NOT MPI_C_FOUND)
            set(HDF5_FOUND FALSE)
        endif()
    endif()
    if(HDF5_FOUND)
        set(MESHIOPLUSPLUS_HDF5_ENABLED TRUE)
    endif()
endif()
if(MESHIOPLUSPLUS_HDF5_ENABLED)
    message(STATUS "meshio++ : HDF5 formats enabled (HDF5 ${HDF5_VERSION})")
    target_compile_definitions(meshioplusplus_core PUBLIC MESHIOPLUSPLUS_HAS_HDF5)
    target_include_directories(meshioplusplus_core PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(meshioplusplus_core PUBLIC ${HDF5_C_LIBRARIES})
    if(HDF5_IS_PARALLEL AND MPI_C_FOUND)
        target_include_directories(meshioplusplus_core PUBLIC ${MPI_C_INCLUDE_DIRS})
        target_link_libraries(meshioplusplus_core PUBLIC MPI::MPI_C)
    endif()
else()
    message(STATUS "meshio++ : HDF5 not found/disabled - med, cgns, h5m, hmf and xdmf-HDF compiled out")
endif()

# --- netCDF (optional): enables the exodus format.
option(KRATOS_MESHIOPLUSPLUS_NETCDF "Enable the netCDF-backed meshio++ formats (exodus)" ON)
set(MESHIOPLUSPLUS_NETCDF_ENABLED FALSE)
if(KRATOS_MESHIOPLUSPLUS_NETCDF)
    find_package(netCDF CONFIG QUIET)
    if(netCDF_FOUND)
        set(MESHIOPLUSPLUS_NETCDF_ENABLED TRUE)
    else()
        find_library(MESHIOPLUSPLUS_NETCDF_LIBRARY netcdf)
        find_path(MESHIOPLUSPLUS_NETCDF_INCLUDE_DIR netcdf.h)
        if(MESHIOPLUSPLUS_NETCDF_LIBRARY AND MESHIOPLUSPLUS_NETCDF_INCLUDE_DIR)
            set(MESHIOPLUSPLUS_NETCDF_ENABLED TRUE)
        endif()
    endif()
endif()
if(MESHIOPLUSPLUS_NETCDF_ENABLED)
    message(STATUS "meshio++ : netCDF formats enabled")
    target_compile_definitions(meshioplusplus_core PUBLIC MESHIOPLUSPLUS_HAS_NETCDF)
    if(netCDF_FOUND)
        target_link_libraries(meshioplusplus_core PUBLIC netCDF::netcdf)
    else()
        target_include_directories(meshioplusplus_core PUBLIC ${MESHIOPLUSPLUS_NETCDF_INCLUDE_DIR})
        target_link_libraries(meshioplusplus_core PUBLIC ${MESHIOPLUSPLUS_NETCDF_LIBRARY})
    endif()
else()
    message(STATUS "meshio++ : netCDF not found/disabled - exodus compiled out")
endif()

# --- Parallel backend: match Kratos (OpenMP by default, sequential fallback).
find_package(OpenMP QUIET COMPONENTS CXX)
if(OpenMP_CXX_FOUND)
    target_compile_definitions(meshioplusplus_core PUBLIC MESHIOPLUSPLUS_PARALLEL_OPENMP)
    target_link_libraries(meshioplusplus_core PUBLIC OpenMP::OpenMP_CXX)
    message(STATUS "meshio++ : parallel backend OPENMP")
else()
    target_compile_definitions(meshioplusplus_core PUBLIC MESHIOPLUSPLUS_PARALLEL_SEQ)
    message(STATUS "meshio++ : parallel backend SEQ")
endif()

message(STATUS "meshioplusplus : ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include")
