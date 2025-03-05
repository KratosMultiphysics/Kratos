cmake_minimum_required(VERSION 3.15...3.27)
project(fcpw)

option(FCPW_USE_ENOKI "Build enoki" ON)
option(FCPW_USE_EIGHT_WIDE_BRANCHING "Use 8 wide branching (default 4)" OFF)
option(FCPW_ENABLE_GPU_SUPPORT "Enable GPU support" OFF)
option(FCPW_BUILD_BINDINGS "Build bindings" OFF)
option(FCPW_BUILD_DEMO "Build demo" OFF)
option(FCPW_BUILD_TESTS "Build tests" OFF)

################################################################################
if(FCPW_BUILD_BINDINGS)
    if (CMAKE_VERSION VERSION_LESS 3.18)
        set(DEV_MODULE Development)
    else()
        set(DEV_MODULE Development.Module)
    endif()
    find_package(Python 3.8
        REQUIRED COMPONENTS Interpreter ${DEV_MODULE}
        OPTIONAL_COMPONENTS Development.SABIModule)
endif()

# submodule check
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if(CMAKE_BUILD_TYPE MATCHES RELEASE)
    add_definitions(-DNDEBUG)
endif()

################################################################################
# add dependencies
get_directory_property(hasParent PARENT_DIRECTORY)

if(hasParent)
    set(FCPW_EIGEN_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/deps/eigen-git-mirror PARENT_SCOPE)
else()
    set(FCPW_EIGEN_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/deps/eigen-git-mirror)
endif()

if(FCPW_ENABLE_GPU_SUPPORT)
    if(hasParent)
        set(FCPW_SLANG_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/deps/slang/include PARENT_SCOPE)
    else()
        set(FCPW_SLANG_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/deps/slang/include)
    endif()

    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_SYSTEM_NAME STREQUAL "Windows")
        message(STATUS "Operating system: ${CMAKE_SYSTEM_NAME}")
    else()
        # unfortunately macOS is not currently officially supported by gfx
        message(FATAL_ERROR "Unsupported operating system: ${CMAKE_SYSTEM_NAME}")
    endif()

    # find the .lib/.so files
    set(SLANG_LIBRARY_PATH ${CMAKE_CURRENT_SOURCE_DIR}/deps/slang/lib)
    find_library(FCPW_SLANG_LIBRARY NAMES slang PATHS ${SLANG_LIBRARY_PATH})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        find_library(FCPW_SLANG_GLSLANG_LIBRARY NAMES slang-glslang PATHS ${SLANG_LIBRARY_PATH})
    else()
        set(FCPW_SLANG_GLSLANG_LIBRARY "")
    endif()
    find_library(FCPW_GFX_LIBRARY NAMES gfx PATHS ${SLANG_LIBRARY_PATH})

    # find the corresponding .dll files
    if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
        set(SLANG_BINARY_PATH ${CMAKE_CURRENT_SOURCE_DIR}/deps/slang/bin)
        find_file(FCPW_SLANG_BINARY NAMES slang.dll PATHS ${SLANG_BINARY_PATH})
        find_file(FCPW_SLANG_GLSLANG_BINARY NAMES slang-glslang.dll PATHS ${SLANG_BINARY_PATH})
        find_file(FCPW_GFX_BINARY NAMES gfx.dll PATHS ${SLANG_BINARY_PATH})
    endif()

    message(STATUS "SLANG LIBRARY: " ${FCPW_SLANG_LIBRARY})
    message(STATUS "SLANG GLSLANG LIBRARY: " ${FCPW_SLANG_GLSLANG_LIBRARY})
    message(STATUS "GFX LIBRARY: " ${FCPW_GFX_LIBRARY})
endif()

if(FCPW_BUILD_TESTS)
    set(TBB_BUILD_SHARED          ON  CACHE BOOL " " FORCE)
    set(TBB_BUILD_STATIC          OFF CACHE BOOL " " FORCE)
    set(TBB_BUILD_TESTS           OFF CACHE BOOL " " FORCE)
    set(TBB_BUILD_TBBMALLOC       OFF CACHE BOOL " " FORCE)
    set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
    add_subdirectory(deps/tbb)
    set_property(TARGET tbb tbb_def_files PROPERTY FOLDER "dependencies")
    if(hasParent)
        set(FCPW_TESTS_TBB_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/deps/tbb/include PARENT_SCOPE)
    else()
        set(FCPW_TESTS_TBB_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/deps/tbb/include)
    endif()
endif()

if(FCPW_USE_ENOKI)
    add_subdirectory(deps/enoki)
    if(hasParent)
        set(FCPW_ENOKI_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/deps/enoki/include PARENT_SCOPE)
    else()
        set(FCPW_ENOKI_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/deps/enoki/include)
    endif()

    # Update the compilation flags
    enoki_set_compile_flags()
    enoki_set_native_flags()
endif()

################################################################################
# build library
add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)

if(FCPW_USE_ENOKI)
    target_compile_definitions(${PROJECT_NAME} INTERFACE -DFCPW_USE_ENOKI)

    # define SIMD width
    string(TOUPPER "${ENOKI_ARCH_FLAGS}" ENOKI_ARCH_FLAGS_UPPER)
    message(STATUS "Enoki Max ISA: " ${ENOKI_ARCH_FLAGS_UPPER})
    if(${ENOKI_ARCH_FLAGS_UPPER} MATCHES "SSE")
        target_compile_definitions(${PROJECT_NAME} INTERFACE -DFCPW_SIMD_WIDTH=4)
    elseif(${ENOKI_ARCH_FLAGS_UPPER} MATCHES "AVX2")
        target_compile_definitions(${PROJECT_NAME} INTERFACE -DFCPW_SIMD_WIDTH=8)
    elseif(${ENOKI_ARCH_FLAGS_UPPER} MATCHES "AVX")
        target_compile_definitions(${PROJECT_NAME} INTERFACE -DFCPW_SIMD_WIDTH=8)
    elseif(${ENOKI_ARCH_FLAGS_UPPER} MATCHES "KNL")
        target_compile_definitions(${PROJECT_NAME} INTERFACE -DFCPW_SIMD_WIDTH=16)
    elseif(${ENOKI_ARCH_FLAGS_UPPER} MATCHES "SKX")
        target_compile_definitions(${PROJECT_NAME} INTERFACE -DFCPW_SIMD_WIDTH=16)
    else()
        target_compile_definitions(${PROJECT_NAME} INTERFACE -DFCPW_SIMD_WIDTH=4)
    endif()

    if(FCPW_USE_EIGHT_WIDE_BRANCHING)
        target_compile_definitions(${PROJECT_NAME} INTERFACE -DFCPW_USE_EIGHT_WIDE_BRANCHING)
    endif()
endif()

################################################################################
# build Python bindings
if(FCPW_BUILD_BINDINGS)
    add_subdirectory(deps/nanobind)
    add_subdirectory(python)
endif()

################################################################################
# build demo and tests
if(FCPW_BUILD_DEMO OR FCPW_BUILD_TESTS)
    set(FCPW_ARGS_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope/deps/args)
    add_subdirectory(deps/polyscope)
endif()

if(FCPW_BUILD_DEMO)
    add_subdirectory(demos)
endif()

if(FCPW_BUILD_TESTS)
    add_subdirectory(tests)
endif()