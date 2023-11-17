cmake_minimum_required(VERSION 3.1)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
    cmake_policy(VERSION 3.12)
endif()

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to RelWithDebInfo")
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type")
endif()

project(VexCL)

set(VEXCL_MASTER_PROJECT OFF)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(VEXCL_MASTER_PROJECT ON)
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#----------------------------------------------------------------------------
# Compile-time options
#----------------------------------------------------------------------------
option(VEXCL_SHOW_KERNELS "Show generated kernels in tests and examples" OFF)
option(VEXCL_CACHE_KERNELS "Cache compiled kernels offline" ON)
option(VEXCL_SHOW_COPIES "Log vector copies to stdout for debugging purposes" OFF)
option(VEXCL_AMD_SI_WORKAROUND "Implement workaround for AMD SI GPUs" OFF)
set(VEXCL_CHECK_SIZES 0 CACHE STRING "Check that expressions have correct sizes")

#----------------------------------------------------------------------------
# Installation options
#----------------------------------------------------------------------------
option(VEXCL_INSTALL_CL_HPP "Install the OpenCL C++ header provided by VexCL" OFF)

#----------------------------------------------------------------------------
# Find Boost
#----------------------------------------------------------------------------
option(Boost_USE_STATIC_LIBS "Use static versions of Boost libraries" OFF)
if (WIN32)
    set(Boost_USE_STATIC_LIBS ON)
endif ()

find_package(Boost REQUIRED COMPONENTS
    chrono
    date_time
    filesystem
    program_options
    system
    thread
    unit_test_framework
    )

# Ensure all targets are available
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/VexCLBoostTargets.cmake")

#----------------------------------------------------------------------------
# Generic target
#----------------------------------------------------------------------------
add_library(Common INTERFACE)
add_library(VexCL::Common ALIAS Common)

target_compile_features(Common INTERFACE
    cxx_auto_type
    cxx_nullptr
    cxx_rvalue_references
    cxx_right_angle_brackets
    cxx_static_assert
    cxx_variadic_templates
    cxx_decltype
)

if (VEXCL_AMD_SI_WORKAROUND)
    target_compile_definitions(Common INTERFACE VEXCL_AMD_SI_WORKAROUND)
endif()

target_include_directories(Common INTERFACE
    ${Boost_INCLUDE_DIRS}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )


target_link_libraries(Common INTERFACE
    Boost::filesystem
    Boost::system
    Boost::thread
    )

if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_link_libraries(Common INTERFACE
        Boost::chrono
        Boost::date_time
        )
endif()

target_compile_options(Common INTERFACE
    # g++
    $<$<CXX_COMPILER_ID:GNU>:$<BUILD_INTERFACE:-Wall>>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-missing-braces>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-deprecated-declarations>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-ignored-attributes>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-unused-local-typedefs>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-variadic-macros>
    # Clang
    $<$<CXX_COMPILER_ID:Clang>:$<BUILD_INTERFACE:-Wall>>
    $<$<CXX_COMPILER_ID:Clang>:-Wno-missing-braces>
    $<$<CXX_COMPILER_ID:Clang>:-Wno-deprecated-declarations>
    $<$<CXX_COMPILER_ID:Clang>:-Wno-ignored-attributes>
    # MSVC
    $<$<CXX_COMPILER_ID:MSVC>:/bigobj>
    $<$<CXX_COMPILER_ID:MSVC>:/wd4003>
    $<$<CXX_COMPILER_ID:MSVC>:/wd4996>
    )

target_compile_definitions(Common INTERFACE
    # MSVC
    $<$<CXX_COMPILER_ID:MSVC>:NOMINMAX>
    $<$<CXX_COMPILER_ID:MSVC>:_VARIADIC_MAX=10>
    )

#----------------------------------------------------------------------------
# Find VexCL backends
#----------------------------------------------------------------------------
find_package(OpenCL)

# Ensure all targets are available
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/VexCLOpenCLTarget.cmake")

if(OpenCL_FOUND)
    set(CL_HPP_TARGET_OPENCL_VERSION 120 CACHE STRING "Target OpenCL version")
    set(CL_HPP_MINIMUM_OPENCL_VERSION 120 CACHE STRING "Minimum OpenCL version")

    add_library(OpenCL INTERFACE)
    add_library(VexCL::OpenCL ALIAS OpenCL)

    target_link_libraries(OpenCL INTERFACE Common OpenCL::OpenCL)
    target_compile_definitions(OpenCL INTERFACE VEXCL_BACKEND_OPENCL)
    target_compile_definitions(OpenCL INTERFACE CL_TARGET_OPENCL_VERSION=${CL_HPP_TARGET_OPENCL_VERSION})
    target_compile_definitions(OpenCL INTERFACE CL_HPP_TARGET_OPENCL_VERSION=${CL_HPP_TARGET_OPENCL_VERSION})
    target_compile_definitions(OpenCL INTERFACE CL_HPP_MINIMUM_OPENCL_VERSION=${CL_HPP_MINIMUM_OPENCL_VERSION})

    find_file(OPENCL_HPP CL/opencl.hpp HINTS ${OpenCL_INCLUDE_DIRS})
    message(STATUS " -- ${OPENCL_HPP} --")
    if(OPENCL_HPP)
        target_compile_definitions(OpenCL INTERFACE VEXCL_HAVE_OPENCL_HPP)
    endif()

    target_compile_options(Common INTERFACE
        $<$<CXX_COMPILER_ID:GNU>:-Wno-catch-value>
        )

    message(STATUS "Found VexCL::OpenCL")

    if ("${Boost_VERSION}" VERSION_GREATER_EQUAL "1.61.0")
        add_library(Compute INTERFACE)
        add_library(VexCL::Compute ALIAS Compute)

        target_link_libraries(Compute INTERFACE Common OpenCL::OpenCL)
        target_compile_definitions(Compute INTERFACE VEXCL_BACKEND_COMPUTE)

        message(STATUS "Found VexCL::Compute")
    endif()
endif()

find_package(CUDA)
if(CUDA_FOUND)
    add_library(CUDA INTERFACE)
    add_library(VexCL::CUDA ALIAS CUDA)

    target_include_directories(CUDA INTERFACE "${CUDA_INCLUDE_DIRS}")
    target_link_libraries(CUDA INTERFACE Common "${CUDA_CUDA_LIBRARY}")
    target_compile_definitions(CUDA INTERFACE VEXCL_BACKEND_CUDA)

    message(STATUS "Found VexCL::CUDA")
endif()

find_path(Boost_DLL NAMES boost/dll PATHS ${Boost_INCLUDE_DIRS})
if (Boost_DLL)
    if(OpenCL_INCLUDE_DIR)
        add_library(JIT INTERFACE)
        target_include_directories(JIT INTERFACE "${OpenCL_INCLUDE_DIR}")
    else()
        include(CheckIncludeFile)
        check_include_file("CL/cl_platform.h" HAVE_OpenCL_PLATFORM_H)
        if(HAVE_OpenCL_PLATFORM_H)
            add_library(JIT INTERFACE)
        endif()
    endif()

    if(NOT TARGET JIT)
        message(WARNING "The JIT interface requires OpenCL headers to be available."
                        "You can download them from https://github.com/KhronosGroup/OpenCL-Headers"
                        "Set OpenCL_INCLUDE_DIR to the location of the headers."
                        "For now, disabling the JIT target.")
    endif()
endif()


if(TARGET JIT)
    add_library(VexCL::JIT ALIAS JIT)

    set(VEXCL_JIT_COMPILER_FLAGS "" CACHE STRING "VexCL JIT compiler flags")
    target_compile_definitions(JIT INTERFACE VEXCL_JIT_COMPILER_FLAGS=${VEXCL_JIT_COMPILER_FLAGS})

    find_package(OpenMP)

    # Have to check several OPENMP_FOUND due to bug in
    # one version of CMake and the docs (fixed in patch release)
    # OpenMP is missing on macOS llvm default, for example
    if(OpenMP_FOUND OR OPENMP_FOUND OR OpenMP_CXX_FOUND)

        # CMake 3.9 FindOpenMP allows correct linking with Clang in more cases
        if(TARGET OpenMP::OpenMP_CXX)
            target_link_libraries(JIT INTERFACE OpenMP::OpenMP_CXX Common)
        else()
            # Clang may need -fopenmp=libiomp5 instead, can't be detected here without CMake 3.9
            target_link_libraries(JIT INTERFACE
                $<$<CXX_COMPILER_ID:GNU>:${OpenMP_CXX_FLAGS}>
                $<$<CXX_COMPILER_ID:Clang>:${OpenMP_CXX_FLAGS}>
                $<$<CXX_COMPILER_ID:Intel>:${OpenMP_CXX_FLAGS}>
                )
            target_compile_options(JIT INTERFACE ${OpenMP_CXX_FLAGS})
        endif()

        set(VEXCL_OMP_FLAGS "${OpenMP_CXX_FLAGS}")

        # We only need to add libraries to link to if this is using a preprocessor only OpenMP flag
        if("${OpenMP_CXX_FLAGS}" MATCHES ".*X(clang|preprocessor).*")
            foreach(item ${OpenMP_CXX_LIBRARIES})
                set(VEXCL_OMP_FLAGS "${VEXCL_OMP_FLAGS} ${item}")
            endforeach()
        endif()

        # Pass the required flags to code
        target_compile_definitions(JIT INTERFACE VEXCL_OMP_FLAGS=${VEXCL_OMP_FLAGS})
    endif()

    target_link_libraries(JIT INTERFACE Common ${CMAKE_DL_LIBS})
    target_compile_definitions(JIT INTERFACE VEXCL_BACKEND_JIT)

    message(STATUS "Found VexCL::JIT")
endif()

include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/VexCLTools.cmake")

#----------------------------------------------------------------------------
if (VEXCL_MASTER_PROJECT)
    option(VEXCL_BUILD_TESTS    OFF)
    option(VEXCL_BUILD_EXAMPLES OFF)

    foreach(target OpenCL CUDA Compute JIT)
        if (TARGET ${target})
            set(VEXCL_BACKEND "${target}" CACHE STRING "Select VexCL backend (OpenCL/CUDA/Compute/JIT/All)")
            break()
        endif()
    endforeach()
    message(STATUS "Selected backend: ${VEXCL_BACKEND}")

    set_property(CACHE VEXCL_BACKEND PROPERTY STRINGS "All" "OpenCL" "CUDA" "Compute" "JIT")

    if("${VEXCL_BACKEND}" STREQUAL "OpenCL")
        add_library(VexCL::Backend ALIAS OpenCL)
    elseif("${VEXCL_BACKEND}" STREQUAL "Compute")
        add_library(VexCL::Backend ALIAS Compute)
    elseif("${VEXCL_BACKEND}" STREQUAL "CUDA")
        add_library(VexCL::Backend ALIAS CUDA)
    elseif("${VEXCL_BACKEND}" STREQUAL "JIT")
        add_library(VexCL::Backend ALIAS JIT)
    endif()

    #------------------------------------------------------------------------
    # Interoperation with Boost.compute
    #------------------------------------------------------------------------
    option(VEXCL_HAVE_BOOST_COMPUTE "Use Boost.Compute algorithms" OFF)
    if (VEXCL_HAVE_BOOST_COMPUTE)
        find_path(BOOST_COMPUTE_INCLUDE boost/compute.hpp)

        add_library(compute_target INTERFACE)
        target_include_directories(compute_target INTERFACE ${BOOST_COMPUTE_INCLUDE})
        target_compile_definitions(compute_target INTERFACE VEXCL_HAVE_BOOST_COMPUTE)
    endif ()

    #------------------------------------------------------------------------
    # Interoperation with clogs
    #------------------------------------------------------------------------
    option(VEXCL_CLOGS "Use clogs algorithms" OFF)
    if (VEXCL_CLOGS)
        find_path(CLOGS_INCLUDE clogs/clogs.h)
        find_library(CLOGS_LIBRARY clogs)

        add_library(clogs_target INTERFACE)
        target_include_directories(clogs_target INTERFACE ${CLOGS_INCLUDE})
        target_compile_definitions(clogs_target INTERFACE VEXCL_HAVE_CLOGS)
        target_link_libraries(clogs_target INTERFACE ${CLOGS_LIBRARY})
    endif ()

    if (VEXCL_BUILD_TESTS)
        enable_testing()
        add_subdirectory(tests)
    endif()

    if (VEXCL_BUILD_EXAMPLES)
        add_subdirectory(examples)
    endif()

    add_subdirectory(docs)

    install(DIRECTORY vexcl DESTINATION include)

    install(TARGETS Common EXPORT VexCLTargets)

    if (TARGET VexCL::OpenCL)
        if(VEXCL_INSTALL_CL_HPP)
            install(DIRECTORY CL DESTINATION include)
        endif()
        install(TARGETS OpenCL EXPORT VexCLTargets)
    endif()

    if (TARGET VexCL::Compute)
        install(TARGETS Compute EXPORT VexCLTargets)
    endif()

    if (TARGET VexCL::CUDA)
        install(TARGETS CUDA EXPORT VexCLTargets)
    endif()

    if (TARGET VexCL::JIT)
        install(TARGETS JIT EXPORT VexCLTargets)
    endif()

    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/VexCLConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake/VexCLConfig.cmake"
        @ONLY
        )

    # Copies needed so that VexCLConfig can find these files
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/VexCLTools.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake/VexCLTools.cmake"
        COPYONLY
        )

    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/VexCLOpenCLTarget.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake/VexCLOpenCLTarget.cmake"
        COPYONLY
        )

    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/VexCLBoostTargets.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake/VexCLBoostTargets.cmake"
        COPYONLY
        )

    export(EXPORT VexCLTargets
        FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/VexCLTargets.cmake"
        NAMESPACE VexCL::
        )

    export(PACKAGE VexCL)

    install(EXPORT VexCLTargets
        FILE VexCLTargets.cmake
        NAMESPACE VexCL::
        DESTINATION share/vexcl/cmake
        )

    install(
      FILES
        ${CMAKE_CURRENT_BINARY_DIR}/cmake/VexCLConfig.cmake
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/VexCLTools.cmake
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/VexCLOpenCLTarget.cmake
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/VexCLBoostTargets.cmake
      DESTINATION
        share/vexcl/cmake
        )
endif()
