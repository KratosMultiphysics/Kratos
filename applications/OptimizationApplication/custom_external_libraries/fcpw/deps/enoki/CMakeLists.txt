cmake_minimum_required (VERSION 2.8.12)
project(enoki CXX)

if (POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW) # MACOSX_RPATH is enabled by default
endif()

set(CMAKE_MACOSX_RPATH ON)

option(ENOKI_CUDA     "Build Enoki CUDA library?" OFF)
option(ENOKI_AUTODIFF "Build Enoki automatic differentation library?" OFF)
option(ENOKI_PYTHON   "Build pybind11 interface to CUDA & automatic differentiation libraries?" OFF)

if (ENOKI_CUDA)
  set(ENOKI_CUDA_COMPUTE_CAPABILITY "50" CACHE STRING "Compute capability as specified by https://developer.nvidia.com/cuda-gpus")
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_${ENOKI_CUDA_COMPUTE_CAPABILITY},code=compute_${ENOKI_CUDA_COMPUTE_CAPABILITY}")

  if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.17.0")
    # CMake 3.17.0 introduces CMAKE_CUDA_RUNTIME which must be used instead of explicit -cudart shared flag
    set(CMAKE_CUDA_RUNTIME_LIBRARY Shared)
  else()
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -cudart shared")
  endif()

  if (NOT WIN32)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -fvisibility=hidden")
  endif()
  enable_language(CUDA)
  add_definitions(-DENOKI_CUDA=1)
endif()

set(ENOKI_MASTER_PROJECT OFF)
if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
  set(ENOKI_MASTER_PROJECT ON)
endif()

option(ENOKI_TEST "Build Enoki test suite?" OFF)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" ENOKI_U_CMAKE_BUILD_TYPE)

set(ENOKI_CXXFLAGS_BACKUP ${CMAKE_CXX_FLAGS})
set(ENOKI_NATIVE_FLAGS "")

if (MSVC)
  set(ENOKI_ARCH_FLAGS_FILE archflags_win32.cpp)
  if (CMAKE_SIZEOF_VOID_P EQUAL 4)
    message(WARNING "Enoki does not support vectorization on 32-bit Windows due to various")
    message(WARNING "platform limitations (unaligned stack, calling conventions don't allow")
    message(WARNING "passing vector registers, etc.). Switching to scalar mode.")
  endif()
else()
  set(ENOKI_ARCH_FLAGS_FILE archflags_unix.cpp)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  if (ENOKI_CUDA AND CMAKE_BUILD_TYPE MATCHES "Release")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DNDEBUG")
  endif()
endif()

# Clang refuses -march=native on ARM machines
if (${CMAKE_SYSTEM_PROCESSOR} MATCHES "armv7")
  set(ENOKI_ARCH_FLAGS -march=armv7-a -mtune=cortex-a7 -mfpu=neon-vfpv4 -mfloat-abi=hard -mfp16-format=ieee)
  set(ENOKI_NATIVE_FLAGS ${ENOKI_ARCH_FLAGS})
elseif (${CMAKE_SYSTEM_PROCESSOR} MATCHES "aarch64")
  set(ENOKI_ARCH_FLAGS -march=armv8-a+simd -mtune=cortex-a53)
  set(ENOKI_NATIVE_FLAGS ${ENOKI_ARCH_FLAGS})
elseif (CMAKE_CXX_COMPILER MATCHES "/em\\+\\+(-[a-zA-Z0-9.])?$")
  # Emscripten
else()
  if (UNIX)
    set(ENOKI_ARCH_FLAGS_LINK_LIBRARIES ${CMAKE_EXE_LINKER_FLAGS})
  endif()
  try_run(
    ENOKI_ARCH_FLAGS_RETVAL ENOKI_ARCH_FLAGS_COMPILE_RESULT
    ${CMAKE_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/resources/${ENOKI_ARCH_FLAGS_FILE}
    COMPILE_OUTPUT_VARIABLE ENOKI_ARCH_FLAGS_MSG
    RUN_OUTPUT_VARIABLE ENOKI_ARCH_FLAGS
    LINK_LIBRARIES ${ENOKI_ARCH_FLAGS_LINK_LIBRARIES}
  )

  if (NOT ${ENOKI_ARCH_FLAGS_COMPILE_RESULT})
    message(FATAL_ERROR "Failed to compile 'archflags' binary: ${ENOKI_ARCH_FLAGS_MSG}")
  endif()
endif()

set(CMAKE_CXX_FLAGS ${ENOKI_CXXFLAGS_BACKUP})

# Optimize for current architecture
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR
    CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  if (NOT ENOKI_NATIVE_FLAGS)
    set(ENOKI_NATIVE_FLAGS -march=native)
  endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(ENOKI_NATIVE_FLAGS -xHost)
elseif (MSVC)
  # MSVC doesn't have a flag equivalent to -march=native
  # The 'archflags' executable provides this functionality
  set(ENOKI_NATIVE_FLAGS ${ENOKI_ARCH_FLAGS})
endif()

if (NOT ENOKI_MASTER_PROJECT)
  set(ENOKI_ARCH_FLAGS   ${ENOKI_ARCH_FLAGS} PARENT_SCOPE)
  set(ENOKI_NATIVE_FLAGS ${ENOKI_NATIVE_FLAGS} PARENT_SCOPE)
endif()

set(ENOKI_HOST "INTEL")
set(ENOKI_TEST_SSE42 OFF)
set(ENOKI_TEST_AVX OFF)
set(ENOKI_TEST_AVX2 OFF)
set(ENOKI_TEST_KNL OFF)
set(ENOKI_TEST_SKX OFF)
set(ENOKI_TEST_ARM OFF)

string(TOUPPER "${ENOKI_ARCH_FLAGS}" ENOKI_ARCH_FLAGS)

if (MSVC OR ENOKI_ARCH_FLAGS MATCHES "SSE")
  set(ENOKI_TEST_SSE42 ON)
endif()

if (ENOKI_ARCH_FLAGS MATCHES "AVX")
  set(ENOKI_TEST_SSE42 ON)
  set(ENOKI_TEST_AVX ON)
endif()

if (ENOKI_ARCH_FLAGS MATCHES "AVX2")
  set(ENOKI_TEST_SSE42 ON)
  set(ENOKI_TEST_AVX ON)
  set(ENOKI_TEST_AVX2 ON)
endif()

if (ENOKI_ARCH_FLAGS MATCHES "KNL")
  set(ENOKI_TEST_SSE42 ON)
  set(ENOKI_TEST_AVX ON)
  set(ENOKI_TEST_AVX2 ON)
  set(ENOKI_TEST_KNL ON)
endif()

if (ENOKI_ARCH_FLAGS MATCHES "SKX")
  set(ENOKI_TEST_SSE42 ON)
  set(ENOKI_TEST_AVX ON)
  set(ENOKI_TEST_AVX2 ON)
  set(ENOKI_TEST_SKX ON)
endif()

if (ENOKI_ARCH_FLAGS MATCHES "ARMV")
  set(ENOKI_HOST "ARM")
  set(ENOKI_TEST_ARM ON)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/resources")

macro(enoki_set_native_flags)
  add_compile_options(${ENOKI_NATIVE_FLAGS})
  string(REPLACE ";" " " ENOKI_NATIVE_FLAGS_STR "${ENOKI_NATIVE_FLAGS}")
  # Some linkers want to know the architecture flags (for LTO)
  if (NOT MSVC)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ENOKI_NATIVE_FLAGS_STR}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ENOKI_NATIVE_FLAGS_STR}")
  endif()
endmacro()

macro(enoki_set_compile_flags)
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR
      CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
      CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    string(TOUPPER "${CMAKE_BUILD_TYPE}" ENOKI_U_CMAKE_BUILD_TYPE)
    if (NOT (${ENOKI_U_CMAKE_BUILD_TYPE} MATCHES "DEB"))
      if (NOT (${CMAKE_CXX_FLAGS} MATCHES "fsanitize"))
        # Don't use stack security features in release mode
        add_compile_options(-fno-stack-protector)

        # In release mode, don't keep the frame pointer in a dedicated register unless needed
        add_compile_options(-fomit-frame-pointer)
      endif()
    endif()

    # Never update the 'errno' variable due to arithmetic exceptions
    add_compile_options(-fno-math-errno)

    if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel")
      # Look for opportunities to fuse additions & multiplications into FMA operations
      add_compile_options(-ffp-contract=fast)
    endif()
  endif()

  # Disable overly aggressive FP optimization in the Intel compiler
  if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    add_compile_options(-fp-model precise)
  endif()

  if (MSVC)
    # Disable buffer security check cookie
    set(Configurations RELEASE RELWITHDEBINFO MINSIZEREL)
    foreach(Configuration ${Configurations})
      string(REPLACE "/GS"  "" CMAKE_CXX_FLAGS_${Configuration} "${CMAKE_CXX_FLAGS_${Configuration}}")
    endforeach()
    string(REPLACE "/GS"  "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    add_compile_options("$<$<NOT:$<CONFIG:Debug>>:/GS->")

    # Enable intrinsic functions
    add_compile_options("$<$<CONFIG:Release>:/Oi>")

    # Honor __forceinline statements even in debug mode, needed to avoid internal compiler errors on MSVC
    string(REPLACE "/Ob0" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
    add_compile_options("$<$<CONFIG:Debug>:/Ob1>")

    # Don't complain about class/struct X needs to have dll-interface to be used by clients of class/struct Y
    add_compile_options("/wd4251")
  endif()
endmacro()

include_directories(include)

set(ENOKI_HEADERS
    ${PROJECT_SOURCE_DIR}/include/enoki/array.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_avx.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_avx2.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_avx512.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_base.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_call.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_enum.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_fallbacks.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_generic.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_idiv.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_intrin.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_kmask.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_masked.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_macro.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_math.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_recursive.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_router.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_sse42.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_static.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_struct.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_traits.h
    ${PROJECT_SOURCE_DIR}/include/enoki/array_utils.h
    ${PROJECT_SOURCE_DIR}/include/enoki/autodiff.h
    ${PROJECT_SOURCE_DIR}/include/enoki/color.h
    ${PROJECT_SOURCE_DIR}/include/enoki/complex.h
    ${PROJECT_SOURCE_DIR}/include/enoki/dynamic.h
    ${PROJECT_SOURCE_DIR}/include/enoki/fwd.h
    ${PROJECT_SOURCE_DIR}/include/enoki/half.h
    ${PROJECT_SOURCE_DIR}/include/enoki/matrix.h
    ${PROJECT_SOURCE_DIR}/include/enoki/morton.h
    ${PROJECT_SOURCE_DIR}/include/enoki/python.h
    ${PROJECT_SOURCE_DIR}/include/enoki/quaternion.h
    ${PROJECT_SOURCE_DIR}/include/enoki/random.h
    ${PROJECT_SOURCE_DIR}/include/enoki/sh.h
    ${PROJECT_SOURCE_DIR}/include/enoki/special.h
    ${PROJECT_SOURCE_DIR}/include/enoki/stl.h
    ${PROJECT_SOURCE_DIR}/include/enoki/transform.h
)

include(CheckCXXSourceRuns)

macro(CHECK_CXX_COMPILER_AND_LINKER_FLAGS _RESULT _CXX_FLAGS _LINKER_FLAGS)
  set(CMAKE_REQUIRED_FLAGS ${_CXX_FLAGS})
  set(CMAKE_REQUIRED_LIBRARIES ${_LINKER_FLAGS})
  set(CMAKE_REQUIRED_QUIET TRUE)
  check_cxx_source_runs("#include <iostream>\nint main(int argc, char **argv) { std::cout << \"test\"; return 0; }" ${_RESULT})
  set(CMAKE_REQUIRED_FLAGS "")
  set(CMAKE_REQUIRED_LIBRARIES "")
endmacro()

# Prefer libc++ in conjunction with Clang
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND NOT CMAKE_CXX_FLAGS MATCHES "-stdlib=libc\\+\\+")
  CHECK_CXX_COMPILER_AND_LINKER_FLAGS(HAS_LIBCPP "-stdlib=libc++" "-stdlib=libc++")
  if (HAS_LIBCPP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -D_LIBCPP_VERSION")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -stdlib=libc++")
    message(STATUS "Enoki: using libc++.")
  else()
    CHECK_CXX_COMPILER_AND_LINKER_FLAGS(HAS_LIBCPP_AND_CPPABI "-stdlib=libc++" "-stdlib=libc++ -lc++abi")
    if (HAS_LIBCPP_AND_CPPABI)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -D_LIBCPP_VERSION")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ -lc++abi")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -stdlib=libc++ -lc++abi")
      message(STATUS "Enoki: using libc++ and libc++abi.")
    else()
      message(FATAL_ERROR "When Clang is used to compile Enoki, libc++ must be available -- GCC's libstdc++ is not supported! (please insteall the libc++ development headers, provided e.g. by the packages 'libc++-dev' and 'libc++abi-dev' on Debian/Ubuntu).")
    endif()
  endif()
endif()

if (ENOKI_TEST)
  enable_testing()
  add_subdirectory(tests)
endif()

add_definitions(-DENOKI_BUILD=1)

if (ENOKI_CUDA)
  if (CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 11)
    include_directories(ext/cub)
  endif()
  set(CMAKE_SHARED_LINKER_FLAGS_BACKUP "${CMAKE_SHARED_LINKER_FLAGS}")
  string(REPLACE "-stdlib=libc++" "" CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")
  add_library(enoki-cuda SHARED
      src/cuda/common.cuh
      src/cuda/common.cu
      src/cuda/horiz.cu
      src/cuda/jit.cu
  )
  target_compile_definitions(enoki-cuda PRIVATE -DENOKI_CUDA_COMPUTE_CAPABILITY=${ENOKI_CUDA_COMPUTE_CAPABILITY} -DTHRUST_IGNORE_CUB_VERSION_CHECK)
  target_link_libraries(enoki-cuda PRIVATE cuda)
  message(STATUS "Enoki: building the CUDA backend.")
endif()

enoki_set_native_flags()
enoki_set_compile_flags()

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
elseif (MSVC)
  add_compile_options("/std:c++17")
endif()

if (ENOKI_AUTODIFF)
  add_definitions(-DENOKI_AUTODIFF=1)
  add_library(enoki-autodiff SHARED
      ${PROJECT_SOURCE_DIR}/include/enoki/autodiff.h
      ${PROJECT_SOURCE_DIR}/src/autodiff/autodiff.cpp
  )
  target_compile_definitions(enoki-autodiff PRIVATE -DENOKI_AUTODIFF_BUILD=1)
  if (ENOKI_CUDA)
      target_link_libraries(enoki-autodiff PRIVATE enoki-cuda)
  endif()
  message(STATUS "Enoki: building the autodiff backend.")
endif()

if (ENOKI_PYTHON)
  set(ENOKI_PYBIND11_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ext/pybind11"
    CACHE STRING "Path containing the 'pybind11' library used to compile Enoki.")
  message(STATUS "Enoki: building the Python plugin.")

  set(CMAKE_CXX_STANDARD 17)
  add_subdirectory(${ENOKI_PYBIND11_DIR} pybind11)

  pybind11_add_module(
      enoki-python-core
      THIN_LTO
      src/python/main.cpp
  )

  pybind11_add_module(
      enoki-python-scalar
      THIN_LTO

      src/python/common.h
      src/python/complex.h
      src/python/matrix.h
      src/python/quat.h
      src/python/scalar.cpp
      src/python/scalar_0d.cpp
      src/python/scalar_1d.cpp
      src/python/scalar_2d.cpp
      src/python/scalar_3d.cpp
      src/python/scalar_4d.cpp
      src/python/scalar_complex.cpp
      src/python/scalar_matrix.cpp
      src/python/scalar_quat.cpp
      src/python/scalar_pcg32.cpp
  )

  pybind11_add_module(
      enoki-python-dynamic
      THIN_LTO

      src/python/common.h
      src/python/complex.h
      src/python/matrix.h
      src/python/quat.h
      src/python/dynamic.cpp
      src/python/dynamic_0d.cpp
      src/python/dynamic_1d.cpp
      src/python/dynamic_2d.cpp
      src/python/dynamic_3d.cpp
      src/python/dynamic_4d.cpp
      src/python/dynamic_complex.cpp
      src/python/dynamic_matrix.cpp
      src/python/dynamic_pcg32.cpp
  )

  set(ENOKI_PYTHON_TARGETS core scalar dynamic)

  if (ENOKI_CUDA)
    pybind11_add_module(
        enoki-python-cuda
        THIN_LTO

        src/python/common.h
        src/python/complex.h
        src/python/matrix.h
        src/python/quat.h
        src/python/cuda.cpp
        src/python/cuda_0d.cpp
        src/python/cuda_1d.cpp
        src/python/cuda_2d.cpp
        src/python/cuda_3d.cpp
        src/python/cuda_4d.cpp
        src/python/cuda_complex.cpp
        src/python/cuda_matrix.cpp
        src/python/cuda_pcg32.cpp
    )
    target_link_libraries(enoki-python-cuda PRIVATE enoki-cuda cuda)
    set(ENOKI_PYTHON_TARGETS ${ENOKI_PYTHON_TARGETS} cuda)
  endif()

  if (ENOKI_CUDA AND ENOKI_AUTODIFF)
    pybind11_add_module(
        enoki-python-cuda-autodiff
        THIN_LTO

        src/python/common.h
        src/python/complex.h
        src/python/matrix.h
        src/python/quat.h
        src/python/cuda_autodiff.cpp
        src/python/cuda_autodiff_0d.cpp
        src/python/cuda_autodiff_1d.cpp
        src/python/cuda_autodiff_2d.cpp
        src/python/cuda_autodiff_3d.cpp
        src/python/cuda_autodiff_4d.cpp
        src/python/cuda_autodiff_complex.cpp
        src/python/cuda_autodiff_matrix.cpp
    )
    target_link_libraries(enoki-python-cuda-autodiff PRIVATE enoki-autodiff enoki-cuda cuda)
    set(ENOKI_PYTHON_TARGETS ${ENOKI_PYTHON_TARGETS} cuda-autodiff)
  endif()

  foreach (TARGET IN LISTS ENOKI_PYTHON_TARGETS)
    string(REPLACE "-" "_" TARGET_U ${TARGET})

    if (${TARGET} MATCHES "core")
      set_target_properties(enoki-python-${TARGET} PROPERTIES
        OUTPUT_NAME "core"
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/enoki
      )
    else()
      set_target_properties(enoki-python-${TARGET} PROPERTIES
        OUTPUT_NAME ${TARGET_U}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/enoki
      )
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "GCC|Clang")
      target_compile_options(enoki-python-${TARGET} PRIVATE -g0 -Os)
    elseif (MSVC)
      target_compile_options(enoki-python-${TARGET} PRIVATE /Os)
    endif()
  endforeach()

  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/enoki/__init__.py
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/resources/__init__.py
    ${CMAKE_CURRENT_BINARY_DIR}/enoki/__init__.py
  )

  add_custom_target(
    enoki-python-init
    ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/enoki/__init__.py
  )
endif()

# Build the documentation
if (ENOKI_MASTER_PROJECT)
  find_package(Sphinx)

  if (Sphinx_FOUND)
    set(SPHINX_INPUT_DIR  "${CMAKE_CURRENT_SOURCE_DIR}/docs")
    set(SPHINX_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

    add_custom_target(mkdoc
        ${SPHINX_EXECUTABLE} -b html "${SPHINX_INPUT_DIR}" "${SPHINX_OUTPUT_DIR}"
        COMMENT "Building HTML documentation with Sphinx"
        USES_TERMINAL)
  endif()
endif()

mark_as_advanced(
  Eigen3_DIR PYBIND11_INSTALL PYBIND11_PYTHON_VERSION PYBIND11_TEST
  PYBIND11_WERROR USE_PYTHON_INCLUDE_DIR
)
