# Finding and including the EIGEN library (header-only).
#
# Eigen is a header-only linear-algebra library. Kratos no longer ships a copy
# under external_libraries/eigen3; instead the include directory is resolved
# from one of the following sources, in order of priority:
#
#   1. EIGEN3_ROOT  - a user-provided root (env var or -DEIGEN3_ROOT=...),
#                     mirroring how BOOST_ROOT is handled in KratosBoost.cmake.
#   2. OS install   - a system-wide Eigen (e.g. Ubuntu/Debian "libeigen3-dev",
#                     which installs headers under /usr/include/eigen3).
#   3. Auto-download - as a last resort, the pinned version (KRATOS_EIGEN_VERSION)
#                     is downloaded and extracted into the build tree via
#                     FetchContent, similar to how gtest/benchmark are fetched.
#
# This module sets:
#   EIGEN3_INCLUDE_DIR - directory that contains the "Eigen" header folder
#
# and makes the include directory available globally (SYSTEM) so that both the
# core and any application can use Eigen. EIGEN_MPL2_ONLY is defined globally so
# every translation unit agrees on the licensed subset of Eigen.

# Version used for the automatic download fallback.
if(NOT DEFINED KRATOS_EIGEN_VERSION)
  set(KRATOS_EIGEN_VERSION "5.0.0")
endif()

# Check if the EIGEN3_ROOT environment variable is defined
if(DEFINED ENV{EIGEN3_ROOT})
  set(EIGEN3_ROOT $ENV{EIGEN3_ROOT})
endif()

# Search for the Eigen headers. The header "Eigen/src/Core/util/Macros.h" is a
# stable, version-independent marker that exists in every Eigen 3.x release.
# HINTS cover the user-provided root; PATHS cover common OS install locations.
find_path(EIGEN3_INCLUDE_DIR
  NAMES Eigen/src/Core/util/Macros.h
  HINTS
    ${EIGEN3_ROOT}
    $ENV{EIGEN3_ROOT}
  PATHS
    /usr/include/eigen3
    /usr/local/include/eigen3
    /opt/eigen3
    /opt/homebrew/include/eigen3
    /opt/homebrew/opt/eigen/include
)

if(EIGEN3_INCLUDE_DIR)
  if(EIGEN3_ROOT AND EXISTS "${EIGEN3_ROOT}/Eigen/src/Core/util/Macros.h")
    message(STATUS "Eigen include directory (from EIGEN3_ROOT): ${EIGEN3_INCLUDE_DIR}")
  else()
    message(STATUS "Eigen include directory (from OS install): ${EIGEN3_INCLUDE_DIR}")
  endif()
else()
  # Eigen was not found on the system. Download and extract the pinned version
  # into the build tree. Eigen is header-only, so we only download + extract and
  # do NOT run its CMakeLists (which would build tests/examples).
  #
  # This uses file(DOWNLOAD) + execute_process(tar) directly (the same pattern
  # as the tetgen download in the root CMakeLists.txt) instead of FetchContent,
  # because the deprecated single-argument FetchContent_Populate form is broken
  # on some CMake versions (e.g. 3.22) and mis-parses the URL as a list.
  message(STATUS "Eigen not found on the system. Downloading Eigen ${KRATOS_EIGEN_VERSION} ...")

  set(EIGEN3_DOWNLOAD_DIR "${CMAKE_CURRENT_BINARY_DIR}/external_libraries/eigen3")
  set(EIGEN3_ARCHIVE "${EIGEN3_DOWNLOAD_DIR}/eigen-${KRATOS_EIGEN_VERSION}.tar.gz")
  set(EIGEN3_URL "https://gitlab.com/libeigen/eigen/-/archive/${KRATOS_EIGEN_VERSION}/eigen-${KRATOS_EIGEN_VERSION}.tar.gz")

  file(MAKE_DIRECTORY "${EIGEN3_DOWNLOAD_DIR}")

  if(EXISTS "${EIGEN3_ARCHIVE}")
    message(STATUS "-- ${EIGEN3_ARCHIVE} found. unpacking...")
  else()
    message(STATUS "-- ${EIGEN3_ARCHIVE} not found, trying to download it from ${EIGEN3_URL}")
    file(DOWNLOAD "${EIGEN3_URL}" "${EIGEN3_ARCHIVE}" STATUS EIGEN3_DOWNLOAD_STATUS)
    list(GET EIGEN3_DOWNLOAD_STATUS 0 EIGEN3_DOWNLOAD_RESULT)
    if(NOT EIGEN3_DOWNLOAD_RESULT EQUAL 0)
      list(GET EIGEN3_DOWNLOAD_STATUS 1 EIGEN3_DOWNLOAD_ERROR)
      message(FATAL_ERROR "Failed to download Eigen ${KRATOS_EIGEN_VERSION} from ${EIGEN3_URL}: ${EIGEN3_DOWNLOAD_ERROR}")
    endif()
    message(STATUS "-- ${EIGEN3_ARCHIVE} downloaded. unpacking...")
  endif()

  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xzf "${EIGEN3_ARCHIVE}"
    WORKING_DIRECTORY "${EIGEN3_DOWNLOAD_DIR}"
    RESULT_VARIABLE EIGEN3_EXTRACT_RESULT
  )
  if(NOT EIGEN3_EXTRACT_RESULT EQUAL 0)
    message(FATAL_ERROR "Failed to extract ${EIGEN3_ARCHIVE}")
  endif()

  # The archive extracts into a top-level "eigen-<version>" folder.
  set(EIGEN3_INCLUDE_DIR "${EIGEN3_DOWNLOAD_DIR}/eigen-${KRATOS_EIGEN_VERSION}")
  if(NOT EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h")
    message(FATAL_ERROR "Eigen ${KRATOS_EIGEN_VERSION} was downloaded but the expected headers were not found at ${EIGEN3_INCLUDE_DIR}")
  endif()
  message(STATUS "Eigen include directory (downloaded): ${EIGEN3_INCLUDE_DIR}")
endif()

# Make the Eigen headers available globally (SYSTEM) so both the core and any
# application can use them, and pin the licensed subset for every translation unit.
include_directories( SYSTEM ${EIGEN3_INCLUDE_DIR} )
add_definitions( -DEIGEN_MPL2_ONLY )
