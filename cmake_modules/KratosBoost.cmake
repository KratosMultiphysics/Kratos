# Finding and including the BOOST library (version should not matter anymore)

# Check if the BOOST_ROOT environment variable is defined
if(DEFINED ENV{BOOST_ROOT})
  set(BOOST_ROOT $ENV{BOOST_ROOT})
  # Attempt to find the Boost library using CMake's find_package
  find_package(Boost)

  # If Boost is found with the specified BOOST_ROOT, set the include and library directories
  if(Boost_FOUND)
    # Check that boost/any.hpp exists in the BOOST_ROOT include directory
    if(EXISTS "${BOOST_ROOT}/boost/any.hpp")
      set(BOOST_INCLUDEDIR "${Boost_INCLUDE_DIRS}")
      set(BOOST_LIBRARYDIR "${Boost_LIBRARY_DIRS}")
    endif(EXISTS "${BOOST_ROOT}/boost/any.hpp")
  else (Boost_FOUND)
    message(FATAL_ERROR "It was not possible to find Boost in your machine with BOOST_ROOT=${BOOST_ROOT}. Please define BOOST_ROOT correctly if you have it installed.")
  endif(Boost_FOUND)
else(DEFINED ENV{BOOST_ROOT})
  # Check if BOOST_INCLUDEDIR and BOOST_LIBRARYDIR are already defined as environment variables
  if(DEFINED ENV{BOOST_INCLUDEDIR})
  set(BOOST_INCLUDEDIR $ENV{BOOST_INCLUDEDIR})
  endif(DEFINED ENV{BOOST_INCLUDEDIR})

  if(DEFINED ENV{BOOST_LIBRARYDIR})
  set(BOOST_LIBRARYDIR $ENV{BOOST_LIBRARYDIR})
  endif(DEFINED ENV{BOOST_LIBRARYDIR})
endif(DEFINED ENV{BOOST_ROOT})

# Check if BOOST_INCLUDEDIR and BOOST_LIBRARYDIR are already defined
if(BOOST_INCLUDEDIR AND BOOST_LIBRARYDIR)
  message(STATUS "Using predefined Boost paths for BOOST_INCLUDEDIR and BOOST_LIBRARYDIR, instead of searching for Boost")

  # If BOOST_ROOT is not defined it is the parent directory of BOOST_INCLUDEDIR and BOOST_LIBRARYDIR
  if(NOT BOOST_ROOT)
    get_filename_component(BOOST_ROOT ${BOOST_INCLUDEDIR} DIRECTORY)
    MESSAGE(STATUS "BOOST_ROOT not defined, setting it to ${BOOST_ROOT}")
  endif(NOT BOOST_ROOT)
else(BOOST_INCLUDEDIR AND BOOST_LIBRARYDIR)
  # Attempt to find the Boost library using CMake's find_package
  find_package(Boost)

  # If Boost is not found, attempt to locate it in the external libraries folder
  if(NOT Boost_FOUND)
    # Set BOOST_ROOT to the expected external library directory within the project source
    set(BOOST_ROOT "${KRATOS_SOURCE_DIR}/external_libraries/boost")

    # Try finding Boost again using the newly set BOOST_ROOT
    find_package(Boost)

    # If Boost is still not found, terminate with an error message
    if(NOT Boost_FOUND)
      message(FATAL_ERROR "It was not possible to find Boost in your machine. Please define BOOST_ROOT correctly if you have it installed.")
    endif(NOT Boost_FOUND)
  endif(NOT Boost_FOUND)
endif()

# Configure Boost usage settings
set(Boost_USE_STATIC_LIBS   OFF) # Use shared (dynamic) libraries instead of static ones
set(Boost_USE_MULTITHREADED  ON) # Enable multithreading support in Boost
set(Boost_REALPATH           ON) # Resolve symbolic links to their real paths

# Ensure Boost include and library directories are properly set
if(BOOST_INCLUDEDIR)
  message(STATUS "Boost include directory manually defined")
elseif(Boost_INCLUDE_DIRS)
  message(STATUS "Boost include directory found with find_package")
  set(BOOST_INCLUDEDIR ${Boost_INCLUDE_DIRS})
else(BOOST_INCLUDEDIR)
  message(STATUS "Boost include deduced from root directory")
  set(BOOST_INCLUDEDIR ${BOOST_ROOT}/include)
endif(BOOST_INCLUDEDIR)
include_directories(SYSTEM ${BOOST_INCLUDEDIR})

if (BOOST_LIBRARYDIR)
  message(STATUS "Boost library directory manually defined")
elseif(Boost_LIBRARY_DIRS)
  message(STATUS "Boost library directory found with find_package")
  set(BOOST_LIBRARYDIR ${Boost_LIBRARY_DIRS})
else(BOOST_LIBRARYDIR)
  message(STATUS "Boost library deduced from root directory")
  set(BOOST_LIBRARYDIR ${BOOST_ROOT}/lib)
endif(BOOST_LIBRARYDIR)

# Display information about the found Boost installation
message(STATUS "Boost include directory: ${BOOST_INCLUDEDIR}")
message(STATUS "Boost library directory: ${BOOST_LIBRARYDIR}")

# Check the Boost version and handle compatibility
if(Boost_VERSION_STRING VERSION_LESS 1.70)
  # Warn the user if Boost version is lower than 1.70
  message(WARNING "Kratos requires at least Boost version 1.70 to enable Boost-related performance improvements")
else(Boost_VERSION_STRING VERSION_LESS 1.70)
  # Define BOOST_UBLAS_MOVE_SEMANTICS for improved performance with Boost UBLAS in newer versions
  add_definitions(-DBOOST_UBLAS_MOVE_SEMANTICS)
endif(Boost_VERSION_STRING VERSION_LESS 1.70)