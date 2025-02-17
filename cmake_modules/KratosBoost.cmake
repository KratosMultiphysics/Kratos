# Finding and including the BOOST library (version should not matter anymore)

# Check if Boost_INCLUDE_DIRS and Boost_LIBRARY_DIRS are already defined
if(Boost_INCLUDE_DIRS AND Boost_LIBRARY_DIRS)
  message(STATUS "Using predefined Boost paths for Boost_INCLUDE_DIRS and Boost_LIBRARY_DIRS, instead of searching for Boost")
else(Boost_INCLUDE_DIRS AND Boost_LIBRARY_DIRS)
  # Check if the BOOST_ROOT environment variable is defined
  # If so, use it to set the BOOST_ROOT variable in CMake
  if(DEFINED ENV{BOOST_ROOT})
    set(BOOST_ROOT $ENV{BOOST_ROOT})
  endif(DEFINED ENV{BOOST_ROOT})

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
if(Boost_INCLUDE_DIRS)
  include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
endif(Boost_INCLUDE_DIRS)

# Display information about the found Boost installation
message(STATUS "Boost Include: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost Linkdir: ${Boost_LIBRARY_DIRS}")

# Check the Boost version and handle compatibility
if(Boost_VERSION_STRING VERSION_LESS 1.70)
  # Warn the user if Boost version is lower than 1.70
  message(WARNING "Kratos requires at least Boost version 1.70 to enable Boost-related performance improvements")
else(Boost_VERSION_STRING VERSION_LESS 1.70)
  # Define BOOST_UBLAS_MOVE_SEMANTICS for improved performance with Boost UBLAS in newer versions
  add_definitions(-DBOOST_UBLAS_MOVE_SEMANTICS)
endif(Boost_VERSION_STRING VERSION_LESS 1.70)
