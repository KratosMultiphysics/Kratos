# Not including current dir (we want to force users to specify the mpi/ in the include path)
set(CMAKE_INCLUDE_CURRENT_DIR OFF)

message("**** configuring Kratos to be compiled with MPI-Support ****")

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## Kratos MPI extension sources
file(GLOB_RECURSE KRATOS_MPI_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/spatial_containers/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/*.cpp
)

## Kratos MPI extension Python interface
file(GLOB_RECURSE KRATOS_MPI_PYTHON_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/python/*.cpp)

## Kratos tests sources. Enabled by default
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_MPI_TEST_UTILITIES_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_utilities/*.cpp)
    file(GLOB_RECURSE KRATOS_MPI_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)

    add_library(KratosMPICoreTestUtilities SHARED ${KRATOS_MPI_TEST_UTILITIES_SOURCES})
    add_executable(KratosMPICoreTest ${KRATOS_MPI_TEST_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/testing/mpi_testing.cpp)

    target_link_libraries(KratosMPICoreTestUtilities PUBLIC KratosMPICore KratosCoreTestUtilities)
    target_link_libraries(KratosMPICoreTest PUBLIC KratosMPICoreTestUtilities)

    # Discover mechanism for mpi tests in not implemented in gtest. We can wrapp the execution in a custom command
    # but installing and calling it would be best for this case.

    # gtest_discover_tests(KratosMPICoreTest)
    install(TARGETS KratosMPICoreTestUtilities DESTINATION libs)
    install(TARGETS KratosMPICoreTest DESTINATION test)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################################

## Define library for the C++ layer
add_library(KratosMPICore SHARED ${KRATOS_MPI_SOURCES})
target_link_libraries(KratosMPICore PUBLIC KratosCore ${MPI_LIBRARIES} )
#set_target_properties(KratosMPICore PROPERTIES COMPILE_DEFINITIONS "KRATOS_MPI_CORE=IMPORT,API")

## Define library for the Python interface layer
pybind11_add_module(KratosMPI MODULE THIN_LTO ${KRATOS_MPI_PYTHON_SOURCES})
# add_library(Kratos SHARED ${KRATOS_PYTHON_SOURCES})
target_link_libraries(KratosMPI PRIVATE KratosMPICore KratosCore )
# set_target_properties(KratosMPI PROPERTIES PREFIX "")

# Install library
install(TARGETS KratosMPICore DESTINATION libs)
install(TARGETS KratosMPI DESTINATION libs)

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosMPICore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosMPI" PARENT_SCOPE)

# Install python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/python_scripts/mpi_module_init.py" DESTINATION "KratosMultiphysics/mpi" RENAME "__init__.py" )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

if(${INSTALL_TESTING_FILES} MATCHES ON )
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION kratos/mpi
    PATTERN "*.git" EXCLUDE
    PATTERN "*.c" EXCLUDE
    PATTERN "*.h" EXCLUDE
    PATTERN "*.cpp" EXCLUDE
    PATTERN "*.hpp" EXCLUDE
  )

  if(${KRATOS_USE_FUTURE} MATCHES ON)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/future/tests DESTINATION kratos/future/mpi
      PATTERN "*.git" EXCLUDE
      PATTERN "*.c" EXCLUDE
      PATTERN "*.h" EXCLUDE
      PATTERN "*.cpp" EXCLUDE
      PATTERN "*.hpp" EXCLUDE
    )
  endif(${KRATOS_USE_FUTURE} MATCHES ON)
endif(${INSTALL_TESTING_FILES} MATCHES ON)
