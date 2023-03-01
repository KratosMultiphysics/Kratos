set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** Configuring KratosTrilinosApplication ****")

include(pybind11Tools)

if(${TRILINOS_EXCLUDE_ML_SOLVER} MATCHES ON)
    add_definitions(-DTRILINOS_EXCLUDE_ML_SOLVER)
    message("**** WARNING: Manually disabled ML-Solver ****")
endif(${TRILINOS_EXCLUDE_ML_SOLVER} MATCHES ON)

if(${TRILINOS_EXCLUDE_AZTEC_SOLVER} MATCHES ON)
    add_definitions(-DTRILINOS_EXCLUDE_AZTEC_SOLVER)
    message("**** WARNING: Manually disabled Aztec-Solver ****")
endif(${TRILINOS_EXCLUDE_AZTEC_SOLVER} MATCHES ON)

if(${TRILINOS_EXCLUDE_AMESOS_SOLVER} MATCHES ON)
    add_definitions(-DTRILINOS_EXCLUDE_AMESOS_SOLVER)
    message("**** WARNING: Manually disabled Amesos-Solver ****")
endif(${TRILINOS_EXCLUDE_AMESOS_SOLVER} MATCHES ON)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

# Find trilinos solver package -- this is our own find package!!
IF(NOT TRILINOS_FOUND)
    MESSAGE( FATAL_ERROR "Could not find Trilinos or one of its packages. Please set the cmake var TRILINOS_ROOT or the environment vars TRILINOS_INCLUDE_DIR and TRILINOS_LIBRARY_DIR, note also that if in your system the library have a prefix, like \"libtrilinos_epetra.so\" instead of \"libepetra.so\" you can specify the prefix by the variable -DTRILINOS_LIBRARY_PREFIX=\"trilinos_\" "    )
ENDIF(NOT TRILINOS_FOUND)

# Amesos2 is currenlty under development and is disabled by default. This can be removed in the future to enable it by default.
OPTION ( TRILINOS_EXCLUDE_AMESOS2_SOLVER "Exclude Amesos2 Solver" ON )

if(${TRILINOS_EXCLUDE_AMESOS2_SOLVER} MATCHES ON)
    add_definitions(-DTRILINOS_EXCLUDE_AMESOS2_SOLVER)
    message("**** WARNING: Manually disabled Amesos2-Solver ****")
elseif(NOT HAVE_AMESOS2)
    add_definitions(-DTRILINOS_EXCLUDE_AMESOS2_SOLVER)
    message("**** WARNING: Amesos2-Solver could not be found and will be disabled ****")
endif(${TRILINOS_EXCLUDE_AMESOS2_SOLVER} MATCHES ON)

# Including the amgcl, which is header only
include_directories("../../external_libraries")

# CHecking MPI is included
IF (NOT ${USE_MPI} OR NOT DEFINED USE_MPI)
    MESSAGE(FATAL_ERROR "TrilinosApplication requires MPI, please ser ON USE_MPI")
ENDIF (NOT ${USE_MPI} OR NOT DEFINED USE_MPI)

# Sources for the Trilinos application core
file(GLOB_RECURSE KRATOS_TRILINOS_APPLICATION_CORE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/trilinos_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/external_includes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_factories/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

# Sources for the Python module
file(GLOB_RECURSE KRATOS_TRILINOS_PYTHON_INTERFACE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

## Kratos tests sources. Enabled by default
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_TRILINOS_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

# ###############################################################
## TrilinosApplication core library (C++ parts)
add_library( KratosTrilinosCore SHARED ${KRATOS_TRILINOS_APPLICATION_CORE_SOURCES} ${KRATOS_TRILINOS_TEST_SOURCES} )
target_link_libraries(KratosTrilinosCore PRIVATE KratosCore KratosMPICore ${TRILINOS_LIBRARIES} )
set_target_properties( KratosTrilinosCore PROPERTIES COMPILE_DEFINITIONS "TRILINOS_APPLICATION=EXPORT,API")

target_include_directories(KratosTrilinosCore SYSTEM PUBLIC ${TRILINOS_INCLUDE_DIR})

## TrilinosApplication python module
pybind11_add_module( KratosTrilinosApplication MODULE THIN_LTO ${KRATOS_TRILINOS_PYTHON_INTERFACE_SOURCES} )
target_link_libraries( KratosTrilinosApplication PRIVATE KratosTrilinosCore)
set_target_properties( KratosTrilinosApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosTrilinosCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosTrilinosApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# Changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosTrilinosApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# Changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosTrilinosApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/TrilinosApplication.py KratosMultiphysics/TrilinosApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# TODO: Move this to the new system. Not sure how it will work on windows but being Trilinos should give no problem.
message(STATUS "Trilinos libs to be installed: " ${TRILINOS_LIBRARIES})
install(FILES ${TRILINOS_LIBRARIES} DESTINATION libs )

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

## Install targets
install(TARGETS KratosTrilinosCore DESTINATION libs )
install(TARGETS KratosTrilinosApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosTrilinosCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosTrilinosApplication" PARENT_SCOPE)
