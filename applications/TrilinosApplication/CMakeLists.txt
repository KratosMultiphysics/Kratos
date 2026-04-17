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

# If manually compiling some dependencies like MUMPS, linking to Fortran may be required
OPTION ( TRILINOS_APPLICATION_LINK_GFORTRAN "Not linking Trilinos application with GFortran" OFF )

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

# Adding Tpetra definitions
if(HAVE_TPETRA)
    add_definitions(-DHAVE_TPETRA)
endif()

## Kratos tests sources. Enabled by default
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_TRILINOS_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
    if(NOT HAVE_TPETRA)
        list(FILTER KRATOS_TRILINOS_TEST_SOURCES EXCLUDE REGEX ".*experimental.*")
    else()
        # Keep all experimental tests out of unity batches to avoid Epetra/Tpetra alias and type collisions.
        set(KRATOS_TRILINOS_EXPERIMENTAL_TEST_SOURCES ${KRATOS_TRILINOS_TEST_SOURCES})
        list(FILTER KRATOS_TRILINOS_EXPERIMENTAL_TEST_SOURCES INCLUDE REGEX ".*[Ee]xperimental.*\\.cpp$")
        if(KRATOS_TRILINOS_EXPERIMENTAL_TEST_SOURCES)
            set_source_files_properties(${KRATOS_TRILINOS_EXPERIMENTAL_TEST_SOURCES} PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
        endif()
    endif()
    kratos_add_gtests(TARGET KratosTrilinosCore SOURCES "${KRATOS_TRILINOS_TEST_SOURCES}" USE_MPI)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## Kratos benchmark sources. Disabled by default
if(${KRATOS_BUILD_BENCHMARK} MATCHES ON)
    file(GLOB_RECURSE KRATOS_BENCHMARK_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/*.cpp
    )

    foreach(file ${KRATOS_BENCHMARK_SOURCES})
        get_filename_component(filename ${file} NAME_WE)
        add_executable(${filename} ${file})
        target_link_libraries(${filename} PUBLIC KratosTrilinosCore benchmark::benchmark)
        set_target_properties(${filename} PROPERTIES COMPILE_DEFINITIONS "KRATOS_BENCHMARK=IMPORT,API")
        install(TARGETS ${filename} DESTINATION benchmark)
    endforeach(file ${KRATOS_BENCHMARK_SOURCES})
endif(${KRATOS_BUILD_BENCHMARK} MATCHES ON)

# ###############################################################
## TrilinosApplication core library (C++ parts)
set(TRILINOS_LINK_LIBS KratosCore KratosMPICore ${TRILINOS_LIBRARIES})
if (${TRILINOS_APPLICATION_LINK_GFORTRAN})
    set(TRILINOS_LINK_LIBS ${TRILINOS_LINK_LIBS} -lgfortran)
endif()
add_library( KratosTrilinosCore SHARED ${KRATOS_TRILINOS_APPLICATION_CORE_SOURCES} )
target_link_libraries(KratosTrilinosCore PUBLIC ${TRILINOS_LINK_LIBS})
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

# Exclude from unity build
set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/linear_solver_trilinos.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
# Conflicts with the Epetra adapter headers and their side-effects, which are included in other TUs that can't be excluded from unity batches due to their own issues, so exclude this TU as well to avoid both issues.
if(HAVE_TPETRA)
    # Uses #define KRATOS_AMGCL_MPI / #undef around amgcl_solver_impl.hpp; macro side-effects and heavy Epetra adapter headers are incompatible with unity batches.
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/external_includes/amgcl_mpi_solver_impl.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
    # The only TU that includes MatrixMarket_Tpetra.hpp, whose mmio_Tpetra.h reuses the MM_IO_H guard already claimed by Kratos' mmio.h with an incompatible signature.
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/auxiliary_matrix_market.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
    # Defines anonymous-namespace helpers (ClearMatrix, SetToZero, TwoNorm, Mult, …) that collide with identically-named helpers in add_trilinos_space_to_python.cpp when merged, causing unresolved overloads.
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_trilinos_space_experimental_to_python.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
endif(HAVE_TPETRA)

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

# A list to store the final set of files to be installed (libraries + links)
set(TRILINOS_LIBRARIES_AND_LINKS "")

foreach(LIB_FILE IN LISTS TRILINOS_LIBRARIES)
    # Get the directory (e.g., /path/to/libs)
    get_filename_component(LIB_DIR "${LIB_FILE}" DIRECTORY)
    
    # Get the base name (e.g., libthing.so.1.2.3)
    get_filename_component(LIB_NAME "${LIB_FILE}" NAME)
    
    # Identify the common prefix for the links (e.g., "libthing.so")
    # This assumes the links are generally named with the base SO name.
    # e.g., for libthing.so.1.2.3, this should find libthing.so, libthing.so.1, etc.
    
    # 1. Try to find the common SONAME part (e.g., "libthing.so")
    # This is a bit of a heuristic, but works for most shared library naming schemes.
    string(REGEX MATCH ".*\\.so(\\.|$)" SONAME_PREFIX "${LIB_NAME}")
    
    if(SONAME_PREFIX)
        # Remove the trailing dot or end of string if present, leaving the SONAME base.
        string(REGEX REPLACE "(\\.so)(\\..*|$)" "\\1" BASE_SONAME "${LIB_NAME}")
    else()
        # Fallback if the naming is non-standard (use the whole name as base)
        set(BASE_SONAME "${LIB_NAME}")
    endif()

    # 2. Use GLOB to find all files in the directory starting with the common prefix
    # e.g., /path/to/libs/libthing.so* will find libthing.so, libthing.so.1, libthing.so.1.2.3
    file(GLOB RELATED_FILES "${LIB_DIR}/${BASE_SONAME}*")

    # 3. Append the discovered files to the master list
    list(APPEND TRILINOS_LIBRARIES_AND_LINKS ${RELATED_FILES})
endforeach()

# Clean up any potential duplicates
list(REMOVE_DUPLICATES TRILINOS_LIBRARIES_AND_LINKS)

# Now, install the full list of files (the libraries + their symbolic links)
message(STATUS "Trilinos libs to be installed: " ${TRILINOS_LIBRARIES_AND_LINKS})
install(FILES ${TRILINOS_LIBRARIES_AND_LINKS} DESTINATION libs )

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