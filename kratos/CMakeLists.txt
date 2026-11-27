set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(${USE_TRIANGLE_NONFREE_TPL} MATCHES ON )
    include_directories( ${KRATOS_SOURCE_DIR}/external_libraries/triangle )
else(${USE_TRIANGLE_NONFREE_TPL} MATCHES ON )
    include_directories(${KRATOS_SOURCE_DIR}/external_libraries/delaunator-cpp/include)
endif(${USE_TRIANGLE_NONFREE_TPL} MATCHES ON )
include_directories( ${KRATOS_SOURCE_DIR}/external_libraries/tinyexpr )

SET(LIST_FOLDER sources conditions constraints containers elements factories geometries includes integration linear_solvers mappers modeler modified_shape_functions operations processes controllers response_functions solving_strategies spaces spatial_containers utilities tensor_adaptors)
SET(LIST_FUTURE future)
SET(LIST_LEGACY legacy)

## Kratos main source code
SET(LIST_OF_SOURCE_EXCEPTIONS
    # kratos_version is excluded because I separated it to another CU so you didn't have to recompile the whole core every time the version changed, just relink it. In order to not pollute with more .dll / .so a separate .o is created and then using during the linking of KratosCore.dll
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/kratos_version.cpp

    # TPL sources are excluded from
    ${CMAKE_CURRENT_SOURCE_DIR}/modeler/cad_tessellation_modeler.cpp
)

## Kratos python interface code
FILE(GLOB_RECURSE KRATOS_PYTHON_SOURCES "python/*.cpp")

kratos_add_sources(KRATOS_CORE_SOURCES SOURCE_DIRS "${LIST_FOLDER}" EXCEPTIONS "${LIST_OF_SOURCE_EXCEPTIONS}")
if(${KRATOS_USE_FUTURE} MATCHES ON)
    add_definitions( -DKRATOS_USE_FUTURE )
    message(WARNING "You are using sources from Kratos::Future. This is not recommended for stability.")
    ## TODO: Change this in the future (no pun intended) and make it part of the interface target in a different folder.
    FILE(GLOB_RECURSE KRATOS_FUTURE_PYTHON_SOURCES "future/python/*.cpp")
    FILE(GLOB_RECURSE KRATOS_FUTURE_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/future/tests/cpp_tests/*.cpp)
    kratos_add_sources(KRATOS_CORE_SOURCES SOURCE_DIRS "${LIST_FUTURE}" EXCEPTIONS "${KRATOS_FUTURE_PYTHON_SOURCES}" "${KRATOS_FUTURE_TEST_SOURCES}")
endif(${KRATOS_USE_FUTURE} MATCHES ON)
if(${KRATOS_USE_LEGACY} MATCHES ON)
    add_definitions( -DKRATOS_USE_LEGACY )
    message(WARNING "You are using sources from Kratos::Legacy. This will be removed soon.")
    ## TODO: Change this in the future and make it part of the interface target in a different folder.
    FILE(GLOB_RECURSE KRATOS_LEGACY_PYTHON_SOURCES "legacy/python/*.cpp")
    FILE(GLOB_RECURSE KRATOS_LEGACY_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/legacy/tests/cpp_tests/*.cpp)
    kratos_add_sources(KRATOS_CORE_SOURCES SOURCE_DIRS "${LIST_LEGACY}" EXCEPTIONS "${KRATOS_LEGACY_PYTHON_SOURCES}" "${KRATOS_LEGACY_TEST_SOURCES}")
endif(${KRATOS_USE_LEGACY} MATCHES ON)

if(${USE_TRIANGLE_NONFREE_TPL} MATCHES ON )
    SET(KRATOS_TPL_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/modeler/cad_tessellation_modeler.cpp
    )
endif(${USE_TRIANGLE_NONFREE_TPL} MATCHES ON )

## Kratos Version
set( KRATOS_VERSION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/kratos_version.cpp
)

## Kratos testing engine sources
file(GLOB_RECURSE KRATOS_CORE_TESTING_ENGINE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/testing/*.cpp)

## Kratos I/O sources
file(GLOB_RECURSE KRATOS_CORE_INPUT_OUTPUT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/input_output/*.cpp)
file(GLOB KRATOS_EXPRESSION_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/expression/*.cpp)
file(GLOB KRATOS_TENSOR_ADAPTOR_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tensor_adaptors/*.cpp)

## Kratos tests sources. Enabled by default
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    include(GoogleTest)

    file(GLOB_RECURSE KRATOS_TEST_UTILITIES_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_utilities/*.cpp
    )

    file(GLOB_RECURSE KRATOS_TEST_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/testing/*.cpp
    )

    add_library(KratosCoreTestUtilities SHARED ${KRATOS_TEST_UTILITIES_SOURCES})
    add_executable(KratosCoreTest ${KRATOS_TEST_SOURCES} ${KRATOS_FUTURE_TEST_SOURCES} ${KRATOS_LEGACY_TEST_SOURCES})

    target_link_libraries(KratosCoreTestUtilities KratosCore GTest::gtest GTest::gmock)
    set_target_properties(KratosCoreTestUtilities PROPERTIES COMPILE_DEFINITIONS "KRATOS_TEST_UTILS=IMPORT,API")
    target_link_libraries(KratosCoreTest PUBLIC KratosCoreTestUtilities GTest::gmock_main)
    set_target_properties(KratosCoreTest PROPERTIES COMPILE_DEFINITIONS "KRATOS_TEST_CORE=IMPORT,API")

    install(TARGETS KratosCoreTestUtilities DESTINATION libs)
    install(TARGETS KratosCoreTest DESTINATION test)

    gtest_discover_tests(KratosCoreTest DISCOVERY_MODE PRE_TEST)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## Kratos benchmark sources. Disabled by default
if(${KRATOS_BUILD_BENCHMARK} MATCHES ON)
    file(GLOB_RECURSE KRATOS_BENCHMARK_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/*.cpp
    )

    foreach(file ${KRATOS_BENCHMARK_SOURCES})
        get_filename_component(filename ${file} NAME_WE)
        add_executable(${filename} ${file})
        target_link_libraries(${filename} PUBLIC KratosCore benchmark::benchmark)
        set_target_properties(${filename} PROPERTIES COMPILE_DEFINITIONS "KRATOS_BENCHMARK=IMPORT,API")
        install(TARGETS ${filename} DESTINATION benchmark)
    endforeach(file ${KRATOS_BENCHMARK_SOURCES})
endif(${KRATOS_BUILD_BENCHMARK} MATCHES ON)

## Define KratosVersion object
add_library(KratosVersion OBJECT ${KRATOS_VERSION_SOURCES})
set_target_properties(KratosVersion PROPERTIES COMPILE_DEFINITIONS "KRATOS_VERSION=IMPORT,API")
target_compile_definitions(KratosVersion PRIVATE
    KRATOS_MAJOR_VERSION=${KratosMultiphysics_MAJOR_VERSION}
    KRATOS_MINOR_VERSION=${KratosMultiphysics_MINOR_VERSION}
    KRATOS_PATCH_VERSION="${KratosMultiphysics_PATCH_VERSION}"
    KRATOS_SHA1_NUMBER="${KratosMultiphysics_SHA1_NUMBER}"
    KRATOS_BRANCH_NAME="${KratosMultiphysics_BRANCH_NAME}"
    KRATOS_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
)

## Define library KratosCore to be included in all of the others
add_library(KratosCore SHARED ${KRATOS_CORE_SOURCES} ${KRATOS_CORE_INPUT_OUTPUT_SOURCES} ${KRATOS_EXPRESSION_SOURCES} ${KRATOS_TENSOR_ADAPTOR_SOURCES} ${KRATOS_TPL_SOURCES} $<TARGET_OBJECTS:KratosVersion>)

# Set the shared external libraries for kratos (gidpost, tyniexpr, json, etc...)
set(KRATOS_EXTERNAL_LIBRARIES gidpost tinyexpr clipper)

# Set the TPL libraries if enabled (triangle, tetgen, etc...)
if(${USE_TRIANGLE_NONFREE_TPL} MATCHES ON )
  set(KRATOS_EXTERNAL_LIBRARIES ${KRATOS_EXTERNAL_LIBRARIES} triangle)
endif(${USE_TRIANGLE_NONFREE_TPL} MATCHES ON )

# older versions of GCC implemented std::filesystem in a separate extention, hence explicit linking
# against this extension is required. See the notes at the end of the official documentation:
# https://en.cppreference.com/w/cpp/filesystem
if(${CMAKE_COMPILER_IS_GNUCXX})
    if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 9.1)
        set(KRATOS_EXTERNAL_LIBRARIES ${KRATOS_EXTERNAL_LIBRARIES} -lstdc++fs)
    endif()
endif()

if(KRATOS_USE_TBB)
    set(KRATOS_EXTERNAL_LIBRARIES ${KRATOS_EXTERNAL_LIBRARIES} TBB::tbb)
endif()
target_link_libraries(KratosCore PUBLIC ${KRATOS_EXTERNAL_LIBRARIES})

set_target_properties(KratosCore PROPERTIES COMPILE_DEFINITIONS "KRATOS_CORE=IMPORT,API")

# Set batch size in the unity build
if(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
endif(CMAKE_UNITY_BUILD MATCHES ON)

target_compile_definitions(KratosCore PRIVATE
    KRATOS_MAJOR_VERSION=${KratosMultiphysics_MAJOR_VERSION}
    KRATOS_MINOR_VERSION=${KratosMultiphysics_MINOR_VERSION}
)

if (${AMGCL_GPGPU} MATCHES ON)
    add_subdirectory(
        ${KRATOS_SOURCE_DIR}/external_libraries/vexcl
        ${KRATOS_BINARY_DIR}/external_libraries/vexcl)

    if ("${AMGCL_GPGPU_BACKEND}" STREQUAL "OpenCL" AND TARGET VexCL::OpenCL)
        target_link_libraries(KratosCore PUBLIC VexCL::OpenCL)
        add_definitions(-DAMGCL_GPGPU)
    elseif ("${AMGCL_GPGPU_BACKEND}" STREQUAL "CUDA" AND TARGET VexCL::CUDA)
        target_link_libraries(KratosCore PUBLIC VexCL::CUDA)
        add_definitions(-DAMGCL_GPGPU)
    else()
        message(WARNING "AMGCL GPGPU backend not found")
    endif()
endif()

## Define library Kratos which defines the basic python interface
pybind11_add_module(Kratos MODULE THIN_LTO ${KRATOS_PYTHON_SOURCES} ${KRATOS_FUTURE_PYTHON_SOURCES} ${KRATOS_LEGACY_PYTHON_SOURCES})
target_link_libraries(Kratos PRIVATE KratosCore)

# Precompile most demanding / widely included headers
if(KRATOS_USE_PCH MATCHES ON)
    target_precompile_headers(KratosCore PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:${CMAKE_CURRENT_SOURCE_DIR}/includes/model_part.h>")
endif(KRATOS_USE_PCH MATCHES ON)

# Pass python version flags ONLY TO THE INTERFACE
target_compile_definitions(Kratos PRIVATE
    PYTHON_VERSION_MAJOR=${PYTHON_VERSION_MAJOR}
    PYTHON_VERSION_MINOR=${PYTHON_VERSION_MINOR}
    PYTHON_VERSION_PATCH=${PYTHON_VERSION_PATCH}
)

if(CMAKE_UNITY_BUILD MATCHES ON)
    ## Kratos cannot be compiled with just one unity build as it takes way too much ram
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/sources/kratos_parameters.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/sources/memory_info.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/utilities/mortar_utilities.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/utilities/exact_mortar_segmentation_utility.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/testing/distributed_test_case.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)

    ## This is giving issues in Windows Unity build
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/sources/system_information.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)

    ## This cannot be compiled with unity build because calculate_discontinuous_distance_to_skin_process.cpp does implicit instantiation before specialization
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/processes/find_global_nodal_entity_neighbours_process.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/linear_solvers/linear_solver_ublas.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)

    ## The following are conflicting due to KratosComponents
    file(GLOB_RECURSE KRATOS_STRATEGIES_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/solving_strategies/*.cpp)
    FOREACH(exception ${KRATOS_STRATEGIES_SOURCES})
        set_source_files_properties (${exception} PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
    ENDFOREACH(exception ${KRATOS_STRATEGIES_SOURCES})

    ## This sources in particular take too much ram by their own
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/python/add_amgcl_solver_to_python.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/python/add_strategies_to_python.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/python/add_geometrical_utilities_to_python.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/python/add_other_utilities_to_python.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/python/add_variable_utils_to_python.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
endif(CMAKE_UNITY_BUILD MATCHES ON)

# Install python scripts
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/ )

# Kratos Testing. Install everything except sources
if(${INSTALL_TESTING_FILES} MATCHES ON )
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION kratos
        PATTERN "*.git" EXCLUDE
        PATTERN "*.c" EXCLUDE
        PATTERN "*.h" EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.hpp" EXCLUDE
    )

    if(${KRATOS_USE_FUTURE} MATCHES ON)
        install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/future/tests DESTINATION kratos/future
                PATTERN "*.git" EXCLUDE
                PATTERN "*.c" EXCLUDE
                PATTERN "*.h" EXCLUDE
                PATTERN "*.cpp" EXCLUDE
                PATTERN "*.hpp" EXCLUDE
            )
    endif(${KRATOS_USE_FUTURE} MATCHES ON)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS Kratos DESTINATION libs)
install(TARGETS KratosCore DESTINATION libs)

# Mpi
if(${USE_MPI} MATCHES ON )
  add_subdirectory(mpi)
endif(${USE_MPI} MATCHES ON )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};Kratos" PARENT_SCOPE)
