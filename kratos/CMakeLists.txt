set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(${USE_TRIANGLE_NONFREE_TPL} MATCHES ON )
    include_directories( ${KRATOS_SOURCE_DIR}/external_libraries/triangle )
else(${USE_TRIANGLE_NONFREE_TPL} MATCHES ON )
    include_directories(${KRATOS_SOURCE_DIR}/external_libraries/delaunator-cpp/include)
endif(${USE_TRIANGLE_NONFREE_TPL} MATCHES ON )
include_directories( ${KRATOS_SOURCE_DIR}/external_libraries/tinyexpr )

SET(ListFolders sources conditions constraints containers elements factories geometries includes integration linear_solvers mappers modeler modified_shape_functions operations processes controllers response_functions solving_strategies spaces spatial_containers utilities)

## Kratos main source code
SET(LIST_OF_SOURCE_EXCEPTIONS
    # kratos_version is excluded because I separated it to another CU so you didn't have to recompile the whole core every time the version changed, just relink it. In order to not pollute with more .dll / .so a separate .o is created and then using during the linking of KratosCore.dll
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/kratos_version.cpp

    # TPL sources are excluded from
    ${CMAKE_CURRENT_SOURCE_DIR}/modeler/cad_tessellation_modeler.cpp
)

FOREACH(subFolder ${ListFolders})
    FILE(GLOB_RECURSE ${subFolder}Sources "${subFolder}/*.cpp" "${subFolder}/*.c") #all .cpp and .c

    # Auxiliar list
    SET(AUXILIAR_CORE_SOURCES "")
    FOREACH(file ${${subFolder}Sources})
        SET(EXCEPTION_CHECK 0)
        FOREACH(exception ${LIST_OF_SOURCE_EXCEPTIONS})
            IF(${file} STREQUAL ${exception})
                SET(EXCEPTION_CHECK 1)
            ENDIF(${file} STREQUAL ${exception})
        ENDFOREACH(exception ${LIST_OF_SOURCE_EXCEPTIONS})
        IF (${EXCEPTION_CHECK} EQUAL 0)
            LIST(APPEND AUXILIAR_CORE_SOURCES ${file})
        ENDIF(${EXCEPTION_CHECK} EQUAL 0)
    ENDFOREACH(file ${${subFolder}Sources})

    # Append to final list
    LIST(APPEND KRATOS_CORE_SOURCES ${AUXILIAR_CORE_SOURCES})
ENDFOREACH(subFolder ${ListFolders})

if(${USE_TRIANGLE_NONFREE_TPL} MATCHES ON )
    SET(KRATOS_TPL_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/modeler/cad_tessellation_modeler.cpp
    )
endif(${USE_TRIANGLE_NONFREE_TPL} MATCHES ON )

## Kratos python interface code
FILE(GLOB_RECURSE KRATOS_PYTHON_SOURCES "python/*.cpp") #all .cpp

## Kratos Version
set( KRATOS_VERSION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/kratos_version.cpp
)

## Kratos testing engine sources
file(GLOB_RECURSE KRATOS_CORE_TESTING_ENGINE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/testing/*.cpp)

## Kratos I/o sources
file(GLOB_RECURSE KRATOS_CORE_INPUT_OUTPUT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/input_output/*.cpp)
file(GLOB KRATOS_EXPRESSION_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/expression/*.cpp)

## Kratos tests sources. Enabled by default
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_TEST_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_utilities/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp
    )
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

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
add_library(KratosCore SHARED ${KRATOS_CORE_SOURCES} ${KRATOS_CORE_TESTING_ENGINE_SOURCES} ${KRATOS_CORE_INPUT_OUTPUT_SOURCES} ${KRATOS_EXPRESSION_SOURCES} ${KRATOS_TEST_SOURCES} ${KRATOS_TPL_SOURCES} $<TARGET_OBJECTS:KratosVersion>)

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
pybind11_add_module(Kratos MODULE THIN_LTO ${KRATOS_PYTHON_SOURCES})
target_link_libraries(Kratos PRIVATE ${KRATOS_TEST_LIBRARIES} KratosCore)

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

    ## this cannot be compiled with unity build because calculate_discontinuous_distance_to_skin_process.cpp does implicit instantiation before specialization
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/processes/find_global_nodal_entity_neighbours_process.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)

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
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS Kratos DESTINATION libs )
install(TARGETS KratosCore DESTINATION libs)

# Mpi
if(${USE_MPI} MATCHES ON )
  add_subdirectory(mpi)
endif(${USE_MPI} MATCHES ON )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};Kratos" PARENT_SCOPE)
