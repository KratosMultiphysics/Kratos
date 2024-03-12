set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMeshingApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${KRATOS_SOURCE_DIR}/external_libraries/triangle )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt )
include_directories( ${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication )

##add subdirectories for compilation of dependent libraries
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt)
set(KRATOS_MESHING_APPLICATION_EXTRA_LIBS tetMeshOpt)

if(${USE_TETGEN_NONFREE_TPL} MATCHES ON)
    include_directories(${TETGEN_INCLUDE})
    add_library(tetgen STATIC ${TETGEN_INCLUDE}/tetgen.cxx ${TETGEN_INCLUDE}/predicates.cxx)
    target_compile_definitions(tetgen PRIVATE -DTETLIBRARY)
endif(${USE_TETGEN_NONFREE_TPL} MATCHES ON)

## generate variables with the sources
file( GLOB KRATOS_MESHING_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/meshing_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/meshing_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/external_includes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_includes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/*.cpp
)

file(GLOB_RECURSE KRATOS_MESHING_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB KRATOS_MESHING_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

if(${INCLUDE_MMG} MATCHES ON)
    file( GLOB KRATOS_MESHING_APPLICATION_CORE ${KRATOS_MESHING_APPLICATION_CORE}
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mmg/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/mmg/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/mmg/*.cpp
    )
    if(${KRATOS_BUILD_TESTING} MATCHES ON)
        file(GLOB_RECURSE KRATOS_MESHING_APPLICATION_TESTING_SOURCES ${KRATOS_MESHING_APPLICATION_TESTING_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/mmg/*.cpp)
    endif(${KRATOS_BUILD_TESTING} MATCHES ON)

    IF (MMG_ROOT)
        set(MMG_BUILDDIR ${MMG_ROOT})
        set(MMG_DIR ${MMG_ROOT})
        set(MMG_LIBDIR ${MMG_ROOT}/lib)
        set(MMG_INCDIR ${MMG_ROOT}/include)
    ENDIF(MMG_ROOT)

    MESSAGE(STATUS "Searching for MMG in MMG_BUILDDIR: ${MMG_BUILDDIR}")
    MESSAGE(STATUS "... MMG_LIBDIR: ${MMG_LIBDIR}")
    MESSAGE(STATUS "... MMG_INCDIR: ${MMG_INCDIR}")

    find_package(MMG)
    find_package(MMG2D)
    find_package(MMG3D)
    find_package(MMGS)

    include_directories(${MMG_INCLUDE_DIRS} ${MMG2D_INCLUDE_DIRS} ${MMG3D_INCLUDE_DIRS} ${MMGS_INCLUDE_DIRS})
    set(ALL_MMG_LIBRARIES ${MMG_LIBRARIES} ${MMG2D_LIBRARIES} ${MMG3D_LIBRARIES} ${MMGS_LIBRARIES})
    set(KRATOS_MESHING_APPLICATION_EXTRA_LIBS ${KRATOS_MESHING_APPLICATION_EXTRA_LIBS} ${ALL_MMG_LIBRARIES})
    FOREACH (FILE ${ALL_MMG_LIBRARIES})
        GET_FILENAME_COMPONENT(RESOLVED_FILE_LIB "${FILE}" REALPATH)
        LIST(APPEND RESOLVED_FILES_LIBS "${RESOLVED_FILE_LIB}")

        # Adding other missing libs (that may be linked)
        GET_FILENAME_COMPONENT(FILE_PATH ${FILE} DIRECTORY)
        GET_FILENAME_COMPONENT(FILE_NAME ${FILE} NAME)
        FILE(GLOB LOCATION_FILE "${FILE_PATH}/${FILE_NAME}*")
        FOREACH (AUX_FILE ${LOCATION_FILE})
            LIST(FIND RESOLVED_FILES_LIBS ${AUX_FILE} FOUND_AUX)
            IF (${AUX_FILE} STREQUAL ${FILE})
                SET(FOUND_AUX 0)
            ENDIF (${AUX_FILE} STREQUAL ${FILE})
            IF (${FOUND_AUX} EQUAL -1)
                LIST(APPEND RESOLVED_FILES_LIBS "${AUX_FILE}")
            ENDIF (${FOUND_AUX} EQUAL -1)
        ENDFOREACH (AUX_FILE ${LOCATION_FILE})

        INSTALL(FILES ${RESOLVED_FILES_LIBS} DESTINATION ${CMAKE_INSTALL_PREFIX}/libs  COMPONENT libraries)
    ENDFOREACH (FILE ${ALL_MMG_LIBRARIES})

    message(STATUS "Compilation with MMG: ${MMG_LIBRARIES} ${MMG_INCLUDE_DIRS}")
    message(STATUS "Compilation with MMG2D: ${MMG2D_LIBRARIES} ${MMG2D_INCLUDE_DIRS}")
    message(STATUS "Compilation with MMG3D: ${MMG3D_LIBRARIES} ${MMG3D_INCLUDE_DIRS}")
    message(STATUS "Compilation with MMGS: ${MMGS_LIBRARIES} ${MMGS_INCLUDE_DIRS}")

    # Add definitions
    add_definitions( -DINCLUDE_MMG )
endif(${INCLUDE_MMG} MATCHES ON)

if(${INCLUDE_PMMG} MATCHES ON)
    IF(${USE_MPI} MATCHES OFF )
        message( FATAL_ERROR "Please compile Kratos in MPI to use ParMmg by setting the flag -DUSE_MPI=ON")
    ENDIF(${USE_MPI} MATCHES OFF )

    file( GLOB KRATOS_MESHING_APPLICATION_CORE ${KRATOS_MESHING_APPLICATION_CORE}
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/parmmg/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/parmmg/*.cpp
    )

    IF(NOT DEFINED INCLUDE_MMG)
        message( FATAL_ERROR "MMG not included. Please add -DINCLUDE_MMG=ON to the configure file in order to compile with ParMMG")
    ELSE(NOT DEFINED INCLUDE_MMG)
        if(${INCLUDE_MMG} MATCHES OFF)
            message( FATAL_ERROR " MMG disabled. Please set -DINCLUDE_MMG=ON in the configure file in order to compile with ParMMG")
        endif(${INCLUDE_MMG} MATCHES OFF)
    ENDIF(NOT DEFINED INCLUDE_MMG)

    IF(NOT DEFINED MMG3D_LIBRARIES)
        message( FATAL_ERROR "MMG3D not found. Please make sure MMG3D libraries are added setting the MMG build path as MMG_ROOT")
    ENDIF(NOT DEFINED MMG3D_LIBRARIES)

    add_definitions( -DINCLUDE_PMMG )

    IF(NOT DEFINED PMMG_LIBRARY_PREFIX)
        SET(PMMG_LIBRARY_PREFIX "lib")
    ENDIF(NOT DEFINED PMMG_LIBRARY_PREFIX)

    IF (PMMG_ROOT)
        IF(NOT DEFINED PMMG_LIBRARY_EXTENSION)
            IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
                IF(EXISTS "${PMMG_ROOT}/lib/${PMMG_LIBRARY_PREFIX}parmmg.dll") # Dynamic library
                    SET(PMMG_LIBRARY_EXTENSION "dll")
                ELSE(EXISTS "${PMMG_ROOT}/lib/${PMMG_LIBRARY_PREFIX}parmmg.dll") # Static library
                    SET(PMMG_LIBRARY_EXTENSION "lib")
                ENDIF(EXISTS "${PMMG_ROOT}/lib/${PMMG_LIBRARY_PREFIX}parmmg.dll")
            ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
                IF(EXISTS "${PMMG_ROOT}/lib/${PMMG_LIBRARY_PREFIX}parmmg.so") # Dynamic library
                    SET(PMMG_LIBRARY_EXTENSION "so")
                ELSE(EXISTS "${PMMG_ROOT}/lib/${PMMG_LIBRARY_PREFIX}parmmg.so") # Static library
                    SET(PMMG_LIBRARY_EXTENSION "a")
                ENDIF(EXISTS "${PMMG_ROOT}/lib/${PMMG_LIBRARY_PREFIX}parmmg.so")
            ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        ENDIF(NOT DEFINED PMMG_LIBRARY_EXTENSION)

        IF(NOT DEFINED PMMG_INCLUDE_DIR)
            SET(PMMG_INCLUDE_DIR ${PMMG_ROOT}/include/)
        ENDIF(NOT DEFINED PMMG_INCLUDE_DIR)
        IF(NOT DEFINED PMMG_LIBRARY)
            SET(PMMG_LIBRARY ${PMMG_ROOT}/lib/${PMMG_LIBRARY_PREFIX}parmmg.${PMMG_LIBRARY_EXTENSION})
        ENDIF(NOT DEFINED PMMG_LIBRARY)
    ENDIF (PMMG_ROOT)

    include_directories(${PMMG_INCLUDE_DIR})
    set(KRATOS_MESHING_APPLICATION_EXTRA_LIBS ${KRATOS_MESHING_APPLICATION_EXTRA_LIBS} KratosMPICore ${PMMG_LIBRARY})
    set(ALL_PMMG_LIBRARIES ${PMMG_LIBRARY})
    FOREACH (FILE ${ALL_PMMG_LIBRARIES})
        GET_FILENAME_COMPONENT(RESOLVED_FILE_LIB "${FILE}" REALPATH)
        LIST(APPEND RESOLVED_FILES_LIBS "${RESOLVED_FILE_LIB}")

        # Adding other missing libs (that may be linked)
        GET_FILENAME_COMPONENT(FILE_PATH ${FILE} DIRECTORY)
        GET_FILENAME_COMPONENT(FILE_NAME ${FILE} NAME)
        FILE(GLOB LOCATION_FILE "${FILE_PATH}/${FILE_NAME}*")
        FOREACH (AUX_FILE ${LOCATION_FILE})
            LIST(FIND RESOLVED_FILES_LIBS ${AUX_FILE} FOUND_AUX)
            IF (${AUX_FILE} STREQUAL ${FILE})
                SET(FOUND_AUX 0)
            ENDIF (${AUX_FILE} STREQUAL ${FILE})
            IF (${FOUND_AUX} EQUAL -1)
                LIST(APPEND RESOLVED_FILES_LIBS "${AUX_FILE}")
            ENDIF (${FOUND_AUX} EQUAL -1)
        ENDFOREACH (AUX_FILE ${LOCATION_FILE})

        INSTALL(FILES ${RESOLVED_FILES_LIBS} DESTINATION ${CMAKE_INSTALL_PREFIX}/libs  COMPONENT libraries)
    ENDFOREACH (FILE ${ALL_PMMG_LIBRARIES})

    message("Compilation with PMMG = " ${PMMG_LIBRARY})
    message("... and MMG3D = " ${MMG3D_LIBRARIES})
endif(${INCLUDE_PMMG} MATCHES ON)

add_library(KratosMeshingCore SHARED ${KRATOS_MESHING_APPLICATION_CORE} ${KRATOS_MESHING_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosMeshingCore PUBLIC KratosCore ${KRATOS_MESHING_APPLICATION_EXTRA_LIBS})

set_target_properties(KratosMeshingCore PROPERTIES COMPILE_DEFINITIONS "MESHING_APPLICATION=EXPORT,API")

if(${USE_TETGEN_NONFREE_TPL} MATCHES ON )
    target_link_libraries(KratosMeshingCore PUBLIC ${TETGEN_LIBRARIES})
endif(${USE_TETGEN_NONFREE_TPL} MATCHES ON )

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosMeshingApplication MODULE THIN_LTO ${KRATOS_MESHING_APPLICATION_SOURCES})
target_link_libraries(KratosMeshingApplication PRIVATE KratosMeshingCore)
set_target_properties(KratosMeshingApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosMeshingCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosMeshingApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosMeshingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosMeshingApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Unity build
if(CMAKE_UNITY_BUILD MATCHES ON)
    ## The MMG is not compatible with unity build
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mmg/mmg_utilities.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/parmmg/pmmg_utilities.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
endif(CMAKE_UNITY_BUILD MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/MeshingApplication.py KratosMultiphysics/MeshingApplication/__init__.py )
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_registry_lists.py KratosMultiphysics/MeshingApplication/python_registry_lists.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

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

# Install targets
install(TARGETS KratosMeshingCore DESTINATION libs )
install(TARGETS KratosMeshingApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosMeshingCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosMeshingApplication" PARENT_SCOPE)
