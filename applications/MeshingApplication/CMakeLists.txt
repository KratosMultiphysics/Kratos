set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMeshingApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${KRATOS_SOURCE_DIR}/external_libraries/triangle )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt )
include_directories( ${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication )

if(${USE_TETGEN_NONFREE_TPL} MATCHES ON )
    add_definitions( -DTETGEN143 )
    add_definitions( -DTETLIBRARY )
    add_definitions( -DUSE_TETGEN_NONFREE_TPL )
    include_directories(${TETGEN_INCLUDES})
endif(${USE_TETGEN_NONFREE_TPL} MATCHES ON )

##add subdirectories for compilation of dependent libraries
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt)

## generate variables with the sources
set( KRATOS_MESHING_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/meshing_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/meshing_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/external_includes/mesh_suites.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/local_refine_geometry_mesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/uniform_refinement_utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/meshing_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_includes/meshing_flags.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/metrics_levelset_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/internal_variables_interpolation_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/metric_fast_init_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/metrics_hessian_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/metrics_error_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/nodal_values_interpolation_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/multiscale_refining_process.cpp
)

set(KRATOS_MESHING_APPLICATION_CORE_PMMG "")
if(${INCLUDE_PMMG} MATCHES ON)
set( KRATOS_MESHING_APPLICATION_CORE_PMMG ${KRATOS_MESHING_APPLICATION_CORE_PMMG}
    #${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/pmmg_utilities.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/pmmg_process.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/custom_io/pmmg_io.cpp
)
endif(${INCLUDE_PMMG} MATCHES ON)


## If we have MMG support, use it
set(KRATOS_MESHING_APPLICATION_CORE_MMG "")
if(${INCLUDE_MMG} MATCHES ON)
set( KRATOS_MESHING_APPLICATION_CORE_MMG ${KRATOS_MESHING_APPLICATION_CORE_MMG}
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mmg_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/mmg_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/mmg_io.cpp
)
endif(${INCLUDE_MMG} MATCHES ON)

file(GLOB_RECURSE KRATOS_MESHING_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_MESHING_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_library(KratosMeshingCore SHARED ${KRATOS_MESHING_APPLICATION_CORE} ${KRATOS_MESHING_APPLICATION_CORE_MMG} ${KRATOS_MESHING_APPLICATION_CORE_PMMG} ${KRATOS_MESHING_APPLICATION_TESTING_SOURCES})

IF(NOT DEFINED MMG_LIBRARY_PREFIX)
    SET(MMG_LIBRARY_PREFIX "lib")
ENDIF(NOT DEFINED MMG_LIBRARY_PREFIX)

if(${INCLUDE_MMG} MATCHES ON)
    # Define includes and libraries with just one include
    IF (MMG_ROOT)
        IF(NOT DEFINED MMG_LIBRARY_EXTENSION)
            IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
                IF(EXISTS "${MMG_ROOT}/lib/${MMG_LIBRARY_PREFIX}mmg.dll") # Dynamic library
                    SET(MMG_LIBRARY_EXTENSION "dll")
                ELSE(EXISTS "${MMG_ROOT}/lib/${MMG_LIBRARY_PREFIX}mmg.dll") # Static library
                    SET(MMG_LIBRARY_EXTENSION "lib")
                ENDIF(EXISTS "${MMG_ROOT}/lib/${MMG_LIBRARY_PREFIX}mmg.dll")
            ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
                IF(EXISTS "${MMG_ROOT}/lib/${MMG_LIBRARY_PREFIX}mmg.so") # Dynamic library
                    SET(MMG_LIBRARY_EXTENSION "so")
                ELSE(EXISTS "${MMG_ROOT}/lib/${MMG_LIBRARY_PREFIX}mmg.so") # Static library
                    SET(MMG_LIBRARY_EXTENSION "a")
                ENDIF(EXISTS "${MMG_ROOT}/lib/${MMG_LIBRARY_PREFIX}mmg.so")
            ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        ENDIF(NOT DEFINED MMG_LIBRARY_EXTENSION)

        IF(NOT DEFINED MMG_INCLUDE_DIR)
            SET(MMG_INCLUDE_DIR ${MMG_ROOT}/include/)
        ENDIF(NOT DEFINED MMG_INCLUDE_DIR)
        IF(NOT DEFINED MMG2D_INCLUDE_DIR)
            SET(MMG2D_INCLUDE_DIR ${MMG_ROOT}/include/mmg/mmg2d/)
        ENDIF(NOT DEFINED MMG2D_INCLUDE_DIR)
        IF(NOT DEFINED MMG3D_INCLUDE_DIR)
            SET(MMG3D_INCLUDE_DIR ${MMG_ROOT}/include/mmg/mmg3d/)
        ENDIF(NOT DEFINED MMG3D_INCLUDE_DIR)
        IF(NOT DEFINED MMGS_INCLUDE_DIR)
            SET(MMGS_INCLUDE_DIR ${MMG_ROOT}/include/mmg/mmgs/)
        ENDIF(NOT DEFINED MMGS_INCLUDE_DIR)
        IF(NOT DEFINED MMG_LIBRARY)
            SET(MMG_LIBRARY ${MMG_ROOT}/lib/${MMG_LIBRARY_PREFIX}mmg.${MMG_LIBRARY_EXTENSION})
        ENDIF(NOT DEFINED MMG_LIBRARY)
        IF(NOT DEFINED MMG2D_LIBRARY)
            SET(MMG2D_LIBRARY ${MMG_ROOT}/lib/${MMG_LIBRARY_PREFIX}mmg2d.${MMG_LIBRARY_EXTENSION})
        ENDIF(NOT DEFINED MMG2D_LIBRARY)
        IF(NOT DEFINED MMG3D_LIBRARY)
            SET(MMG3D_LIBRARY ${MMG_ROOT}/lib/${MMG_LIBRARY_PREFIX}mmg3d.${MMG_LIBRARY_EXTENSION})
        ENDIF(NOT DEFINED MMG3D_LIBRARY)
        IF(NOT DEFINED MMGS_LIBRARY)
            SET(MMGS_LIBRARY ${MMG_ROOT}/lib/${MMG_LIBRARY_PREFIX}mmgs.${MMG_LIBRARY_EXTENSION})
        ENDIF(NOT DEFINED MMGS_LIBRARY)
    ENDIF (MMG_ROOT)

    # Add definitions
    add_definitions( -DINCLUDE_MMG )

    INCLUDE_DIRECTORIES(${MMG_INCLUDE_DIR})
    INCLUDE_DIRECTORIES(${MMG2D_INCLUDE_DIR})
    INCLUDE_DIRECTORIES(${MMG3D_INCLUDE_DIR})
    INCLUDE_DIRECTORIES(${MMGS_INCLUDE_DIR})

    target_link_libraries(KratosMeshingCore PUBLIC KratosCore tetMeshOpt ${MMG_LIBRARY} ${MMG2D_LIBRARY} ${MMG3D_LIBRARY} ${MMGS_LIBRARY})

    message("Installed MMG = " ${MMG_LIBRARY})
    install(FILES ${MMG_LIBRARY} DESTINATION libs)
    message("Installed MMG2D = " ${MMG2D_LIBRARY})
    install(FILES ${MMG2D_LIBRARY} DESTINATION libs)
    message("Installed MMG3D = " ${MMG3D_LIBRARY})
    install(FILES ${MMG3D_LIBRARY} DESTINATION libs)
    message("Installed MMGS = " ${MMGS_LIBRARY})
    install(FILES ${MMGS_LIBRARY} DESTINATION libs)
else(${INCLUDE_MMG} MATCHES ON)
    target_link_libraries(KratosMeshingCore PUBLIC KratosCore tetMeshOpt)
endif(${INCLUDE_MMG} MATCHES ON)

target_link_libraries(KratosMeshingCore PUBLIC KratosCore triangle tetMeshOpt)
set_target_properties(KratosMeshingCore PROPERTIES COMPILE_DEFINITIONS "MESHING_APPLICATION=EXPORT,API")

IF(NOT DEFINED PMMG_LIBRARY_PREFIX)
    SET(PMMG_LIBRARY_PREFIX "lib")
ENDIF(NOT DEFINED PMMG_LIBRARY_PREFIX)

if(${INCLUDE_PMMG} MATCHES ON)
    IF(NOT DEFINED INCLUDE_MMG)
        message( FATAL_ERROR "MMG not included. Please add -DINCLUDE_MMG=ON to the configure file in order to compile with ParMMG")
    ELSE(NOT DEFINED INCLUDE_MMG)
        if(${INCLUDE_MMG} MATCHES OFF)
            message( FATAL_ERROR " MMG disabled. Please set -DINCLUDE_MMG=ON in the configure file in order to compile with ParMMG")
        endif(${INCLUDE_MMG} MATCHES OFF)
    ENDIF(NOT DEFINED PMMG_LIBRARY_PREFIX)


    add_definitions( -DINCLUDE_PMMG )
    ## find Metis for ParMmg
    if(NOT DEFINED USE_METIS_5)
        message("USE_METIS_5 not defined. Setting to ON")
        SET (USE_METIS_5 ON)
    endif(NOT DEFINED USE_METIS_5)

    if(${USE_METIS_5} MATCHES ON)
        message("using metis (and not Parmetis) version 5")
        add_definitions(-DKRATOS_USE_METIS_5)
        find_package(METIS REQUIRED)
        include_directories(${METIS_INCLUDE_DIR})
        link_directories(${METIS_LIBRARIES})
        install(FILES ${METIS_LIBRARIES} DESTINATION libs)
    else(${USE_METIS_5} MATCHES ON)
        message("using PARMETIS < 4 as default. Set USE_METIS_5=ON to compile with METIS v5 instead- note that current compile option is deprecated and will be removed soon")
        find_package(PARMETIS REQUIRED)
        include_directories(${PARMETIS_INCLUDE_DIR})
        link_directories(${PARMETIS_LIBRARIES})
    endif(${USE_METIS_5} MATCHES ON)
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
        IF (MMG_ROOT)
            IF(NOT DEFINED MMG3D_INCLUDE_DIR)
                SET(MMG3D_INCLUDE_DIR ${MMG_ROOT}/include/mmg/mmg3d/)
            ENDIF(NOT DEFINED MMG3D_INCLUDE_DIR)
            IF(NOT DEFINED MMG3D_LIBRARY)
                SET(MMG3D_LIBRARY ${MMG_ROOT}/lib/${MMG_LIBRARY_PREFIX}mmg3d.${MMG_LIBRARY_EXTENSION})
            ENDIF(NOT DEFINED MMG3D_LIBRARY)
        ENDIF (MMG_ROOT)
    ENDIF (PMMG_ROOT)

    INCLUDE_DIRECTORIES(${PMMG_INCLUDE_DIR})
    INCLUDE_DIRECTORIES(${MMG3D_INCLUDE_DIR})

    message("Installed PMMG = " ${PMMG_LIBRARY})
    install(FILES ${PMMG_LIBRARY} DESTINATION libs)
    message("... with MMG3D = " ${MMG3D_LIBRARY})
    install(FILES ${MMG3D_LIBRARY} DESTINATION libs)

    target_link_libraries(KratosMeshingCore PUBLIC KratosCore KratosMPICore triangle tetMeshOpt ${METIS_LIBRARY} ${PMMG_LIBRARY} ${MMG3D_LIBRARY})
endif(${INCLUDE_PMMG} MATCHES ON)

set_target_properties(KratosMeshingCore PROPERTIES COMPILE_DEFINITIONS "MESHING_APPLICATION=EXPORT,API")

if(${USE_TETGEN_NONFREE_TPL} MATCHES ON )
    target_link_libraries(KratosMeshingCore PUBLIC ${TETGEN_LIBRARIES})
endif(${USE_TETGEN_NONFREE_TPL} MATCHES ON )

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosMeshingApplication MODULE THIN_LTO ${KRATOS_MESHING_APPLICATION_SOURCES})
target_link_libraries(KratosMeshingApplication PRIVATE KratosMeshingCore)
set_target_properties(KratosMeshingApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosMeshingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosMeshingApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    ## The MMG is not compatibel with cotire
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/mmg_process.cpp PROPERTIES COTIRE_EXCLUDED TRUE)
    cotire(KratosMeshingCore)
    cotire(KratosMeshingApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/MeshingApplication.py KratosMultiphysics/MeshingApplication/__init__.py )

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
