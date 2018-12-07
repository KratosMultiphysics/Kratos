set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMeshingApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/external_libraries/triangle )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication ) 

if(${USE_TETGEN_NONFREE_TPL} MATCHES ON )
    add_definitions( -DTETGEN143 )
    add_definitions( -DTETLIBRARY )
    add_definitions( -DUSE_TETGEN_NONFREE_TPL )
    include_directories(${TETGEN_INCLUDES})
    #include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetgen/tetgen1.4.3 )
    #add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetgen/tetgen1.4.3 )
endif(${USE_TETGEN_NONFREE_TPL} MATCHES ON )

##add subdirectories for compilation of dependent libraries
#add_subdirectory(${CMAKE_SOURCE_DIR}/external_libraries/triangle)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/triangle)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt)

## generate variables with the sources
if(${INCLUDE_MMG} MATCHES ON)
set( KRATOS_MESHING_APPLICATION_CORE  
    ${CMAKE_CURRENT_SOURCE_DIR}/meshing_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/external_includes/mesh_suites.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/local_refine_geometry_mesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/uniform_refinement_utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_includes/meshing_flags.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/metrics_levelset_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/internal_variables_interpolation_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/metric_fast_init_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/metrics_hessian_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/metrics_error_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/nodal_values_interpolation_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/integration_values_extrapolation_to_nodes_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/multiscale_refining_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/mmg_process.cpp
)
else(${INCLUDE_MMG} MATCHES ON)
set( KRATOS_MESHING_APPLICATION_CORE  
    ${CMAKE_CURRENT_SOURCE_DIR}/meshing_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/external_includes/mesh_suites.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/local_refine_geometry_mesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/uniform_refinement_utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_includes/meshing_flags.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/metrics_levelset_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/internal_variables_interpolation_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/metric_fast_init_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/metrics_hessian_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/metrics_error_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/nodal_values_interpolation_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/integration_values_extrapolation_to_nodes_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/multiscale_refining_process.cpp
)
endif(${INCLUDE_MMG} MATCHES ON)

set (KRATOS_MESHING_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/kratos_meshing_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_meshers_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
)

if(${KRATOS_BUILD_TESTING} MATCHES ON) 
    file(GLOB_RECURSE KRATOS_MESHING_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp) 
endif(${KRATOS_BUILD_TESTING} MATCHES ON) 

add_library(KratosMeshingCore SHARED ${KRATOS_MESHING_APPLICATION_CORE} ${KRATOS_MESHING_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosMeshingCore PUBLIC KratosCore triangle tetMeshOpt)
set_target_properties(KratosMeshingCore PROPERTIES COMPILE_DEFINITIONS "MESHING_APPLICATION=EXPORT,API")

if(${INCLUDE_MMG} MATCHES ON)
    add_definitions( -DINCLUDE_MMG )

    INCLUDE_DIRECTORIES(${MMG_INCLUDE_DIR})
    INCLUDE_DIRECTORIES(${MMG2D_INCLUDE_DIR})
    INCLUDE_DIRECTORIES(${MMG3D_INCLUDE_DIR})
    INCLUDE_DIRECTORIES(${MMGS_INCLUDE_DIR})

    target_link_libraries(KratosMeshingCore KratosCore triangle tetMeshOpt ${MMG_LIBRARY} ${MMG2D_LIBRARY} ${MMG3D_LIBRARY} ${MMGS_LIBRARY})
endif(${INCLUDE_MMG} MATCHES ON)

if(${USE_TETGEN_NONFREE_TPL} MATCHES ON )
    target_link_libraries(KratosMeshingCore ${TETGEN_LIBRARIES})
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

if(USE_COTIRE MATCHES ON)
    ## The MMG is not compatibel with cotire
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/mmg_process.cpp PROPERTIES COTIRE_EXCLUDED TRUE)
    cotire(KratosMeshingApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosMeshingCore DESTINATION libs )
install(TARGETS KratosMeshingApplication DESTINATION libs )

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)


# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MeshingApplication.py" DESTINATION KratosMultiphysics )
