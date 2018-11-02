set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosDemStructuresCouplingApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/DEM_application )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive  ")

## generate variables with the sources
set(KRATOS_DEM_STRUCTURES_COUPLING_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/dem_structures_coupling_application.cpp

    ## CONDITIONS
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/surface_load_from_DEM_condition_3d.cpp
)

set(KRATOS_DEM_STRUCTURES_COUPLING_APPLICATION_PYTHON_INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/dem_structures_coupling_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
)


add_library(KratosDemStructuresCouplingCore SHARED ${KRATOS_DEM_STRUCTURES_COUPLING_APPLICATION_CORE})
target_link_libraries(KratosDemStructuresCouplingCore PUBLIC KratosDEMCore KratosStructuralMechanicsCore KratosCore)
set_target_properties(KratosDemStructuresCouplingCore PROPERTIES COMPILE_DEFINITIONS "DEM_STRUCTURES_COUPLING_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosDemStructuresCouplingApplication MODULE THIN_LTO ${KRATOS_DEM_STRUCTURES_COUPLING_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosDemStructuresCouplingApplication PUBLIC KratosDemStructuresCouplingCore)
set_target_properties(KratosDemStructuresCouplingApplication PROPERTIES PREFIX "")


if(${ACTIVATE_DEBUG_MACRO} MATCHES ON) #MSI: Flag defined for debug Macro
    add_definitions(-DDEBUG_MACRO)
endif(${ACTIVATE_DEBUG_MACRO} MATCHES ON)

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosDemStructuresCouplingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosDemStructuresCouplingApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_TESTING_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosDemStructuresCouplingApplication subdir inc_dirs = ${inc_dirs}")
# Add Cotire
if(USE_COTIRE MATCHES ON)
  cotire(KratosDemStructuresCouplingCore)
  cotire(KratosDemStructuresCouplingApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosDemStructuresCouplingCore DESTINATION libs )
install(TARGETS KratosDemStructuresCouplingApplication DESTINATION libs )

# Add to the KratosMultiphysics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/DemStructuresCouplingApplication.py" DESTINATION KratosMultiphysics )
