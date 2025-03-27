set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosDEMFEMVolumeCouplingApplication ****")

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/DEMApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/MappingApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/CoSimulationApplication)

# # This should replace the older kratos_add_dependency(...) calls
# kratos_add_dependency_to_application(
#     DEMFEM_VOLUME_COUPLING_APPLICATION       # <-- Your application "internal" name
#     DEM_APPLICATION
#     STRUCTURAL_MECHANICS_APPLICATION
#     MAPPING_APPLICATION
#     COSIMULATION_APPLICATION
# )


################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${KRATOS_SOURCE_DIR}/applications/DEMApplication )
include_directories( ${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication )

## generate variables with the sources
set(KRATOS_DEMFEM_VOLUME_COUPLING_APPLICATION_CORE

    ${CMAKE_CURRENT_SOURCE_DIR}/DEMFEM_volume_coupling_application.h
    ${CMAKE_CURRENT_SOURCE_DIR}/DEMFEM_volume_coupling_application.cpp
    
    ## ELEMENTS FEM
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/volume_coupling_element.h
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/volume_coupling_element.cpp

    ## ELEMENTS DEM
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/volume_coupling_particle.h
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/volume_coupling_particle.cpp
       
)

set(KRATOS_DEMFEM_VOLUME_COUPLING_APPLICATION_PYTHON_INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/DEMFEM_volume_coupling_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
)


add_library(KratosDEMFEMVolumeCouplingCore SHARED ${KRATOS_DEMFEM_VOLUME_COUPLING_APPLICATION_CORE})
target_link_libraries(KratosDEMFEMVolumeCouplingCore PUBLIC KratosDEMCore KratosStructuralMechanicsCore KratosCore)
set_target_properties(KratosDEMFEMVolumeCouplingCore PROPERTIES COMPILE_DEFINITIONS "DEMFEM_VOLUME_COUPLING_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosDEMFEMVolumeCouplingApplication MODULE THIN_LTO ${KRATOS_DEMFEM_VOLUME_COUPLING_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosDEMFEMVolumeCouplingApplication PUBLIC KratosDEMFEMVolumeCouplingCore)
set_target_properties(KratosDEMFEMVolumeCouplingApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosDEMFEMVolumeCouplingCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosDEMFEMVolumeCouplingApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

if(${ACTIVATE_DEBUG_MACRO} MATCHES ON) #MSI: Flag defined for debug Macro
    add_definitions(-DDEBUG_MACRO)
endif(${ACTIVATE_DEBUG_MACRO} MATCHES ON)

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosDEMFEMVolumeCouplingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosDEMFEMVolumeCouplingApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/DEMFEM_volume_coupling_application.py KratosMultiphysics/DEMFEMVolumeCouplingApplication/__init__.py)

# Install python files
#get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
#kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
# if(${INSTALL_TESTING_FILES} MATCHES ON )
#     get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
#     install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
#         PATTERN "*.git" EXCLUDE
#         PATTERN "*.c" EXCLUDE
#         PATTERN "*.h" EXCLUDE
#         PATTERN "*.cpp" EXCLUDE
#         PATTERN "*.hpp" EXCLUDE
#   )
# endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosDEMFEMVolumeCouplingCore DESTINATION libs )
install(TARGETS KratosDEMFEMVolumeCouplingApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosDEMFEMVolumeCouplingCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosDEMFEMVolumeCouplingApplication" PARENT_SCOPE)