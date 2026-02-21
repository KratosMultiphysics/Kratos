set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFemToDemApplication ****")

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/PfemFluidDynamicsApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/DEMApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/DemStructuresCouplingApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/ConstitutiveLawsApplication)

################### PYBIND11
include(pybind11Tools)

include_directories(${KRATOS_SOURCE_DIR}/kratos )
include_directories(${KRATOS_SOURCE_DIR}/applications/DEMApplication)
include_directories(${KRATOS_SOURCE_DIR}/applications/PfemFluidDynamicsApplication)
include_directories(${KRATOS_SOURCE_DIR}/applications/DelaunayMeshingApplication)
include_directories(${KRATOS_SOURCE_DIR}/applications/DemStructuresCouplingApplication)
include_directories(${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication)
include_directories(${KRATOS_SOURCE_DIR}/applications/ConstitutiveLawsApplication)

file(GLOB_RECURSE KRATOS_FEM_TO_DEM_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/fem_to_dem_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fem_to_dem_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

file(GLOB_RECURSE KRATOS_FEM_TO_DEM_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosFemToDemApplicationCore SHARED ${KRATOS_FEM_TO_DEM_APPLICATION_CORE})
target_link_libraries(KratosFemToDemApplicationCore PUBLIC KratosDEMCore KratosCore KratosPfemFluidDynamicsCore KratosConstitutiveLawsCore KratosStructuralMechanicsCore)
set_target_properties(KratosFemToDemApplicationCore PROPERTIES COMPILE_DEFINITIONS "FEM_TO_DEM_APPLICATION=EXPORT,API")

## define library Kratos which defines the basic python interface
pybind11_add_module(KratosFemToDemApplication MODULE THIN_LTO ${KRATOS_FEM_TO_DEM_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosFemToDemApplication PRIVATE KratosFemToDemApplicationCore)
set_target_properties(KratosFemToDemApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosFemToDemApplicationCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosFemToDemApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosFemToDemApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosFemToDemApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/FemToDemApplication.py KratosMultiphysics/FemToDemApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    kratos_python_install_tests(
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tests"
        DESTINATION "applications/${CURRENT_DIR_NAME}" 
        INSTALL_PYTHON_USING_LINKS "${INSTALL_PYTHON_USING_LINKS}"
        EXCLUSION_PATTERNS "\.git$" "\.c$" "\.h$" "\.cpp$" "\.hpp$"
    )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosFemToDemApplicationCore DESTINATION libs )
install(TARGETS KratosFemToDemApplication DESTINATION libs )


# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosFemToDemApplicationCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosFemToDemApplication" PARENT_SCOPE)