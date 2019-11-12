set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPfemApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../DelaunayMeshingApplication )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../SolidMechanicsApplication )

## Set a default value for triangle and tetgen in case the user does not set it
## to avoid problems with the define linkage block. By default we will compile the lib

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/DelaunayMeshingApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/SolidMechanicsApplication)

if(NOT DEFINED ${INCLUDE_TRIANGLE})
  set(INCLUDE_TRIANGLE ON)
  set(TRIANGLE_INSTALL_DIR ${CMAKE_SOURCE_DIR}/external_libraries/triangle/)
  include_directories( ${TRIANGLE_INSTALL_DIR} )
endif(NOT DEFINED ${INCLUDE_TRIANGLE})

if(NOT DEFINED ${INCLUDE_TETGEN})
  set(INCLUDE_TETGEN ON)
  set(TETGEN_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../DelaunayMeshingApplication/external_modules/tetgen/)
  include_directories( ${TETGEN_INSTALL_DIR} )
  link_directories( ${TETGEN_INSTALL_DIR} )
endif(NOT DEFINED ${INCLUDE_TETGEN})

## Generate variables with the sources
set( KRATOS_PFEM_APPLICATION_CORE
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_application_variables.cpp

	##elements
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fluid_elements/fluid_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fluid_elements/updated_lagrangian_segregated_fluid_element.cpp

)

## Generate variables with the sources
set( KRATOS_PFEM_APPLICATION_PYTHON_INTERFACE
  # Custom python
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/pfem_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
)

add_library(KratosPfemCore SHARED ${KRATOS_PFEM_APPLICATION_CORE})
target_link_libraries(KratosPfemCore PUBLIC KratosCore KratosDelaunayMeshingCore KratosSolidMechanicsCore)
set_target_properties(KratosPfemCore PROPERTIES COMPILE_DEFINITIONS "PFEM_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosPfemApplication MODULE THIN_LTO ${KRATOS_PFEM_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosPfemApplication PRIVATE KratosPfemCore)
set_target_properties(KratosPfemApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPfemApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosPfemApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosPfemCore)
    cotire(KratosPfemApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/PfemApplication.py KratosMultiphysics/PfemApplication/__init__.py )

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
install(TARGETS KratosPfemCore DESTINATION libs )
install(TARGETS KratosPfemApplication DESTINATION libs )