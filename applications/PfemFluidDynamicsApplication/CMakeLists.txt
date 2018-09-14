set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPfemFluidDynamicsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../DelaunayMeshingApplication )

if(NOT (${DELAUNAY_MESHING_APPLICATION} MATCHES "ON"))
    message(SEND_ERROR "The Delaunay Meshing Application must be compiled before this one. Set DELAUNAY_MESHING_APPLICATION=ON in the configure file.")
endif(NOT (${DELAUNAY_MESHING_APPLICATION} MATCHES "ON"))

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


## generate variables with the sources
set( KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_CORE
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_fluid_dynamics_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_fluid_dynamics_application_variables.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/two_step_updated_lagrangian_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/two_step_updated_lagrangian_V_P_implicit_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/two_step_updated_lagrangian_V_P_implicit_nodally_integrated_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/two_step_updated_lagrangian_V_P_implicit_solid_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/two_step_updated_lagrangian_V_P_implicit_nodally_integrated_solid_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/two_step_updated_lagrangian_V_P_implicit_nodally_integrated_fluid_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_V_implicit_solid_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/two_step_updated_lagrangian_V_P_implicit_fluid_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/two_step_updated_lagrangian_V_P_explicit_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/two_step_updated_lagrangian_V_P_explicit_solid_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_V_explicit_solid_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/two_step_updated_lagrangian_V_P_explicit_fluid_element.cpp
)

set( KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_PYTHON_INTERFACE
        # Custom python
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/pfem_fluid_dynamics_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_modelers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_bounding_to_python.cpp
)

add_library(KratosPfemFluidDynamicsCore SHARED ${KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_CORE})
target_link_libraries(KratosPfemFluidDynamicsCore PUBLIC KratosCore KratosDelaunayMeshingCore)
set_target_properties(KratosPfemFluidDynamicsCore PROPERTIES COMPILE_DEFINITIONS "PFEM_FLUID_DYNAMICS_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosPfemFluidDynamicsApplication MODULE THIN_LTO ${KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosPfemFluidDynamicsApplication PRIVATE KratosPfemFluidDynamicsCore) 
set_target_properties(KratosPfemFluidDynamicsApplication PROPERTIES PREFIX "")

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPfemFluidDynamicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosPfemFluidDynamicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

if(USE_COTIRE MATCHES ON)
    cotire(KratosPfemFluidDynamicsCore)
    cotire(KratosPfemFluidDynamicsApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosPfemFluidDynamicsCore DESTINATION libs )
install(TARGETS KratosPfemFluidDynamicsApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PfemFluidDynamicsApplication.py" DESTINATION KratosMultiphysics )
