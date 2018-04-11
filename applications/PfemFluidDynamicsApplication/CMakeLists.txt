set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPfemFluidDynamicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
#include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../SolidMechanicsApplication )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../PfemApplication )

#if(NOT (${SOLID_MECHANICS_APPLICATION} MATCHES "ON"))
#    message(SEND_ERROR "The Solid Mechanics Application must be compiled before this one. Set SOLID_MECHANICS_APPLICATION=ON in the configure file.")
#endif(NOT (${SOLID_MECHANICS_APPLICATION} MATCHES "ON"))

if(NOT (${PFEM_APPLICATION} MATCHES "ON"))
    message(SEND_ERROR "The Pfem Base Application must be compiled before this one. Set PFEM_APPLICATION=ON in the configure file.")
endif(NOT (${PFEM_APPLICATION} MATCHES "ON"))

if(NOT DEFINED ${INCLUDE_TRIANGLE})
  set(INCLUDE_TRIANGLE ON)
  set(TRIANGLE_INSTALL_DIR ${CMAKE_SOURCE_DIR}/external_libraries/triangle/)
  include_directories( ${TRIANGLE_INSTALL_DIR} )
endif(NOT DEFINED ${INCLUDE_TRIANGLE})

if(NOT DEFINED ${INCLUDE_TETGEN})
  set(INCLUDE_TETGEN ON)
  set(TETGEN_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../PfemApplication/external_modules/tetgen/)
  include_directories( ${TETGEN_INSTALL_DIR} )
  link_directories( ${TETGEN_INSTALL_DIR} )
endif(NOT DEFINED ${INCLUDE_TETGEN})


## generate variables with the sources
set( KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_fluid_dynamics_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_fluid_dynamics_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/pfem_fluid_dynamics_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_modelers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_bounding_to_python.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/two_step_updated_lagrangian_V_P_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/two_step_updated_lagrangian_V_P_solid_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_V_solid_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/two_step_updated_lagrangian_V_P_fluid_element.cpp

)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosPfemFluidDynamicsApplication SHARED ${KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_SOURCES})
#####***target_link_libraries(KratosPfemFluidDynamicsApplication KratosCore)
#target_link_libraries(KratosPfemFluidDynamicsApplication KratosCore SolidMechanicsApplication KratosPfemApplication)
target_link_libraries(KratosPfemFluidDynamicsApplication KratosCore KratosPfemApplication)
set_target_properties(KratosPfemFluidDynamicsApplication PROPERTIES PREFIX "")
set_target_properties(KratosPfemFluidDynamicsApplication PROPERTIES COMPILE_DEFINITIONS "PFEM_FLUID_DYNAMICS_APPLICATION=EXPORT,API")
install(TARGETS KratosPfemFluidDynamicsApplication DESTINATION libs )


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

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PfemFluidDynamicsApplication.py" DESTINATION KratosMultiphysics )
