set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFluidTransportApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

#dependecies on other applications
include_directories( ${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication )

# generate variables with the sources
set( KRATOS_FLUID_TRANSPORT_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/fluid_transport_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fluid_transport_application_variables.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/fluid_transport_python_application.cpp
  
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/steady_convection_diffusion_FIC_element.cpp
)

# define library Kratos which defines the basic python interface
add_library(KratosFluidTransportApplication SHARED ${KRATOS_FLUID_TRANSPORT_APPLICATION_SOURCES})
target_link_libraries(KratosFluidTransportApplication KratosCore)
set_target_properties(KratosFluidTransportApplication PROPERTIES PREFIX "")
install(TARGETS KratosFluidTransportApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosFluidTransportApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosFluidTransportApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/FluidTransportApplication.py" DESTINATION KratosMultiphysics )
