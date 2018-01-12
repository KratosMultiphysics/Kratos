set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosStabilizedCFDApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_STABILIZED_CFD_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/stabilized_cfd_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/stabilized_cfd_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/stabilized_cfd_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/turbulence_statistics_container.cpp
)

# define library Kratos which defines the basic python interface
add_library(KratosStabilizedCFDApplication SHARED ${KRATOS_STABILIZED_CFD_APPLICATION_SOURCES})
target_link_libraries(KratosStabilizedCFDApplication KratosCore )
set_target_properties(KratosStabilizedCFDApplication PROPERTIES PREFIX "")
install(TARGETS KratosStabilizedCFDApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosStabilizedCFDApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosStabilizedCFDApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/StabilizedCFDApplication.py" DESTINATION KratosMultiphysics )
