set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosCoSimulationApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_CO_SIMULATION_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/co_simulation_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/co_simulation_python_application.cpp
)

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosCoSimulationApplication MODULE THIN_LTO ${KRATOS_CO_SIMULATION_APPLICATION_SOURCES})
target_link_libraries(KratosCoSimulationApplication KratosCore )
set_target_properties(KratosCoSimulationApplication PROPERTIES PREFIX "")
install(TARGETS KratosCoSimulationApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosCoSimulationApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosCoSimulationApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(DIRECTORY DESTINATION "KratosMultiphysics/CoSimulationApplication")
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/CoSimulationApplication.py" DESTINATION "KratosMultiphysics/CoSimulationApplication" RENAME "__init__.py")

