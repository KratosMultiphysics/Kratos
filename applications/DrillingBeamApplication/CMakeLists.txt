set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosDrillingBeamApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos  )

## generate variables with the sources
set( KRATOS_DRILLING_BEAM_APPLICATION_CORE
	${CMAKE_CURRENT_SOURCE_DIR}/drilling_beam_application.cpp
)

set(KRATOS_DRILLING_BEAM_APPLICATION_PYTHON_INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/drilling_beam_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
)

add_library(KratosDrillingBeamApplicationCore SHARED ${KRATOS_DRILLING_BEAM_APPLICATION_CORE})
target_link_libraries(KratosDrillingBeamApplicationCore PUBLIC KratosCore)
set_target_properties(KratosDrillingBeamApplicationCore PROPERTIES COMPILE_DEFINITIONS "DRILLING_BEAM_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosDrillingBeamApplication MODULE THIN_LTO ${KRATOS_DRILLING_BEAM_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosDrillingBeamApplication PUBLIC KratosDrillingBeamApplicationCore KratosCore )
set_target_properties(KratosDrillingBeamApplication PROPERTIES PREFIX "")


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosDrillingBeamApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

install(TARGETS KratosDrillingBeamApplicationCore DESTINATION libs )
install(TARGETS KratosDrillingBeamApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/DrillingBeamApplication.py" DESTINATION "KratosMultiphysics/DrillingBeamApplication" RENAME "__init__.py" )
