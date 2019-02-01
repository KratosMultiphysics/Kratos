set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosSolversApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories(${CMAKE_SOURCE_DIR}/kratos)

# generate variables with the sources
set( KRATOS_SOLVERS_APPLICATION_CORE
	${CMAKE_CURRENT_SOURCE_DIR}/solvers_application.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/solvers_application_variables.cpp

	# Time integration methods
	${CMAKE_CURRENT_SOURCE_DIR}/custom_solvers/solution_local_flags.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_solvers/time_integration_methods/static_step_rotation_method.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_solvers/time_integration_methods/newmark_step_rotation_method.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_solvers/time_integration_methods/simo_step_rotation_method.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_solvers/time_integration_methods/emc_step_rotation_method.cpp
)

## generate variables with the testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
  file(GLOB_RECURSE KRATOS_SOLVERS_APPLICATION_TESTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## Generate variables with the sources
set( KRATOS_SOLVERS_APPLICATION_PYTHON_INTERFACE
        # Custom python
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/solvers_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
)

add_library(KratosSolversCore SHARED ${KRATOS_SOLVERS_APPLICATION_CORE})
target_link_libraries(KratosSolversCore PUBLIC KratosCore)
set_target_properties(KratosSolversCore PROPERTIES COMPILE_DEFINITIONS "SOLVERS_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosSolversApplication MODULE THIN_LTO ${KRATOS_SOLVERS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosSolversApplication PRIVATE KratosSolversCore)
set_target_properties(KratosSolversApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosSolversApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set_target_properties(KratosSolversApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install testing files
if(${INSTALL_TESTING_FILES} MATCHES ON)
  get_filename_component(CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component(CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

if(USE_COTIRE MATCHES ON)
  cotire(KratosSolversCore)
  cotire(KratosSolversApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosSolversCore DESTINATION libs )
install(TARGETS KratosSolversApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/SolversApplication.py" DESTINATION "KratosMultiphysics/SolversApplication" RENAME "__init__.py")
