set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFSIApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_FSI_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/fsi_application.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/fsi_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_mappers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_convergence_accelerators_to_python.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/AdvancedNMPointsMapper.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosFSIApplication MODULE THIN_LTO ${KRATOS_FSI_APPLICATION_SOURCES})

target_link_libraries(KratosFSIApplication PRIVATE KratosCore )
set_target_properties(KratosFSIApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosFSIApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosFSIApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add Cotire
if(USE_COTIRE MATCHES ON)
  cotire(KratosFSIApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosFSIApplication DESTINATION libs )

if(${INSTALL_TESTING_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)
# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/FSIApplication.py" DESTINATION KratosMultiphysics )
