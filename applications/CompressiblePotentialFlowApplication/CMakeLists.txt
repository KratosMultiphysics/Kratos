set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosCompressiblePotentialFlowApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/compressible_potential_flow_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/compressible_potential_flow_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/compressible_potential_flow_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
)

## define compressible potential application test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
 	file(GLOB_RECURSE KRATOS_COMPRESSIBLE_POTENTIAL_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosCompressiblePotentialFlowApplication MODULE THIN_LTO ${KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION_SOURCES} ${KRATOS_COMPRESSIBLE_POTENTIAL_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosCompressiblePotentialFlowApplication PRIVATE KratosCore )
set_target_properties(KratosCompressiblePotentialFlowApplication PROPERTIES PREFIX "")
install(TARGETS KratosCompressiblePotentialFlowApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosCompressiblePotentialFlowApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosCompressiblePotentialFlowApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/CompressiblePotentialFlowApplication.py" DESTINATION KratosMultiphysics )
