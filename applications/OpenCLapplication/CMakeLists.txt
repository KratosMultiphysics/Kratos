set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosOpenCLApplication ****")

#find mkl solver package -- this is our own find package!!
find_package(OPENCL REQUIRED)
include_directories(${OPENCL_INCLUDE_DIRS})
link_directories(${OPENCL_LIBRARIES})

add_definitions( -DVIENNACL_EXPERIMENTAL_DOUBLE_PRECISION_WITH_STREAM_SDK )

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries )

#include test solvers
add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/opencl_linear_solver )
add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/vcl_solver )

# message("OPENCL_INCLUDE_DIRS = ${OPENCL_INCLUDE_DIRS}")
# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("applications subdir inc_dirs = ${inc_dirs}")

## generate variables with the sources
set( KRATOS_OPENCL_APPLICATION_SOURCES    
	${CMAKE_CURRENT_SOURCE_DIR}/OpenCLapplication.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/opencl_python_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_linear_solvers_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_utilities_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/pugixml/src/pugixml.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_spatial_containers_to_python.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosOpenCLApplication SHARED ${KRATOS_OPENCL_APPLICATION_SOURCES})
target_link_libraries(KratosOpenCLApplication KratosCore ${OPENCL_LIBRARIES} rt )
set_target_properties(KratosOpenCLApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosOpenCLApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

install(TARGETS KratosOpenCLApplication DESTINATION libs )

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/OpenCLApplication.py" DESTINATION KratosMultiphysics )
