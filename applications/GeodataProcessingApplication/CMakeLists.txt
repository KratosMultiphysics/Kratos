set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosGeodataProcessingApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_GEODATA_PROCESSING_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/geodata_processing_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/geodata_processing_application_variables.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/geodata_processing_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/hello.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/extrusion_height_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/cleaning_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fill_cfd_modelpart_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/building_utilities.cpp
)

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosGeodataProcessingApplication MODULE THIN_LTO ${KRATOS_GEODATA_PROCESSING_APPLICATION_SOURCES})
target_link_libraries(KratosGeodataProcessingApplication PRIVATE KratosCore )
set_target_properties(KratosGeodataProcessingApplication PROPERTIES PREFIX "")

if(USE_COTIRE MATCHES ON)
    cotire(KratosGeodataProcessingApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosGeodataProcessingApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosGeodataProcessingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosGeodataProcessingApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )


# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# # Add to the KratosMultiphisics Python module
# install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/GeodataProcessingApplication.py" DESTINATION "KratosMultiphysics/GeodataProcessingApplication" RENAME "__init__.py")
