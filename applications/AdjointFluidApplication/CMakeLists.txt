set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosAdjointFluidApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_ADJOINT_FLUID_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/adjoint_fluid_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/adjoint_fluid_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_schemes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_response_functions_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/adjoint_fluid_python_application.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/vms_adjoint_element.cpp
    )

###############################################################
## define library which defines the basic python interface
pybind11_add_module(KratosAdjointFluidApplication MODULE THIN_LTO ${KRATOS_ADJOINT_FLUID_APPLICATION_SOURCES})
target_link_libraries(KratosAdjointFluidApplication PRIVATE KratosCore ${HDF5_CXX_LIBRARIES})
set_target_properties(KratosAdjointFluidApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosAdjointFluidApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so 
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosAdjointFluidApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

install(TARGETS KratosAdjointFluidApplication DESTINATION libs)

if(${INSTALL_PYTHON_FILES} MATCHES ON)
    get_filename_component(CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/AdjointFluidApplication.py" DESTINATION KratosMultiphysics )
