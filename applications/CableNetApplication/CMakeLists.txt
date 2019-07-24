set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosCableNetApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_CABLE_NET_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/cable_net_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cable_net_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/cable_net_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_test_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/custom_test_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_test_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/custom_test_condition.cpp



  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/sliding_cable_element_3D.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/ring_element_3D.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/weak_coupling_slide.cpp
)

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosCableNetApplication MODULE THIN_LTO ${KRATOS_CABLE_NET_APPLICATION_SOURCES})
target_link_libraries(KratosCableNetApplication PRIVATE KratosCore )
set_target_properties(KratosCableNetApplication PROPERTIES PREFIX "")

if(USE_COTIRE MATCHES ON)
    cotire(KratosCableNetApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosCableNetApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosCableNetApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosCableNetApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/CableNetApplication.py" DESTINATION "KratosMultiphysics/CableNetApplication" RENAME "__init__.py")
