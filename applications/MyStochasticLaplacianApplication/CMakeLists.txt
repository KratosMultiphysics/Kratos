set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMyStochasticLaplacianApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_MY_STOCHASTIC_LAPLACIAN_APPLICATION_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/my_stochastic_laplacian_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/my_stochastic_laplacian_application_variables.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/my_stochastic_laplacian_python_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/my_stochastic_laplacian_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/my_stochastic_laplacian_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_source_condition.cpp
)

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosMyStochasticLaplacianApplication MODULE THIN_LTO ${KRATOS_MY_STOCHASTIC_LAPLACIAN_APPLICATION_SOURCES})
target_link_libraries(KratosMyStochasticLaplacianApplication PRIVATE KratosCore )
set_target_properties(KratosMyStochasticLaplacianApplication PROPERTIES PREFIX "")

if(USE_COTIRE MATCHES ON)
    cotire(KratosMyStochasticLaplacianApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosMyStochasticLaplacianApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosMyStochasticLaplacianApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosMyStochasticLaplacianApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MyStochasticLaplacianApplication.py" DESTINATION KratosMultiphysics )
