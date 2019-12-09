set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosULFApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_ULF_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/ULF_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ULF_application_variables.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/kratos_ulf_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_io_to_python.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_fluid.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_fluid3D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_fluid_inc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_fluid3D_inc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/ulf_frac2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/ulf_frac3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/ulf_axisym.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fluid_2dGLS_expl.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/Point_Neumann2D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/Point_Neumann3D.cpp
 	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/Point_Neumann_Axisym.cpp
 	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/surface_tension.cpp

    #${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/ulf_frac2d_swimming.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/ulf_frac3d_swimming.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosULFApplication MODULE ${KRATOS_ULF_APPLICATION_SOURCES} )
target_link_libraries(KratosULFApplication PRIVATE KratosCore )
set_target_properties(KratosULFApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosULFApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosULFApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosULFApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ULFApplication.py KratosMultiphysics/ULFApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
        PATTERN "*.git" EXCLUDE
        PATTERN "*.c" EXCLUDE
        PATTERN "*.h" EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.hpp" EXCLUDE
  )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

## Install targets
install(TARGETS KratosULFApplication DESTINATION libs )
