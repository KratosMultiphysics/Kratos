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

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosULFApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so 
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosULFApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

install(TARGETS KratosULFApplication DESTINATION libs )


if(${INSTALL_TESTING_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)


if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)
# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ULFApplication.py" DESTINATION KratosMultiphysics )
