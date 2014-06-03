set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMultiScaleApplication ****")

######################################################################################
#
# external libraries configuration
#
######################################################################################



######################################################################################
#
# includes
#
######################################################################################

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

######################################################################################
#
# sources
#
######################################################################################

set( KRATOS_MULTISCALE_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/multiscale_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/j2_constitutive_law_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/isotropic_damage_plane_stress_2d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/scalar_damage_interface_2d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/scalar_damage_interface_3d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_interface_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rve_condition_base.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rve_corner_condition_2D4N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rve_corner_condition_3D8N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rve_periodic_condition_DX_2D2N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rve_periodic_condition_DY_2D2N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rve_periodic_condition_DX_2D3N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rve_periodic_condition_DY_2D3N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rve_periodic_condition_DX_3D2N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rve_periodic_condition_DY_3D2N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rve_periodic_condition_DZ_3D2N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/rve_boundary.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/rve_boundary_2D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/rve_boundary_3D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/multiscale_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_constitutive_laws_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_conditions_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_linear_solvers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_io_to_python.cpp
)

######################################################################################
#
# link
#
######################################################################################

add_library(KratosMultiScaleApplication SHARED ${KRATOS_MULTISCALE_APPLICATION_SOURCES})

target_link_libraries( KratosMultiScaleApplication 
                       KratosCore )

######################################################################################
#
# properties
#
######################################################################################

set_target_properties(KratosMultiScaleApplication PROPERTIES PREFIX "")

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosMultiScaleApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

######################################################################################
#
# install
#
######################################################################################

install( TARGETS 
           KratosMultiScaleApplication 
		 DESTINATION 
		   libs )

# install( FILES 
		   # ${OTHER_FILES_TO_INSTALL_HERE}
		 # DESTINATION 
		   # libs )

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MultiScaleApplication.py" DESTINATION KratosMultiphysics )
