set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMultiScaleApplication ****")

######################################################################################
#
# external libraries configuration
#
######################################################################################

if(${MULTISCALE_APPLICATION_USE_EIGEN} MATCHES ON)

	add_definitions(-DMULTISCALE_APPLICATION_USE_EIGEN)

	message("Configuring Eigen")
	message(${MULTISCALE_APPLICATION_EIGEN_DIR})
	
endif(${MULTISCALE_APPLICATION_USE_EIGEN} MATCHES ON)

if(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)

	add_definitions(-DMULTISCALE_APPLICATION_USE_MUMPS)

	message("Configuring OpenBLAS")
	message(${MULTISCALE_APPLICATION_OPENBLAS_DIR})
	set( OPENBLAS_LIB_DIR  ${MULTISCALE_APPLICATION_OPENBLAS_DIR} )
	set( OPENBLAS_LIB ${OPENBLAS_LIB_DIR}/libopenblas.lib )
	set( OPENBLAS_INSTALL_FILES
		${OPENBLAS_LIB_DIR}/libopenblas.dll
		${OPENBLAS_LIB_DIR}/libgcc_s_sjlj_64-1.dll
		${OPENBLAS_LIB_DIR}/libgfortran_64-3.dll
		${OPENBLAS_LIB_DIR}/libquadmath_64-0.dll
		${OPENBLAS_LIB_DIR}/pthreadGC2_64.dll
	)

	message("Configuring MUMPS")
	message(${MULTISCALE_APPLICATION_MUMPS_INCLUDE_DIR})
	message(${MULTISCALE_APPLICATION_MUMPS_LIB_DIR})
	set( MUMPS_INCLUDE_DIR  ${MULTISCALE_APPLICATION_MUMPS_INCLUDE_DIR} )
	set( MUMPS_LIB_DIR      ${MULTISCALE_APPLICATION_MUMPS_LIB_DIR} )
	set( MUMPS_LIB ${MUMPS_LIB_DIR}/libmumps.lib )
	set( MUMPS_INSTALL_FILES ${MUMPS_LIB_DIR}/libmumps.dll )

else(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)

	message("MUMPS and OpenBLAS will not be configured")

endif(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)

######################################################################################
#
# includes
#
######################################################################################

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

if(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)
	message("include MUMPS")
	include_directories( ${MUMPS_INCLUDE_DIR} )
else(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)
	message("MUMPS will not be included")
endif(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)

if(${MULTISCALE_APPLICATION_USE_EIGEN} MATCHES ON)
	message("include Eigen")
	include_directories( ${MULTISCALE_APPLICATION_EIGEN_DIR} )
else(${MULTISCALE_APPLICATION_USE_EIGEN} MATCHES ON)
	message("Eigen will not be included")
endif(${MULTISCALE_APPLICATION_USE_EIGEN} MATCHES ON)

######################################################################################
#
# sources
#
######################################################################################

set( KRATOS_MULTISCALE_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/multiscale_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/j2_constitutive_law_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/damage_iso_plane_stress_2d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/damage_tc_plane_stress_2d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/damage_tc_3d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/scalar_damage_interface_2d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/scalar_damage_interface_3d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/plastic_damage_interface_2d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_interface_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/opt_triangle_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/eas_quad_element_v2.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/q4ristab_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/ebst_element_2d3n.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/agq4_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/periodic_condition_lm_2D2N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rve_weak_periodic_condition_2D2N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/rve_geometry_descriptor.cpp
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

if(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)

	message("link also MUMPS and OpenBLAS")
	target_link_libraries( KratosMultiScaleApplication 
						   KratosCore
						   ${OPENBLAS_LIB}
						   ${MUMPS_LIB} )

else(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)

	message("MUMPS and OpenBLAS will not be linked")
	target_link_libraries( KratosMultiScaleApplication 
						   KratosCore )

endif(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)



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

if(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)

	message("install also MUMPS and OpenBLAS")
	install( FILES 
			  ${OPENBLAS_INSTALL_FILES}
			  ${MUMPS_INSTALL_FILES}
			DESTINATION 
			  libs )

else(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)

	message("MUMPS and OpenBLAS will not be installed")

endif(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)



if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MultiScaleApplication.py" DESTINATION KratosMultiphysics )
