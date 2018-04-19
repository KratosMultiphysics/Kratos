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

	message("MUMPS will not be configured")
	if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
		message("OpenBLAS will be configured")
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
	else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
		message("BLAS will be configured")
		message("${MULTISCALE_APPLICATION_BLAS_DIR}")
		set( BLAS_LIB_DIR  ${MULTISCALE_APPLICATION_BLAS_DIR})
		set( BLAS_INSTALL_FILES ${BLAS_LIB_DIR}/libblas.so)
	endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

endif(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)

######################################################################################
#
# includes
#
######################################################################################

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

if(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)
	message("MUMPS will be included")
	include_directories( ${MUMPS_INCLUDE_DIR} )
else(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)
	message("MUMPS will not be included")
endif(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)

if(${MULTISCALE_APPLICATION_USE_EIGEN} MATCHES ON)
	message("Eigen will be included")
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
	${CMAKE_CURRENT_SOURCE_DIR}/multiscale_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/conv_diff_constitutive_law_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/conv_diff_anisotropic_3d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/conv_diff_plane_stress_2d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/conv_diff_anisotropic_2d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/linear_elastic_no_tension_plane_stress_2d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/linear_elastic_thick_shell_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/j2_constitutive_law_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/damage_iso_plane_stress_2d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/damage_tc_plane_stress_2d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/damage_tc_3d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/interpolated_constitutive_law_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/interpolated_constitutive_law_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/scalar_damage_interface_2d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/scalar_damage_interface_3d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/conv_diff_interface_2d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/conv_diff_interface_3d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/plastic_damage_interface_2d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_interface_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_thick_interface_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/opt_triangle_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/eas_quad_element_v2.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/q4ristab_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/convdiff_interface_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/convdiff_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/ebst_element_2d3n.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/agq4_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_elastic_link_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/periodic_condition_lm_2D2N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rve_weak_periodic_condition_2D2N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rve_weak_rotation_constraint_condition_3D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rve_weak_tshear_constraint_condition_3D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/rve_geometry_descriptor.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/rve_predictor_calculator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/multiscale_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_constitutive_laws_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_conditions_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_linear_solvers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_io_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/simplexnoise.cpp
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
						   ${BLAS_LIBRARIES}
						   ${MUMPS_LIB} )

else(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)

	message("MUMPS will not be linked")
	target_link_libraries( KratosMultiScaleApplication
						   KratosCore
						   ${BLAS_LIBRARIES}
						   ${LAPACK_LIBRARIES}
						   )

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

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosMultiScaleApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

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

	message("MUMPS will not be installed")
	if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
		message("install OPENBLAS")
		install( FILES
				${OPENBLAS_INSTALL_FILES}
				DESTINATION
				  libs )
	else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
		message("install BLAS")
		install( FILES
				${BLAS_INSTALL_FILES}
				DESTINATION
				  libs )
	endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")


endif(${MULTISCALE_APPLICATION_USE_MUMPS} MATCHES ON)



if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MultiScaleApplication.py" DESTINATION KratosMultiphysics )


INCLUDE (CheckLibraryExists)
set(CMAKE_REQUIRED_LIBRARIES  ${LAPACK_LIBRARIES})
message("Finding required functions from LAPACK")
message("LAPACK library: " ${LAPACK_LIBRARIES})

foreach( lapack_lib ${LAPACK_LIBRARIES})
	message("checking dgels in library: " ${lapack_lib})
    CHECK_LIBRARY_EXISTS(${lapack_lib} "dgels" "" CHECK_HAVE_DGELS)
	if(CHECK_HAVE_DGELS)
		message("found dgels")
		add_definitions(-DMULTISCALE_APPLICATION_LAPACK_HAS_DGELS)
		add_definitions(-DMULTISCALE_APPLICATION_LAPACK_DGELS=dgels)
	else()
		message("checking DGELS in library: " ${lapack_lib})
        CHECK_LIBRARY_EXISTS(${lapack_lib} "DGELS" "" CHECK_HAVE_DGELS_1)
		if(CHECK_HAVE_DGELS_1)
			message("found DGELS")
			add_definitions(-DMULTISCALE_APPLICATION_LAPACK_HAS_DGELS)
			add_definitions(-DMULTISCALE_APPLICATION_LAPACK_DGELS=DGELS)
		else()
			message("checking dgels_ in library: " ${lapack_lib})
            CHECK_LIBRARY_EXISTS(${lapack_lib} "dgels_" "" CHECK_HAVE_DGELS_2)
			if(CHECK_HAVE_DGELS_2)
				message("found dgels_")
				add_definitions(-DMULTISCALE_APPLICATION_LAPACK_HAS_DGELS)
				add_definitions(-DMULTISCALE_APPLICATION_LAPACK_DGELS=dgels_)
			else()
				message("Cannot find DGELS")
			endif()
		endif()
	endif()
endforeach(lapack_lib)

