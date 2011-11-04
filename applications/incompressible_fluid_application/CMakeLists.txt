set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosIncompressibleFluidApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_INCOMPRESSIBLE_FLUID_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/incompressible_fluid_application.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/kratos_incompressible_fluid_python_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_io_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_edgebased_levelset_solver_to_python.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/proj_dirichlet_cond.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/fluid3d_neumann.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/no_slip_condition_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/monolithic2d_neumann.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fluid_2d_split.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fluid_3d.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fluid_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fluid_2dcoupled.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fluid_3dcoupled.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/NDfluid_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/NDfluid_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/NDfluid_2d_CrankNicolson.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fluid_2dGLS.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fluid_2dGLS_expl.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fluid_3dGLS_expl.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fluid_2dGLS_expl_comp.cpp   
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fluid_3dGLS_expl_comp.cpp   
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/asgs_2d.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/asgs_3d.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/asgs_pr_dc_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/asgs_compressible_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/asgs_compressible_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/asgs_comp_pr_dc_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/asgs_comp_pr_dc_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/nonewtonian_asgs_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/nonewtonian_asgs_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/bingham_nonnewtonian_asgs_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/bingham_nonnewtonian_asgs_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/explicit_asgs_compressible_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/explicit_asgs_comp_pr_dc_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/explicit_asgs_compressible_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/explicit_asgs_comp_pr_dc_3d.cpp     
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosIncompressibleFluidApplication SHARED ${KRATOS_INCOMPRESSIBLE_FLUID_APPLICATION_SOURCES})
target_link_libraries(KratosIncompressibleFluidApplication KratosCore )
set_target_properties(KratosIncompressibleFluidApplication PROPERTIES PREFIX "")

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosIncompressibleFluidApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

install(TARGETS KratosIncompressibleFluidApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")
