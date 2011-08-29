set(CMAKE_INCLUDE_CURRENT_DIR ON)


# ## add the include directories needed for the compilation
# message("current dir = " ${CMAKE_CURRENT_SOURCE_DIR})
# include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/sources )
# include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/container )
# include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/python )


## generate variables with the sources
set( KRATOS_CORE_SOURCES     
 	${CMAKE_CURRENT_SOURCE_DIR}/sources/global_variables.cpp;
 	${CMAKE_CURRENT_SOURCE_DIR}/sources/variables.cpp;
 	${CMAKE_CURRENT_SOURCE_DIR}/sources/datafile_io.cpp; 
     	${CMAKE_CURRENT_SOURCE_DIR}/sources/integration_rules.cpp;
 	${CMAKE_CURRENT_SOURCE_DIR}/containers/variables_list_data_value_container.cpp;
 	${CMAKE_CURRENT_SOURCE_DIR}/sources/mmio.c;
 	${CMAKE_CURRENT_SOURCE_DIR}/sources/timer.cpp;
 	${CMAKE_CURRENT_SOURCE_DIR}/sources/serializer.cpp;
     )

set( KRATOS_PYTHON_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_banded_matrix_to_python.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_model_part_to_python.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_compressed_matrix_to_python.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_constant_matrix_to_python.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_points_to_python.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_containers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_constitutive_law_to_python.cpp       
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_processes_to_python.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_coordinate_matrix_to_python.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_quadratures_to_python.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_geometries_to_python.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_scalar_matrix_to_python.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_hermitian_matrix_to_python.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_sparse_matrix_to_python.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_identity_matrix_to_python.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_integration_points_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_symmetric_matrix_to_python.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_io_to_python.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_triangular_matrix_to_python.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_kernel_to_python.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_vector_to_python.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_kratos_application_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_zero_matrix_to_python.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_matrix_to_python.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_mesh_to_python.cpp			
	${CMAKE_CURRENT_SOURCE_DIR}/python/vectorial_spaces_python.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_properties_to_python.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_node_to_python.cpp			
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_linear_solvers_to_python.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_strategies_to_python.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_utilities_to_python.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_process_info_to_python.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_parallel_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_parallel_linear_solvers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_matrix_market_interface_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/python/add_serializer_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/kratos_python.cpp
     )

###############################################################
## define library KratosCore to be included in all of the others
add_library(KratosCore SHARED ${KRATOS_CORE_SOURCES})
target_link_libraries(KratosCore ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} gidpost )
install(TARGETS KratosCore DESTINATION libs )

###############################################################
## define library Kratos which defines the basic python interface
add_library(Kratos SHARED ${KRATOS_PYTHON_SOURCES})
target_link_libraries(Kratos KratosCore )
set_target_properties(Kratos PROPERTIES PREFIX "")
install(TARGETS Kratos DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("subdir inc_dirs = ${inc_dirs}")