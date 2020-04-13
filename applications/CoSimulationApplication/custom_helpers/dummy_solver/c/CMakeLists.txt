project( DummySolverC C)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../../../co_sim_io/c )
link_directories(${CMAKE_SOURCE_DIR}/libs)

add_executable (dummy_solver_c dummy_solver.c)
target_link_libraries(dummy_solver_c co_sim_c_io)

# using PARENT_DIRECTORY to check whether installed in Kratos or standalone
get_directory_property( IN_KRATOS PARENT_DIRECTORY )
if(IN_KRATOS)
    # if compiled from inside of Kratos, then install to the "libs" folder
    install(TARGETS dummy_solver_c DESTINATION libs )
endif()