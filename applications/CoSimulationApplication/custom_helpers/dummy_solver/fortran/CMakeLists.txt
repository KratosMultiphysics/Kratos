project( DummySolverFortran Fortran )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../../../co_sim_io/fortran )
link_directories(${CMAKE_SOURCE_DIR}/libs)

add_executable (dummy_solver_fortran ${CMAKE_CURRENT_SOURCE_DIR}/dummy_solver.f90)
target_link_libraries(dummy_solver_fortran co_sim_c_io)

# using PARENT_DIRECTORY to check whether installed in Kratos or standalone
get_directory_property( IN_KRATOS PARENT_DIRECTORY )
if(IN_KRATOS)
    # if compiled from inside of Kratos, then install to the "libs" folder
    install(TARGETS dummy_solver_fortran DESTINATION libs )
endif()