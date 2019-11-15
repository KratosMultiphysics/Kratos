project( DummySolverFortran Fortran )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../../../co_simulation_io/fortran )
link_directories(${CMAKE_SOURCE_DIR}/libs)

add_executable (dummy_solver_fortran ${CMAKE_CURRENT_SOURCE_DIR}/dummy_solver.f90)
target_link_libraries(dummy_solver_fortran co_sim_c_io)
