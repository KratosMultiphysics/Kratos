project( DummySolverFortran Fortran )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../../../co_simulation_io/fortran )

add_executable (dummy_solver_fortran ${CMAKE_CURRENT_SOURCE_DIR}/dummy_solver.f90)
