project( DummySolverFortran Fortran )

set(KRATOS_CO_SIMULATION_DUMMY_SOLVER_FORTRAN_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/dummy_solver.f90)

add_executable (dummy_solver_fortran ${KRATOS_CO_SIMULATION_DUMMY_SOLVER_FORTRAN_SOURCES})
