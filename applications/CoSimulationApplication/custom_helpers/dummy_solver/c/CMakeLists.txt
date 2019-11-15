project( DummySolverC )

set(KRATOS_CO_SIMULATION_DUMMY_SOLVER_C_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/dummy_solver.c)

add_executable (dummy_solver_c ${KRATOS_CO_SIMULATION_DUMMY_SOLVER_C_SOURCES})
