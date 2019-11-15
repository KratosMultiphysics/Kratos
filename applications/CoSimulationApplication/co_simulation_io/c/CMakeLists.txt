project( CoSimCIO C )

set(KRATOS_CO_SIM_C_IO_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/co_sim_c_io.cpp)

add_library (co_sim_c_io SHARED ${KRATOS_CO_SIM_C_IO_SOURCES})
