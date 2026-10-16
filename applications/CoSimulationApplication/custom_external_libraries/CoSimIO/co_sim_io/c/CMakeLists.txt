if(WIN32)
  # To enable WINDOWS_EXPORT_ALL_SYMBOLS to automatically
  # add export decorators for all classes, structs and functions
  # see https://cmake.org/cmake/help/latest/prop_tgt/WINDOWS_EXPORT_ALL_SYMBOLS.html#prop_tgt:WINDOWS_EXPORT_ALL_SYMBOLS
  cmake_minimum_required(VERSION 3.4)
endif(WIN32)

message("Configuring CoSimIO for C")

add_library (co_sim_io_c SHARED co_sim_io_c.cpp)
set_target_properties(co_sim_io_c PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)

target_link_libraries( co_sim_io_c co_sim_io )

install(TARGETS co_sim_io_c DESTINATION bin)

if (CO_SIM_IO_BUILD_MPI)
    add_library (co_sim_io_c_mpi SHARED co_sim_io_c_mpi.cpp)

    target_link_libraries(co_sim_io_c_mpi co_sim_io_c co_sim_io_mpi)

    install(TARGETS co_sim_io_c_mpi DESTINATION bin)
endif()
