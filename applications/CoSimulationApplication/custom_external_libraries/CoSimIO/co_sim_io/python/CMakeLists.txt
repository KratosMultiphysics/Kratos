message("Configuring CoSimIO for Python")

# TODO make this configureable
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../external_libraries/pybind11 dummy)

pybind11_add_module(PyCoSimIO co_sim_io_python.cpp)

target_link_libraries( PyCoSimIO PRIVATE co_sim_io )

install(TARGETS PyCoSimIO DESTINATION bin)
file(WRITE "${CMAKE_INSTALL_PREFIX}/bin/CoSimIO/__init__.py" "from PyCoSimIO import *\nfrom PyCoSimIO import __version__\n")

set(CoSimIOMPIInitFile "${CMAKE_INSTALL_PREFIX}/bin/CoSimIO/mpi/__init__.py")
set(mpi4pyInterfaceInitFile "${CMAKE_INSTALL_PREFIX}/bin/CoSimIO/mpi/mpi4pyInterface/__init__.py")

# dummy init files that give proper errors in case something related to MPI was not compiled
# These files will be overwritten if the corresponding option is enabled, hence here we can write then unconditionally
file(WRITE ${CoSimIOMPIInitFile} "raise Exception('CoSimIO was compiled without MPI support! (use \"CO_SIM_IO_BUILD_MPI\" to enable it)')\n")
file(WRITE ${mpi4pyInterfaceInitFile} "raise Exception('The mpi4py interface was not compiled! (use \"CO_SIM_IO_BUILD_PYTHON_MPI4PY_INTERFACE\" to enable it)')\n")

if (CO_SIM_IO_BUILD_MPI)
    pybind11_add_module(PyCoSimIOMPI co_sim_io_python_mpi.cpp)

    target_link_libraries( PyCoSimIOMPI PRIVATE co_sim_io co_sim_io_mpi )

    install(TARGETS PyCoSimIOMPI DESTINATION bin)
    file(WRITE ${CoSimIOMPIInitFile} "from PyCoSimIOMPI import *\n")

    OPTION ( CO_SIM_IO_BUILD_PYTHON_MPI4PY_INTERFACE "Building the interface to mpi4py MPI communicators" OFF )
    if (CO_SIM_IO_BUILD_PYTHON_MPI4PY_INTERFACE)
        message("Building the interface to mpi4py MPI communicators")
        pybind11_add_module(PyCoSimIOMPI_mpi4pyInterface mpi4py_interface.cpp)
        target_link_libraries( PyCoSimIOMPI_mpi4pyInterface PRIVATE ${MPI_LIBRARIES} )
        install(TARGETS PyCoSimIOMPI_mpi4pyInterface DESTINATION bin)
        file(WRITE ${mpi4pyInterfaceInitFile} "from PyCoSimIOMPI_mpi4pyInterface import *\n")
    endif()
endif()
