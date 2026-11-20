//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// Exposure of the CoSimIO MPI interface to Python

// System includes
// pybind includes
#include <pybind11/pybind11.h>

// CoSimIO includes
#include "co_sim_io_mpi.hpp"
#include "mpi_comm_holder.hpp"

PYBIND11_MODULE(PyCoSimIOMPI, m)
{
    namespace py = pybind11;

    m.def("ConnectMPI", [](const CoSimIO::Info& I_Info)
        { return CoSimIO::ConnectMPI(I_Info, MPI_COMM_WORLD); });

    m.def("ConnectMPI", [](const CoSimIO::Info& I_Info, const CoSimIO::MPICommHolder& holder)
        {
            return CoSimIO::ConnectMPI(I_Info, holder.GetMPIComm());
        });

    py::class_<CoSimIO::MPICommHolder, std::shared_ptr<CoSimIO::MPICommHolder>>(m,"MPICommHolder");
}
