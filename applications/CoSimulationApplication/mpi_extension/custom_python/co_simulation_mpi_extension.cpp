//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "mpi/includes/mpi_data_communicator.h"

// CoSimIO
#include "custom_external_libraries/CoSimIO/co_sim_io/co_sim_io_mpi.hpp"

namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosCoSimulationMPIExtension,m)
{
    namespace py = pybind11;

    auto m_co_sim_io = m.def_submodule("CoSimIO");

    m_co_sim_io.def("ConnectMPI",    [](const CoSimIO::Info& Settings, const DataCommunicator& rDataComm){
        KRATOS_ERROR_IF_NOT(rDataComm.IsDistributed()) << "ConnectMPI requires a distributed DataCommunicator!" << std::endl;
        KRATOS_ERROR_IF_NOT(rDataComm.IsDefinedOnThisRank()) << "This rank is not part of this MPI_Comm!" << std::endl;
        MPI_Comm the_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rDataComm);
        return CoSimIO::ConnectMPI(Settings, the_mpi_comm);
    });
}

}
}
