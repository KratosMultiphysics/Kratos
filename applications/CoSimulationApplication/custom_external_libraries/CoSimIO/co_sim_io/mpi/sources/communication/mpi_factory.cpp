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

// System includes

// Project includes
#include "mpi/includes/communication/mpi_factory.hpp"

#include "mpi/includes/communication/mpi_inter_communication.hpp"

namespace CoSimIO {
namespace Internals {

CommunicationFactory::CommCreateFctsType MPICommunicationFactory::GetCommunicationCreateFunctions() const
{
    CO_SIM_IO_TRY

    auto fcts = CommunicationFactory::GetCommunicationCreateFunctions();

#ifdef CO_SIM_IO_BUILD_MPI_COMMUNICATION
    fcts["mpi_inter"] = [](
        const Info& I_Settings,
        const std::shared_ptr<DataCommunicator> pDataComm){
            return CoSimIO::make_unique<MPIInterCommunication>(I_Settings, pDataComm);};
#else
    fcts["mpi_inter"] = [](
        const Info& I_Settings,
        const std::shared_ptr<DataCommunicator> pDataComm){
            CO_SIM_IO_ERROR << "Communication via MPI must be enabled at compile time with \"CO_SIM_IO_BUILD_MPI_COMMUNICATION\"!" << std::endl;
            return nullptr;};

#endif // CO_SIM_IO_BUILD_MPI_COMMUNICATION

    return fcts;

    CO_SIM_IO_CATCH
}

} // namespace Internals
} // namespace CoSimIO
