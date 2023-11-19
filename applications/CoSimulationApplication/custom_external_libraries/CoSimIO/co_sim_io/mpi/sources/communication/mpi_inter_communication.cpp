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
#include "mpi/includes/communication/mpi_inter_communication.hpp"
#include "mpi/includes/mpi_data_communicator.hpp"

namespace CoSimIO {
namespace Internals {

#ifdef CO_SIM_IO_BUILD_MPI_COMMUNICATION

namespace {

template<typename TMPIDataType>
int ReceiveSize(
    MPI_Comm Comm,
    TMPIDataType DataType,
    const int Rank)
{
    int size;
    MPI_Status status;
    MPI_Probe(Rank, 0, Comm, &status);
    MPI_Get_count(&status, DataType, &size);
    return size;
}

}

MPIInterCommunication::MPIInterCommunication(
    const Info& I_Settings,
    std::shared_ptr<DataCommunicator> I_DataComm)
    : Communication(I_Settings, I_DataComm)
{
    CO_SIM_IO_ERROR_IF_NOT(I_DataComm->IsDistributed()) << "MPI communication only works with a MPIDataCommunicator!" << std::endl;
}

MPIInterCommunication::~MPIInterCommunication()
{
    if (GetIsConnected()) {
        CO_SIM_IO_INFO("CoSimIO") << "Warning: Disconnect was not performed, attempting automatic disconnection!" << std::endl;
        Info tmp;
        Disconnect(tmp);
    }
}

Info MPIInterCommunication::ConnectDetail(const Info& I_Info)
{
    CO_SIM_IO_TRY

    if (!GetIsPrimaryConnection()) {
        mPortName = GetPartnerInfo().Get<Info>("communication_settings").Get<std::string>("port_name");
    }

    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1 && GetDataCommunicator().Rank()==0) << "Using MPI-port: " << mPortName << std::endl;

    MPI_Comm my_comm = MPIDataCommunicator::GetMPICommunicator(GetDataCommunicator());

    if (GetIsPrimaryConnection()) {
        MPI_Comm_accept(mPortName.c_str(), MPI_INFO_NULL, 0, my_comm, &mInterComm); // todo check return code
    } else {
        MPI_Comm_connect(mPortName.c_str(), MPI_INFO_NULL, 0, my_comm, &mInterComm); // todo check return code
    }

    return Info(); // TODO use

    CO_SIM_IO_CATCH
}

Info MPIInterCommunication::DisconnectDetail(const Info& I_Info)
{
    CO_SIM_IO_TRY

    MPI_Comm_disconnect(&mInterComm); // todo check return code

    if (GetIsPrimaryConnection() && GetDataCommunicator().Rank()==0) {
        MPI_Close_port(mPortName.c_str()); // todo check return code
    }

    return Info(); // TODO use

    CO_SIM_IO_CATCH
}
void MPIInterCommunication::PrepareConnection(const Info& I_Info)
{
    CO_SIM_IO_TRY

    if (GetIsPrimaryConnection()) {
        const auto& r_data_comm = GetDataCommunicator();
        mPortName.resize(MPI_MAX_PORT_NAME);

        if (r_data_comm.Rank()==0) {
            MPI_Open_port(MPI_INFO_NULL, &mPortName.front()); // todo check return code
        }

        r_data_comm.Broadcast(mPortName, 0);

        mPortName.erase(mPortName.find('\0'));
    }

    CO_SIM_IO_CATCH
}

Info MPIInterCommunication::GetCommunicationSettings() const
{
    CO_SIM_IO_TRY

    Info info;

    if (GetIsPrimaryConnection() && GetDataCommunicator().Rank() == 0) {
        info.Set("port_name", mPortName);
    }

    return info;

    CO_SIM_IO_CATCH
}

double MPIInterCommunication::SendString(
    const Info& I_Info,
    const std::string& rData)
{
    CO_SIM_IO_TRY

    const auto start_time(std::chrono::steady_clock::now());

    MPI_Send(
        rData.data(),
        rData.size(),
        MPI_CHAR,
        GetDataCommunicator().Rank(),
        0,
        mInterComm); // todo check return code

    return Utilities::ElapsedSeconds(start_time);

    CO_SIM_IO_CATCH
}

double MPIInterCommunication::ReceiveString(
    const Info& I_Info,
    std::string& rData)
{
    CO_SIM_IO_TRY

    const int size = ReceiveSize(mInterComm, MPI_CHAR, GetDataCommunicator().Rank()); // serves also as synchronization for time measurement
    rData.resize(size);

    const auto start_time(std::chrono::steady_clock::now());

    MPI_Recv(
        &(rData.front()),
        rData.size(),
        MPI_CHAR,
        GetDataCommunicator().Rank(),
        0,
        mInterComm,
        MPI_STATUS_IGNORE); // todo check return code

    return Utilities::ElapsedSeconds(start_time);

    CO_SIM_IO_CATCH
}

double MPIInterCommunication::SendDataContainer(
    const Info& I_Info,
    const Internals::DataContainer<double>& rData)
{
    CO_SIM_IO_TRY

    const auto start_time(std::chrono::steady_clock::now());

    MPI_Send(
        rData.data(),
        rData.size(),
        MPI_DOUBLE,
        GetDataCommunicator().Rank(),
        0,
        mInterComm); // todo check return code

    return Utilities::ElapsedSeconds(start_time);

    CO_SIM_IO_CATCH
}

double MPIInterCommunication::ReceiveDataContainer(
    const Info& I_Info,
    Internals::DataContainer<double>& rData)
{
    CO_SIM_IO_TRY

    const int size = ReceiveSize(mInterComm, MPI_DOUBLE, GetDataCommunicator().Rank()); // serves also as synchronization for time measurement
    rData.resize(size);

    const auto start_time(std::chrono::steady_clock::now());

    MPI_Recv(
        rData.data(),
        rData.size(),
        MPI_DOUBLE,
        GetDataCommunicator().Rank(),
        0,
        mInterComm,
        MPI_STATUS_IGNORE); // todo check return code

    return Utilities::ElapsedSeconds(start_time);

    CO_SIM_IO_CATCH
}

#endif // CO_SIM_IO_BUILD_MPI_COMMUNICATION

} // namespace Internals
} // namespace CoSimIO
