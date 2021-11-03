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
#include "includes/define.hpp" // for "CO_SIM_IO_COMPILED_IN_WINDOWS"

#ifdef CO_SIM_IO_COMPILED_IN_WINDOWS

#else
    #include <fcntl.h>
    #include <sys/stat.h>
    #include <sys/types.h>
    #include <unistd.h>
#endif

// External includes

// Project includes
#include "includes/communication/pipe_communication.hpp"
#include "includes/stream_serializer.hpp"
#include "includes/utilities.hpp"

namespace CoSimIO {
namespace Internals {

PipeCommunication::PipeCommunication(
    const Info& I_Settings,
    std::shared_ptr<DataCommunicator> I_DataComm)
    : Communication(I_Settings, I_DataComm)
{
}

Info PipeCommunication::ConnectDetail(const Info& I_Info)
{
    mpPipe = std::make_shared<BidirectionalPipe>(
        GetCommunicationDirectory(),
        GetConnectionName() + "_r" + std::to_string(GetDataCommunicator().Rank()),
        GetIsPrimaryConnection());

    return Info(); // TODO use
}

Info PipeCommunication::DisconnectDetail(const Info& I_Info)
{
    mpPipe->Close();
    return Info(); // TODO use
}

Info PipeCommunication::ImportInfoImpl(const Info& I_Info)
{
    Info imported_info;
    mpPipe->Receive(imported_info);
    return imported_info;
}

Info PipeCommunication::ExportInfoImpl(const Info& I_Info)
{
    mpPipe->Send(I_Info);
    return Info(); // TODO use
}

Info PipeCommunication::ImportDataImpl(
    const Info& I_Info,
    Internals::DataContainer<double>& rData)
{
    mpPipe->Receive(rData);
    return Info(); // TODO use
}

Info PipeCommunication::ExportDataImpl(
    const Info& I_Info,
    const Internals::DataContainer<double>& rData)
{
    mpPipe->Send(rData);
    return Info(); // TODO use
}

Info PipeCommunication::ImportMeshImpl(
    const Info& I_Info,
    ModelPart& O_ModelPart)
{
    mpPipe->Receive(O_ModelPart);
    return Info(); // TODO use
}

Info PipeCommunication::ExportMeshImpl(
    const Info& I_Info,
    const ModelPart& I_ModelPart)
{
    mpPipe->Send(I_ModelPart);
    return Info(); // TODO use
}

void PipeCommunication::DerivedHandShake()
{
    CO_SIM_IO_ERROR_IF(GetMyInfo().Get<std::string>("operating_system") != GetPartnerInfo().Get<std::string>("operating_system")) << "Pipe communication cannot be used between different operating systems!" << std::endl;
}


PipeCommunication::BidirectionalPipe::BidirectionalPipe(
    const fs::path& rPipeDir,
    const fs::path& rBasePipeName,
    const bool IsPrimary)
{
    mPipeNameWrite = mPipeNameRead = rPipeDir / rBasePipeName;

    #ifdef CO_SIM_IO_COMPILED_IN_WINDOWS
    CO_SIM_IO_ERROR << "Pipe communication is not yet implemented for Windows!" << std::endl;
    #else
    if (IsPrimary) {
        mPipeNameWrite += "_p2s";
        mPipeNameRead  += "_s2p";

        CO_SIM_IO_ERROR_IF(mkfifo(mPipeNameWrite.c_str(), 0666) != 0) << "Pipe " << mPipeNameWrite << " could not be created!" << std::endl;
        CO_SIM_IO_ERROR_IF(mkfifo(mPipeNameRead.c_str(), 0666) != 0) << "Pipe " << mPipeNameRead << " could not be created!" << std::endl;
        CO_SIM_IO_ERROR_IF((mPipeHandleWrite = open(mPipeNameWrite.c_str(), O_WRONLY)) < 0) << "Pipe " << mPipeNameWrite << " could not be opened!" << std::endl;
        CO_SIM_IO_ERROR_IF((mPipeHandleRead = open(mPipeNameRead.c_str(), O_RDONLY)) < 0) << "Pipe " << mPipeNameRead << " could not be opened!" << std::endl;

    } else {
        mPipeNameWrite += "_s2p";
        mPipeNameRead  += "_p2s";

        Utilities::WaitUntilPathExists(mPipeNameWrite); // created last hence wait for it

        CO_SIM_IO_ERROR_IF((mPipeHandleRead = open(mPipeNameRead.c_str(), O_RDONLY)) < 0) << "Pipe " << mPipeNameRead << " could not be opened!" << std::endl;
        CO_SIM_IO_ERROR_IF((mPipeHandleWrite = open(mPipeNameWrite.c_str(), O_WRONLY)) < 0) << "Pipe " << mPipeNameWrite << " could not be opened!" << std::endl;
    }
    #endif
}

void PipeCommunication::BidirectionalPipe::Write(const std::string& rData)
{
    #ifndef CO_SIM_IO_COMPILED_IN_WINDOWS
    SendSize(rData.size());
    write(mPipeHandleWrite, rData.c_str(), rData.size());
    #endif
}

void PipeCommunication::BidirectionalPipe::Read(std::string& rData)
{
    #ifndef CO_SIM_IO_COMPILED_IN_WINDOWS
    std::size_t received_size = ReceiveSize();
    rData.resize(received_size);
    read(mPipeHandleRead, &(rData.front()), received_size); // using front as other methods that access the underlying char are const
    #endif
}

void PipeCommunication::BidirectionalPipe::Close()
{
    #ifndef CO_SIM_IO_COMPILED_IN_WINDOWS
    close(mPipeHandleWrite);
    close(mPipeHandleRead);
    #endif
}

void PipeCommunication::BidirectionalPipe::SendSize(const std::uint64_t Size)
{
    #ifndef CO_SIM_IO_COMPILED_IN_WINDOWS
    write(mPipeHandleWrite, &Size, sizeof(Size));
    #endif
}

std::uint64_t PipeCommunication::BidirectionalPipe::ReceiveSize()
{
    #ifndef CO_SIM_IO_COMPILED_IN_WINDOWS
    std::uint64_t imp_size_u;
    read(mPipeHandleRead, &imp_size_u, sizeof(imp_size_u));
    return imp_size_u;
    #else
    return 0;
    #endif

}

} // namespace Internals
} // namespace CoSimIO
