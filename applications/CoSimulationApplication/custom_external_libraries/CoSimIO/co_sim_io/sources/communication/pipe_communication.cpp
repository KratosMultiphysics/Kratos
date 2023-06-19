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

namespace {

int GetPipeBufferSize(const Info& I_Info)
{
    int default_buffer_size = 8192;

    #ifdef CO_SIM_IO_COMPILED_IN_LINUX
    default_buffer_size = 65536;
    #endif

    return I_Info.Get<int>("buffer_size", default_buffer_size);
}

} // anonymous namespace

PipeCommunication::PipeCommunication(
    const Info& I_Settings,
    std::shared_ptr<DataCommunicator> I_DataComm)
    : Communication(I_Settings, I_DataComm),
      mBufferSize(GetPipeBufferSize(I_Settings))
{
}

PipeCommunication::~PipeCommunication()
{
    if (GetIsConnected()) {
        CO_SIM_IO_INFO("CoSimIO") << "Warning: Disconnect was not performed, attempting automatic disconnection!" << std::endl;
        Info tmp;
        Disconnect(tmp);
    }
}

Info PipeCommunication::ConnectDetail(const Info& I_Info)
{
    CO_SIM_IO_TRY

    CO_SIM_IO_INFO_IF("CoSimIO", GetDataCommunicator().IsDistributed() && GetDataCommunicator().Rank()==0) << "Warning: Connection was done with MPI, but pipe based communication works only within the same machine. Communicating between different compute nodes in a distributed memory machine when does not work, it will hang!" << std::endl;

    mpPipe = std::make_shared<BidirectionalPipe>(
        GetCommunicationDirectory(),
        GetConnectionName() + "_r" + std::to_string(GetDataCommunicator().Rank()),
        GetIsPrimaryConnection(),
        GetPipeBufferSize(I_Info),
        GetEchoLevel());

    return Info(); // TODO use

    CO_SIM_IO_CATCH
}

Info PipeCommunication::DisconnectDetail(const Info& I_Info)
{
    mpPipe->Close();
    return Info(); // TODO use
}

void PipeCommunication::DerivedHandShake() const
{
    CO_SIM_IO_ERROR_IF(GetMyInfo().Get<std::string>("operating_system") != GetPartnerInfo().Get<std::string>("operating_system")) << "Pipe communication cannot be used between different operating systems!" << std::endl;

    const std::size_t my_use_buffer_size = GetMyInfo().Get<Info>("communication_settings").Get<std::size_t>("buffer_size");
    const std::size_t partner_buffer_size = GetPartnerInfo().Get<Info>("communication_settings").Get<std::size_t>("buffer_size");
    CO_SIM_IO_ERROR_IF(my_use_buffer_size != partner_buffer_size) << "Mismatch in buffer_size!\nMy buffer_size: " << my_use_buffer_size << "\nPartner buffer_size: " << partner_buffer_size << std::noboolalpha << std::endl;
}

Info PipeCommunication::GetCommunicationSettings() const
{
    CO_SIM_IO_TRY

    Info info;
    info.Set("buffer_size", mBufferSize);
    return info;

    CO_SIM_IO_CATCH
}


PipeCommunication::BidirectionalPipe::BidirectionalPipe(
    const fs::path& rPipeDir,
    const fs::path& rBasePipeName,
    const bool IsPrimary,
    const int BufferSize,
    const int EchoLevel) : mBufferSize(BufferSize)
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

    // if possible try to resize the pipes to the buffer size
    // not resizing still works but leads to communication in more chunks
    #if defined(F_GETPIPE_SZ) && defined(F_SETPIPE_SZ) // some old kernels don't define this
    const int pipe_buffer_size_read = fcntl(mPipeHandleRead, F_GETPIPE_SZ);
    const int pipe_buffer_size_write = fcntl(mPipeHandleWrite, F_GETPIPE_SZ);

    const int max_pipe_buffer = std::max(pipe_buffer_size_read, pipe_buffer_size_write);
    if (pipe_buffer_size_read  < max_pipe_buffer) {fcntl(mPipeHandleRead,  F_SETPIPE_SZ, BufferSize);}
    if (pipe_buffer_size_write < max_pipe_buffer) {fcntl(mPipeHandleWrite, F_SETPIPE_SZ, BufferSize);}

    if (BufferSize > max_pipe_buffer) {
        CO_SIM_IO_INFO_IF("CoSimIO", EchoLevel>0) << "Requested buffer size (" << BufferSize << ") is larger than pipe buffer size (" << max_pipe_buffer << "). Attempting to increase size of pipe buffer" << std::endl;
        fcntl(mPipeHandleRead,  F_SETPIPE_SZ, BufferSize);
        fcntl(mPipeHandleWrite, F_SETPIPE_SZ, BufferSize);

        const int new_pipe_buffer_size = fcntl(mPipeHandleRead, F_GETPIPE_SZ);
        CO_SIM_IO_ERROR_IF(new_pipe_buffer_size != fcntl(mPipeHandleWrite, F_GETPIPE_SZ)) << "Different buffer sizes after changing size, this should not happen!" << std::endl;

        // not comparing equal, as pipe buffer size is multiple of getpagesize()
        CO_SIM_IO_INFO_IF("CoSimIO", new_pipe_buffer_size == BufferSize && EchoLevel>0) << "Resizing pipe buffer was successful! Pipe buffer size is now " << new_pipe_buffer_size << std::endl;
        CO_SIM_IO_INFO_IF("CoSimIO", new_pipe_buffer_size < BufferSize) << "Resizing pipe buffer was not successful! Pipe buffer size is " << new_pipe_buffer_size << " even though " << BufferSize << " was requested!" << std::endl;
    }

    const int final_pipe_buffer_size = fcntl(mPipeHandleRead, F_GETPIPE_SZ);
    if (BufferSize >= final_pipe_buffer_size) {mBufferSize = final_pipe_buffer_size-1;} // crashes if same size or larger!
    #endif

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
    const ssize_t bytes_written = write(mPipeHandleWrite, &Size, sizeof(Size));
    CO_SIM_IO_ERROR_IF(bytes_written < 0) << "Error in writing to Pipe!" << std::endl;
    #endif
}

std::uint64_t PipeCommunication::BidirectionalPipe::ReceiveSize()
{
    #ifndef CO_SIM_IO_COMPILED_IN_WINDOWS
    std::uint64_t imp_size_u;
    const ssize_t bytes_read = read(mPipeHandleRead, &imp_size_u, sizeof(imp_size_u));
    CO_SIM_IO_ERROR_IF(bytes_read < 0) << "Error in reading from Pipe!" << std::endl;
    return imp_size_u;
    #else
    return 0;
    #endif
}

double PipeCommunication::SendString(
    const Info& I_Info,
    const std::string& rData)
{
    return mpPipe->Write(rData, 1);
}

double PipeCommunication::ReceiveString(
    const Info& I_Info,
    std::string& rData)
{
    return mpPipe->Read(rData, 1);
}

double PipeCommunication::SendDataContainer(
    const Info& I_Info,
    const Internals::DataContainer<double>& rData)
{
    return mpPipe->Write(rData, sizeof(double));
}

double PipeCommunication::ReceiveDataContainer(
    const Info& I_Info,
    Internals::DataContainer<double>& rData)
{
    return mpPipe->Read(rData, sizeof(double));
}

} // namespace Internals
} // namespace CoSimIO
