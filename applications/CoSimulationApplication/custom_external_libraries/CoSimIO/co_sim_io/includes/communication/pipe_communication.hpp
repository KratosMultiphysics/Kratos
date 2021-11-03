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

#ifndef CO_SIM_IO_PIPE_COMMUNICATION_INCLUDED
#define CO_SIM_IO_PIPE_COMMUNICATION_INCLUDED

// System includes
#include <unordered_map>

// Project includes
#include "communication.hpp"

namespace CoSimIO {
namespace Internals {

class CO_SIM_IO_API PipeCommunication : public Communication
{
public:
    PipeCommunication(
        const Info& I_Settings,
        std::shared_ptr<DataCommunicator> I_DataComm);

    ~PipeCommunication() override
    {
        if (GetIsConnected()) {
            CO_SIM_IO_INFO("CoSimIO") << "Warning: Disconnect was not performed, attempting automatic disconnection!" << std::endl;
            Info tmp;
            Disconnect(tmp);
        }
    }

private:
class BidirectionalPipe
{
public:
    BidirectionalPipe(
        const fs::path& rPipeDir,
        const fs::path& rBasePipeName,
        const bool IsPrimary);

    void Write(const std::string& rData);

    void Read(std::string& rData);

    template<class TObjectType>
    void Send(const TObjectType& rObject)
    {
        StreamSerializer serializer;
        serializer.save("object", rObject);

        Write(serializer.GetStringRepresentation());
    }

    template<class TObjectType>
    void Receive(TObjectType& rObject)
    {
        std::string buffer;
        Read(buffer);
        StreamSerializer serializer(buffer);
        serializer.load("object", rObject);
    }

    void Close();


private:
    int mPipeHandleWrite;
    int mPipeHandleRead;

    fs::path mPipeNameWrite;
    fs::path mPipeNameRead;

    void SendSize(const std::uint64_t Size);

    std::uint64_t ReceiveSize();
};

    std::shared_ptr<BidirectionalPipe> mpPipe;

    std::string GetCommunicationName() const override {return "pipe";}

    Info ConnectDetail(const Info& I_Info) override;

    Info DisconnectDetail(const Info& I_Info) override;

    Info ImportInfoImpl(const Info& I_Info) override;

    Info ExportInfoImpl(const Info& I_Info) override;

    Info ImportDataImpl(
        const Info& I_Info,
        Internals::DataContainer<double>& rData) override;

    Info ExportDataImpl(
        const Info& I_Info,
        const Internals::DataContainer<double>& rData) override;

    Info ImportMeshImpl(
        const Info& I_Info,
        ModelPart& O_ModelPart) override;

    Info ExportMeshImpl(
        const Info& I_Info,
        const ModelPart& I_ModelPart) override;

    void DerivedHandShake() override;
};

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_PIPE_COMMUNICATION_INCLUDED
