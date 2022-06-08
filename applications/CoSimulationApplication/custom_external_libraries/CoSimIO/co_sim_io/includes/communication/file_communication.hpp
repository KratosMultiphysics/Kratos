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

#ifndef CO_SIM_IO_FILE_COMMUNICATION_INCLUDED
#define CO_SIM_IO_FILE_COMMUNICATION_INCLUDED

// System includes

// Project includes
#include "communication.hpp"
#include "includes/filesystem_inc.hpp"

namespace CoSimIO {
namespace Internals {

#ifdef CO_SIM_IO_COMPILED_IN_WINDOWS
constexpr bool USE_AUX_FILE_FOR_FILE_AVAILABILITY = true;
#else
constexpr bool USE_AUX_FILE_FOR_FILE_AVAILABILITY = false;
#endif

class CO_SIM_IO_API FileCommunication : public Communication
{
public:
    FileCommunication(
        const Info& I_Settings,
        std::shared_ptr<DataCommunicator> I_DataComm);

    ~FileCommunication() override;

private:
    bool mUseAuxFileForFileAvailability = USE_AUX_FILE_FOR_FILE_AVAILABILITY;
    const bool mUseFileSerializer = true;

    std::string GetCommunicationName() const override {return "file";}

    void DerivedHandShake() const override;

    Info GetCommunicationSettings() const override;

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

    template<class TObjectType>
    Info GenericSendWithFileSerializer(
        const Info& I_Info,
        const TObjectType& rObj);

    template<class TObjectType>
    Info GenericReceiveWithFileSerializer(
        const Info& I_Info,
        TObjectType& rObj);

    template<typename T>
    double GenericSend(
        const Info& I_Info,
        const T& rData,
        const int SizeOfData);

    template<typename T>
    double GenericReceive(
        const Info& I_Info,
        T& rData,
        const int SizeOfData);

    double SendString(
        const Info& I_Info,
        const std::string& rData) override;

    double ReceiveString(
        const Info& I_Info,
        std::string& rData) override;

    double SendDataContainer(
        const Info& I_Info,
        const Internals::DataContainer<double>& rData) override;

    double ReceiveDataContainer(
        const Info& I_Info,
        Internals::DataContainer<double>& rData) override;
};

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_FILE_COMMUNICATION_INCLUDED
