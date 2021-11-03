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

class CO_SIM_IO_API FileCommunication : public Communication
{
public:
    FileCommunication(
        const Info& I_Settings,
        std::shared_ptr<DataCommunicator> I_DataComm);

    ~FileCommunication() override;

private:
    std::string GetCommunicationName() const override {return "file";}

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
};

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_FILE_COMMUNICATION_INCLUDED
