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

#ifndef CO_SIM_IO_COMMUNICATION_FACTORY_INCLUDED
#define CO_SIM_IO_COMMUNICATION_FACTORY_INCLUDED

// System includes
#include <unordered_map>
#include <functional>

// Project includes
#include "communication.hpp"
#include "includes/info.hpp"
#include "includes/data_communicator.hpp"

namespace CoSimIO {
namespace Internals {

class CommunicationFactory
{
public:
    std::unique_ptr<Communication> CO_SIM_IO_API Create(
        const Info& I_Settings,
        const std::shared_ptr<DataCommunicator> pDataComm) const;

protected:
    using CommCreateFctType = std::function<std::unique_ptr<Communication>(const Info&, const std::shared_ptr<DataCommunicator>)>;
    using CommCreateFctsType = std::unordered_map<std::string, CommCreateFctType>;

    virtual CommCreateFctsType CO_SIM_IO_API GetCommunicationCreateFunctions() const;

private:
    virtual bool CO_SIM_IO_API IsMPI() const {return false;}
};

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_COMMUNICATION_FACTORY_INCLUDED
