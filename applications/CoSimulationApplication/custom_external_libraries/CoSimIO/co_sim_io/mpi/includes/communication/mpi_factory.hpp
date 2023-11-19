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

#ifndef CO_SIM_IO_MPI_COMMUNICATION_FACTORY_INCLUDED
#define CO_SIM_IO_MPI_COMMUNICATION_FACTORY_INCLUDED

// System includes

// Project includes
#include "includes/communication/factory.hpp"

namespace CoSimIO {
namespace Internals {

class MPICommunicationFactory : public CommunicationFactory
{

protected:
    using CommCreateFctType = CommunicationFactory::CommCreateFctType;
    using CommCreateFctsType = CommunicationFactory::CommCreateFctsType;

    CommCreateFctsType GetCommunicationCreateFunctions() const override;

private:
    bool IsMPI() const override {return true;}
};

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_MPI_COMMUNICATION_FACTORY_INCLUDED
