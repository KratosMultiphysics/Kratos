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

#ifndef CO_SIM_IO_MPI_COMMUNICATION_H_INCLUDED
#define CO_SIM_IO_MPI_COMMUNICATION_H_INCLUDED

// System includes
#include "mpi.h"

// Project includes
#include "communication.hpp"

namespace CoSimIO {
namespace Internals {

class MPICommunication : public Communication
{
public:
    explicit MPICommunication(const std::string& rName, SettingsType& rSettings, const bool IsConnectionMaster)
        : Communication(rName, rSettings, IsConnectionMaster)
    {
       static_assert(false,"MPI Communication is not implemented yet");
       /*
        Note to self:
        If I directly use the buffer of the sender, then I have to ensure it can be reused when returning from the sending function. This can be achieved with two variants:
        - Use blocking communication. This might block the sender for some time
        - Use "buffered" non-blocking communication: => before calling Isend copy the values to a local buffer. Then wait for the send to complete => this is somehow nasty, could prob only be done in a separate thread??? Or I somehow save it internally and check if in the next time an IO function is called if the send is completed ...
       */
    }
};

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_MPI_COMMUNICATION_H_INCLUDED
