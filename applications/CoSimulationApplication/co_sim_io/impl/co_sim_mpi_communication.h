// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifndef KRATOS_CO_SIM_MPI_COMM_H_INCLUDED
#define KRATOS_CO_SIM_MPI_COMM_H_INCLUDED

// System includes
#include "mpi.h"

// Project includes
#include "co_sim_communication.h"

namespace CoSimIO {
namespace Internals {

class CoSimMPICommunication : public CoSimCommunication
{
public:
    explicit CoSimMPICommunication(const std::string& rName, SettingsType& rSettings, const bool IsConnectionMaster)
        : CoSimCommunication(rName, rSettings, IsConnectionMaster)
    {
       KRATOS_CO_SIM_ERROR << "MPI Communication is not implemented yet" << std::endl;
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

#endif /* KRATOS_CO_SIM_MPI_COMM_H_INCLUDED */