// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher
//

#ifndef KRATOS_CO_SIM_SOCKETS_COMM_H_INCLUDED
#define KRATOS_CO_SIM_SOCKETS_COMM_H_INCLUDED

// System includes

// Project includes
#include "co_sim_comm.h"

namespace CoSim {

class SocketsComm : public CoSimComm
{
public:
    SocketsComm(SettingsType& rSettings);

};

} // namespace CoSim

#endif /* KRATOS_CO_SIM_SOCKETS_COMM_H_INCLUDED */