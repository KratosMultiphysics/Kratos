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

// Project includes
extern "C" {
#include "co_sim_c_io.h"
}
#include "../co_sim_io.h"


void CoSimIO_Connect(const char* pName)
{
    CoSimIO::Connect(pName);
}

void CoSimIO_Disconnect(const char* pName)
{
    CoSimIO::Disconnect(pName);
}

void CoSimIO_ImportData(const char* pName)
{
    // CoSimIO::ImportData(pName);
}