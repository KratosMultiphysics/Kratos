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

#pragma once

#include "co_simulation_application.h"
#include "mpi/testing/mpi_testing.h"

namespace Kratos::Testing {

class KratosCoSimulationMPIFastSuite : public KratosMPICoreFastSuite {
public:
  KratosCoSimulationMPIFastSuite();

private:
  KratosCoSimulationApplication::Pointer mpCoSimulationApp;
};

} // namespace Kratos::Testing
