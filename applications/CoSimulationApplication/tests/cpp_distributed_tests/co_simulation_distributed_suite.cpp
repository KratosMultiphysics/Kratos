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

// External includes
#include <gmock/gmock.h>
#include <gtest/gtest.h>

// Project includes
#include "co_simulation_distributed_suite.h"

namespace Kratos::Testing {

KratosCoSimulationMPIFastSuite::KratosCoSimulationMPIFastSuite()
    : KratosMPICoreFastSuite() {
  mpCoSimulationApp = std::make_shared<KratosCoSimulationApplication>();
  this->ImportApplicationIntoKernel(mpCoSimulationApp);
}

} // namespace Kratos::Testing