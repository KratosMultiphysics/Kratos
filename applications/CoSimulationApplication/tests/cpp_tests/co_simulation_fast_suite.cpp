//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Richard Faasse
//

#include "co_simulation_fast_suite.h"

namespace Kratos::Testing {
KratosCoSimulationFastSuite::KratosCoSimulationFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("CoSimulationApplication")) {
    mpCoSimulationApp = std::make_shared<KratosCoSimulationApplication>();
    this->mKernel.ImportApplication(mpCoSimulationApp);
  }
}

} // namespace Kratos::Testing
