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

#include "shallow_water_fast_suite.h"

namespace Kratos::Testing {
KratosShallowWaterFastSuite::KratosShallowWaterFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("ShallowWaterApplication")) {
    mpShallowWaterApp = std::make_shared<KratosShallowWaterApplication>();
    this->mKernel.ImportApplication(mpShallowWaterApp);
  }
}

} // namespace Kratos::Testing
