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

#include "fluid_dynamics_fast_suite.h"

namespace Kratos::Testing {
KratosFluidDynamicsFastSuite::KratosFluidDynamicsFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("FluidDynamicsApplication")) {
    mpFluidDynamicsApp = std::make_shared<KratosFluidDynamicsApplication>();
    this->mKernel.ImportApplication(mpFluidDynamicsApp);
  }
}

} // namespace Kratos::Testing
