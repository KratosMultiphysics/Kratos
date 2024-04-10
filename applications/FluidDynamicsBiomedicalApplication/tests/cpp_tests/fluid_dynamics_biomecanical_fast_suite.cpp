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

#include "fluid_dynamics_biomecanical_fast_suite.h"

namespace Kratos::Testing {
KratosFluidDynamicsBiomecanicalFastSuite::KratoFluidDynamicsBiomecanicalFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("FluidDynamicsBiomecanicalApplication")) {
    mpFluidDynamicsBiomecanicalApp = std::make_shared<KratosFluidDynamicsBiomecanicalApplication>();
    this->mKernel.ImportApplication(mpFluidDynamicsBiomecanicalApp);
  }
}

} // namespace Kratos::Testing
