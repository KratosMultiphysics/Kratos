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

#include "convection_diffusion_fast_suite.h"

namespace Kratos::Testing {
KratosConvectionDiffusionFastSuite::KratosConvectionDiffusionFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("ConvectionDiffusionApplication")) {
    mpConvectionDiffusionApp = std::make_shared<KratosConvectionDiffusionApplication>();
    this->mKernel.ImportApplication(mpConvectionDiffusionApp);
  }
}

} // namespace Kratos::Testing
