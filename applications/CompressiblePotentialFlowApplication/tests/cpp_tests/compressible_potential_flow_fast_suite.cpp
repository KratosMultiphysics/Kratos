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

#include "compressible_potential_flow_fast_suite.h"

namespace Kratos::Testing {
KratosCompressiblePotentialFlowFastSuite::KratosCompressiblePotentialFlowFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("CompressiblePotentialFlowApplication")) {
    mpCompressiblePotentialFlowApp = std::make_shared<KratosCompressiblePotentialFlowApplication>();
    this->mKernel.ImportApplication(mpCompressiblePotentialFlowApp);
  }
}

} // namespace Kratos::Testing
