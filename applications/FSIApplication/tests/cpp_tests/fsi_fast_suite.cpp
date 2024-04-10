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

#include "fsi_fast_suite.h"

namespace Kratos::Testing {
KratosFSIFastSuite::KratosFSIFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("FSIApplication")) {
    mpFSIApp = std::make_shared<KratosFSIApplication>();
    this->mKernel.ImportApplication(mpFSIApp);
  }
}

} // namespace Kratos::Testing
