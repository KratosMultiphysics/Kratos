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

#include "rans_fast_suite.h"

namespace Kratos::Testing {
KratosRANSFastSuite::KratosRANSFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("RANSApplication")) {
    mpRANSApp = std::make_shared<KratosRANSApplication>();
    this->mKernel.ImportApplication(mpRANSApp);
  }
}

} // namespace Kratos::Testing
