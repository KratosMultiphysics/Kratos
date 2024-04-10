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

#include "exaqute_sandbox_fast_suite.h"

namespace Kratos::Testing {
KratosExaquteSanboxFastSuite::KratosExaquteSanboxFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("ExaquteSanboxApplication")) {
    mpStructuralApp = std::make_shared<KratosExaquteSanboxApplication>();
    this->mKernel.ImportApplication(mpExaquteSanboxApp);
  }
}

} // namespace Kratos::Testing
