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

#include "iga_fast_suite.h"

namespace Kratos::Testing {
KratosIgaFastSuite::KratosIgaFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("IgaApplication")) {
    mpIgaApp = std::make_shared<KratosIgaApplication>();
    this->mKernel.ImportApplication(mpIgaApp);
  }
}

} // namespace Kratos::Testing
