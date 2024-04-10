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

#include "chsarp_wrapper_fast_suite.h"

namespace Kratos::Testing {
KratosCSharpWrapperFastSuite::KratosCSharpWrapperFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("CSharpWrapperApplication")) {
    mpStructuralApp = std::make_shared<KratosCSharpWrapperApplication>();
    this->mKernel.ImportApplication(mpStructuralApp);
  }
}

} // namespace Kratos::Testing
