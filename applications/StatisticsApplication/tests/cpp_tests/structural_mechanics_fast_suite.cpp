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

#include "structural_mechanics_fast_suite.h"

namespace Kratos::Testing {
KratosStructuralMechanicsFastSuite::KratosStructuralMechanicsFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("StructuralMechanicsApplication")) {
    mpStructuralApp = std::make_shared<KratosStructuralMechanicsApplication>();
    this->mKernel.ImportApplication(mpStructuralApp);
  }
}

} // namespace Kratos::Testing
