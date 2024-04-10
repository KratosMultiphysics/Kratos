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

#include "contact_structural_mechanics_fast_suite.h"

namespace Kratos::Testing {
KratosContactStructuralMechanicsFastSuite::KratosContactStructuralMechanicsFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("ContactStructuralMechanicsApplication")) {
    mpContactStructuralMechanicsApp = std::make_shared<KratosContactStructuralMechanicsApplication>();
    this->mKernel.ImportApplication(mpContactStructuralMechanicsApp);
  }
}

} // namespace Kratos::Testing
