// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
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
