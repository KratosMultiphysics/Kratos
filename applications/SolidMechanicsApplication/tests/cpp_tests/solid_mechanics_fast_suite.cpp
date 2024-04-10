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

#include "solid_mechanics_fast_suite.h"

namespace Kratos::Testing {
KratosSolidMechanicsFastSuite::KratosSolidMechanicsFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("SolidMechanicsApplication")) {
    mpSolidApp = std::make_shared<KratosSolidMechanicsApplication>();
    this->mKernel.ImportApplication(mpSolidApp);
  }
}

} // namespace Kratos::Testing
