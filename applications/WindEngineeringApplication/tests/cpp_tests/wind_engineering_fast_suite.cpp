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

#include "wind_engineering_fast_suite.h"

namespace Kratos::Testing {
KratosWindEngineeringsFastSuite::KratosWindEngineeringFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("WindEngineeringApplication")) {
    mpWindEngineeringApp = std::make_shared<KratosWindEngineeringsApplication>();
    this->mKernel.ImportApplication(mpWindEngineeringApp);
  }
}

} // namespace Kratos::Testing
