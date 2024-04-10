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

#include "constitutive_laws_fast_suite.h"

namespace Kratos::Testing {
KratosConstitutiveLawsFastSuite::KratosConstitutiveLawsFastSuite()
    : KratosCoreFastSuite() {
  if (!this->mKernel.IsImported("ConstitutiveLawsApplication")) {
    mpConstitutiveLawsApp = std::make_shared<KratosConstitutiveLawsApplication>();
    this->mKernel.ImportApplication(mpConstitutiveLawsApp);
  }
}

} // namespace Kratos::Testing
