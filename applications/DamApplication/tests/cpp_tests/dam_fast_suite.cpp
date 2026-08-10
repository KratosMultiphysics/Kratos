// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    DamApplication developers
//

// External includes
#include <gmock/gmock.h>
#include <gtest/gtest.h>

// Project includes
#include "dam_fast_suite.h"

namespace Kratos::Testing {

KratosDamFastSuite::KratosDamFastSuite()
    : KratosCoreFastSuite() {
  mpStructuralMechanicsApp = std::make_shared<KratosStructuralMechanicsApplication>();
  this->ImportApplicationIntoKernel(mpStructuralMechanicsApp);
  mpDamApp = std::make_shared<KratosDamApplication>();
  this->ImportApplicationIntoKernel(mpDamApp);
}

} // namespace Kratos::Testing
