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

// External includes
#include <gmock/gmock.h>
#include <gtest/gtest.h>

// Project includes
#include "structural_mechanics_fast_suite.h"

namespace Kratos::Testing {

KratosStructuralMechanicsFastSuite::KratosStructuralMechanicsFastSuite()
    : KratosCoreFastSuite() {
  mpStructuralMechanicsApp = std::make_shared<KratosStructuralMechanicsApplication>();
  this->ImportApplicationIntoKernel(mpStructuralMechanicsApp);
}

} // namespace Kratos::Testing