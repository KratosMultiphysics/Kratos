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

#pragma once

#include "dam_application.h"
#include "structural_mechanics_application.h"
#include "testing/testing.h"

namespace Kratos::Testing {

class KratosDamFastSuite : public KratosCoreFastSuite {
public:
  KratosDamFastSuite();

private:
  KratosStructuralMechanicsApplication::Pointer mpStructuralMechanicsApp;
  KratosDamApplication::Pointer mpDamApp;
};

} // namespace Kratos::Testing
