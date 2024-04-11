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

#pragma once

#include "structural_mechanics_application.h"
#include "testing/testing.h"

namespace Kratos::Testing {

class KratosStructuralMechanicsFastSuite : public KratosCoreFastSuite {
public:
  KratosStructuralMechanicsFastSuite();

private:
  KratosStructuralMechanicsApplication::Pointer mpStructuralMechanicsApp;
};

} // namespace Kratos::Testing
