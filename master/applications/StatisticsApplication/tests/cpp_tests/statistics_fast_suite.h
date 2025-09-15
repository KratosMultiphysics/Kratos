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

#pragma once

#include "structural_mechanics_application.h"
#include "testing/testing.h"

namespace Kratos::Testing {

class KratosStructuralMechanicsFastSuite : public KratosCoreFastSuite {
public:
  KratosStructuralMechanicsFastSuite();

private:
  KratosStructuralMechanicsApplication::Pointer mpStructuralApp;
  //  KratosLinearSolversApplication::Pointer mpLinearSolversApp;
};

} // namespace Kratos::Testing
