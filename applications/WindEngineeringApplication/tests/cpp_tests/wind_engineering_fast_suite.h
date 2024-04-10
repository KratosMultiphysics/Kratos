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

#include "wind_engineering_application.h"
#include "testing/testing.h"

namespace Kratos::Testing {

class KratosWindEngineeringFastSuite : public KratosCoreFastSuite {
public:
  KratosWindEngineeringFastSuite();

private:
  KratosWindEngineeringApplication::Pointer mpWindEngineeringApp;
  //  KratosLinearSolversApplication::Pointer mpLinearSolversApp;
};

} // namespace Kratos::Testing
