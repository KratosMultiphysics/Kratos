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

#include "co_simulation_application.h"
#include "testing/testing.h"

namespace Kratos::Testing {

class KratosCoSimulationFastSuite : public KratosCoreFastSuite {
public:
  KratosCoSimulationFastSuite();

private:
  KratosCoSimulationApplication::Pointer mpCoSimulationApp;
  //  KratosLinearSolversApplication::Pointer mpLinearSolversApp;
};

} // namespace Kratos::Testing
