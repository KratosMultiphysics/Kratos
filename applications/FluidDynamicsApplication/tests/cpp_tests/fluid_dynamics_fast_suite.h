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

#include "fluid_dynamics_application.h"
#include "testing/testing.h"

namespace Kratos::Testing {

class KratosFluidDynamicsFastSuite : public KratosCoreFastSuite {
public:
  KratosFluidDynamicsFastSuite();

private:
  KratosFluidDynamicsApplication::Pointer mpFluidDynamicsApp;
  //  KratosLinearSolversApplication::Pointer mpLinearSolversApp;
};

} // namespace Kratos::Testing
