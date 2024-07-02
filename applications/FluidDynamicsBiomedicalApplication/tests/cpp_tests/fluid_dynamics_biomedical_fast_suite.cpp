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
//                   Carlos A. Roig
//

#pragma once

#include "testing/testing.h"
#include "fluid_dynamics_Biomedical_application.h"

namespace Kratos::Testing {

class FluidDynamicsBiomedicalApplicationFastSuite : public KratosCoreFastSuite {
public:
  FluidDynamicsBiomedicalApplicationFastSuite();

private:
  KratosFluidDynamicsBiomedicalApplication::Pointer mpFluidDynamicsBiomedicalApp;
};

} // namespace Kratos::Testing
