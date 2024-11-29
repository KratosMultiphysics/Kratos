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

#include "testing/testing.h"
#include "convection_diffusion_application.h"

namespace Kratos::Testing {

class KratosConvectionDiffusionFastSuite : public KratosCoreFastSuite {
public:
  KratosConvectionDiffusionFastSuite();

private:
  KratosConvectionDiffusionApplication::Pointer mpConvectionDiffusionApp;
  //  KratosLinearSolversApplication::Pointer mpLinearSolversApp;
};

} // namespace Kratos::Testing
