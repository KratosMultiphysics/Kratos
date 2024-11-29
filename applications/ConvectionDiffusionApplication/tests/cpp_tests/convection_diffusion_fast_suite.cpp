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

// Project includes
#include "convection_diffusion_fast_suite.h"

namespace Kratos::Testing 
{

KratosConvectionDiffusionFastSuite::KratosConvectionDiffusionFastSuite()
    : KratosCoreFastSuite() 
{
    mpConvectionDiffusionApp = std::make_shared<KratosConvectionDiffusionApplication>();
    this->ImportApplicationIntoKernel(mpConvectionDiffusionApp);
}

} // namespace Kratos::Testing