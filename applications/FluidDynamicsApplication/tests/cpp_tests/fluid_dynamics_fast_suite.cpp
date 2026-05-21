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
#include "fluid_dynamics_fast_suite.h"

namespace Kratos::Testing 
{

FluidDynamicsApplicationFastSuite::FluidDynamicsApplicationFastSuite()
    : KratosCoreFastSuite() 
{
    mpFluidDynamicsApp = std::make_shared<KratosFluidDynamicsApplication>();
    this->ImportApplicationIntoKernel(mpFluidDynamicsApp);
}

} // namespace Kratos::Testing