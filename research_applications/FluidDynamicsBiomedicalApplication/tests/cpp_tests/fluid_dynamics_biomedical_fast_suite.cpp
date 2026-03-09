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
#include "tests/cpp_tests/fluid_dynamics_biomedical_fast_suite.h"

namespace Kratos::Testing 
{

FluidDynamicsBiomedicalApplicationFastSuite::FluidDynamicsBiomedicalApplicationFastSuite()
    : KratosCoreFastSuite() 
{
    mpFluidDynamicsBiomedicalApp = std::make_shared<KratosFluidDynamicsBiomedicalApplication>();
    this->ImportApplicationIntoKernel(mpFluidDynamicsBiomedicalApp);
}

} // namespace Kratos::Testing