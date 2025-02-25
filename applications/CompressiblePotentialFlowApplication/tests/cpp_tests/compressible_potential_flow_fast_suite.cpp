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
#include "compressible_potential_flow_fast_suite.h"

namespace Kratos::Testing 
{

CompressiblePotentialFlowApplicationFastSuite::CompressiblePotentialFlowApplicationFastSuite()
    : KratosCoreFastSuite() 
{
    mpCompressiblePotentialFlowApp = std::make_shared<KratosCompressiblePotentialFlowApplication>();
    this->ImportApplicationIntoKernel(mpCompressiblePotentialFlowApp);
}

} // namespace Kratos::Testing