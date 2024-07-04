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

// Project includes
#include "tests/cpp_tests/meshing_fast_suite.h"

namespace Kratos::Testing 
{

KratosMeshingApplicationFastSuite::KratosMeshingApplicationFastSuite()
    : KratosCoreFastSuite() 
{
    mpMeshingApp = std::make_shared<KratosMeshingApplication>();
    this->ImportApplicationIntoKernel(mpMeshingApp);
}

} // namespace Kratos::Testing