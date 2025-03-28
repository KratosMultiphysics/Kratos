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
#include "tests/cpp_tests/mesh_moving_fast_suite.h"

namespace Kratos::Testing 
{

MeshMovingApplicationFastSuite::MeshMovingApplicationFastSuite()
    : KratosCoreFastSuite() 
{
    mpMeshMovingApp = std::make_shared<KratosMeshMovingApplication>();
    this->ImportApplicationIntoKernel(mpMeshMovingApp);
}

} // namespace Kratos::Testing