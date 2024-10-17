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
#include "tests/cpp_tests/mapping_fast_suite.h"

namespace Kratos::Testing 
{

KratosMappingApplicationSerialTestSuite::KratosMappingApplicationSerialTestSuite()
    : KratosCoreFastSuite() 
{
    mpMappingApp = std::make_shared<KratosMappingApplication>();
    this->ImportApplicationIntoKernel(mpMappingApp);
}

} // namespace Kratos::Testing