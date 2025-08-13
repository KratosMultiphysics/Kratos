//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//                   Jordi Cotela Dalmau
//
//

// System includes

// External includes

// Project includes
#include "tests/test_utilities/test_suite.h"

namespace Kratos::Testing 
{
    
void KratosCoreFastSuite::SetUp() { 
    // TODO: Control the log level of the tests
    // mCoutBuffer = std::cout.rdbuf();
    // mCerrBuffer = std::cerr.rdbuf();

    // std::cout.rdbuf(mStream.rdbuf());
    // std::cerr.rdbuf(mStream.rdbuf());   
}

void KratosCoreFastSuite::TearDown() { 
    // TODO: Control the log level of the tests
    // std::cout.rdbuf(mCoutBuffer);
    // std::cerr.rdbuf(mCerrBuffer); 
}

} // namespace Kratos::Testing