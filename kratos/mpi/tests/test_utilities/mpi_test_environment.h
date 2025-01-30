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
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "tests/test_utilities/test_environment.h"

namespace Kratos::Testing
{

/*
 * This Fixture creates a new kernel instance for kratos, so the test is able to interact with the database.
 * Its called this way to that all tests belong to a existing kernel fixture
*/
class KRATOS_API(KRATOS_MPI_CORE) KratosMpiTestEnv : public ::testing::Environment 
{
    public:
        ~KratosMpiTestEnv() override {}

        void SetUp() override;
        void TearDown() override;
};

} // namespace Kratos::Testing