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
#include <gtest/gtest.h>

// Project includes

class KRATOS_API(KRATOS_TEST_UTILS) KratosTestEnv : public ::testing::Environment
{
    public:
        KratosTestEnv();
        ~KratosTestEnv() override {}
        void SetUp() override;
        void TearDown() override;
};
