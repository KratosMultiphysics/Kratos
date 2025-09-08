//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

#pragma once

#define KRATOS_SKIP_TEST GTEST_SKIP()
#define KRATOS_SKIP_TEST_IF(conditional) if (conditional) GTEST_SKIP()
#define KRATOS_SKIP_TEST_IF_NOT(conditional) if (!(conditional)) GTEST_SKIP()
