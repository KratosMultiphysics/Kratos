//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "meshioplusplus_application.h"
#include "testing/testing.h"

namespace Kratos::Testing
{

/**
 * @brief Test fixture for the MeshioPlusPlusApplication.
 * @details The application itself is imported by the custom main
 * (meshioplusplus_fast_suite.cpp) through the shared application initializer list, so this
 * only names the suite the test cases belong to.
 */
class KratosMeshioPlusPlusFastSuite : public KratosCoreFastSuite
{
};

} // namespace Kratos::Testing
