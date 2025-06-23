//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//
//

// System includes
#include <numeric>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/expect.h"
#include "includes/registry.h"
#include "processes/process.h"
#include "processes/output_process.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(RegistryItemGetValue, KratosCoreFastSuite)
{
    auto pProcess = Registry::GetValue<Process>("Processes.KratosMultiphysics.Process.Prototype");

    KRATOS_EXPECT_TRUE((std::is_base_of<Process,decltype(pProcess)>::value))
}

KRATOS_TEST_CASE_IN_SUITE(RegistryItemGetValueDerivedAsBase, KratosCoreFastSuite)
{
    auto pProcess = Registry::GetValue<Process>("Processes.KratosMultiphysics.OutputProcess.Prototype");

    KRATOS_EXPECT_TRUE((std::is_same<Process,decltype(pProcess)>::value))

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(Registry::GetValue<Process>("Processes.KratosMultiphysics.output_process.Prototype"), "");
}


KRATOS_TEST_CASE_IN_SUITE(RegistryItemGetValueDerivedasDerived, KratosCoreFastSuite)
{
    auto pProcess = Registry::GetValueAs<Process, OutputProcess>("Processes.KratosMultiphysics.OutputProcess.Prototype");

    KRATOS_EXPECT_TRUE((std::is_same<OutputProcess,decltype(pProcess)>::value))
}

}  // namespace Kratos::Testing.
