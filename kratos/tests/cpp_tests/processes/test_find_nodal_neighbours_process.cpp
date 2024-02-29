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

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/kratos_flags.h"
#include "utilities/cpp_tests_utilities.h"

/* Processes */
#include "processes/find_global_nodal_neighbours_for_entities_process.h"
#include "processes/find_nodal_neighbours_process.h"

namespace Kratos::Testing
{
/**
* Checks the correct work of the SPR metric process
* Test triangle
*/
KRATOS_TEST_CASE_IN_SUITE(FindNodalNeighboursProcess1, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main",2);

    auto& process_info = r_model_part.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    CppTestsUtilities::Create2DGeometry(r_model_part, "Element2D3N");

    FindNodalNeighboursProcess process(r_model_part);
    process.Execute();

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(1)->GetValue(NEIGHBOUR_NODES).size(), 3);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(1)->GetValue(NEIGHBOUR_ELEMENTS).size(), 2);

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(2)->GetValue(NEIGHBOUR_NODES).size(), 3);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(2)->GetValue(NEIGHBOUR_ELEMENTS).size(), 2);

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(3)->GetValue(NEIGHBOUR_NODES).size(), 5);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(3)->GetValue(NEIGHBOUR_ELEMENTS).size(), 4);

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(4)->GetValue(NEIGHBOUR_NODES).size(), 2);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(4)->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(5)->GetValue(NEIGHBOUR_NODES).size(), 3);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(5)->GetValue(NEIGHBOUR_ELEMENTS).size(), 2);

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(6)->GetValue(NEIGHBOUR_NODES).size(), 2);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(6)->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
}

/**
* Checks the correct work of the nodal SPR compute
* Test tetrahedra
*/
KRATOS_TEST_CASE_IN_SUITE(FindNodalNeighboursProcess2, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main",2);

    auto& process_info = r_model_part.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    CppTestsUtilities::Create3DGeometry(r_model_part, "Element3D4N");

    FindNodalNeighboursProcess process(r_model_part);
    process.Execute();

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(1)->GetValue(NEIGHBOUR_NODES).size(), 4);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(1)->GetValue(NEIGHBOUR_ELEMENTS).size(), 2);

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(2)->GetValue(NEIGHBOUR_NODES).size(), 4);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(2)->GetValue(NEIGHBOUR_ELEMENTS).size(), 2);

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(3)->GetValue(NEIGHBOUR_NODES).size(), 7);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(3)->GetValue(NEIGHBOUR_ELEMENTS).size(), 6);

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(4)->GetValue(NEIGHBOUR_NODES).size(), 5);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(4)->GetValue(NEIGHBOUR_ELEMENTS).size(), 3);

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(5)->GetValue(NEIGHBOUR_NODES).size(), 4);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(5)->GetValue(NEIGHBOUR_ELEMENTS).size(), 2);

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(6)->GetValue(NEIGHBOUR_NODES).size(), 9);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(6)->GetValue(NEIGHBOUR_ELEMENTS).size(), 9);

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(7)->GetValue(NEIGHBOUR_NODES).size(), 6);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(7)->GetValue(NEIGHBOUR_ELEMENTS).size(), 5);

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(8)->GetValue(NEIGHBOUR_NODES).size(), 8);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(8)->GetValue(NEIGHBOUR_ELEMENTS).size(), 7);

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(9)->GetValue(NEIGHBOUR_NODES).size(), 7);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(9)->GetValue(NEIGHBOUR_ELEMENTS).size(), 6);

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(10)->GetValue(NEIGHBOUR_NODES).size(), 4);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(10)->GetValue(NEIGHBOUR_ELEMENTS).size(), 2);

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(11)->GetValue(NEIGHBOUR_NODES).size(), 4);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(11)->GetValue(NEIGHBOUR_ELEMENTS).size(), 2);

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(12)->GetValue(NEIGHBOUR_NODES).size(), 4);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(12)->GetValue(NEIGHBOUR_ELEMENTS).size(), 2);
}

/**
* Checks the correct work of the nodal SPR compute
* Test Triangle
*/
KRATOS_TEST_CASE_IN_SUITE(FindNodalNeighboursProcess2_Conditions, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main",2);

    auto& process_info = r_model_part.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    CppTestsUtilities::Create2DGeometry(r_model_part, "SurfaceCondition3D3N", true, false);

    FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType> process(
        r_model_part,
        NEIGHBOUR_CONDITION_NODES);
    process.Execute();

    KRATOS_EXPECT_EQ(r_model_part.pGetNode(1)->GetValue(NEIGHBOUR_CONDITION_NODES).size(), 3);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(2)->GetValue(NEIGHBOUR_CONDITION_NODES).size(), 3);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(3)->GetValue(NEIGHBOUR_CONDITION_NODES).size(), 5);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(4)->GetValue(NEIGHBOUR_CONDITION_NODES).size(), 2);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(5)->GetValue(NEIGHBOUR_CONDITION_NODES).size(), 3);
    KRATOS_EXPECT_EQ(r_model_part.pGetNode(6)->GetValue(NEIGHBOUR_CONDITION_NODES).size(), 2);
}

}  // namespace Kratos::Testing.
