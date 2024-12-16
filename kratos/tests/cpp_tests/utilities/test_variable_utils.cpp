//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                  Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "utilities/variable_utils.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(VariableUtilsSetHistoricalVariablesToZero, KratosCoreFastSuite)
{
    // Set auxilary nodal structure
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    r_test_model_part.AddNodalSolutionStepVariable(PRESSURE);
    r_test_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_test_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    r_test_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_test_model_part.AddNodalSolutionStepVariable(DEFORMATION_GRADIENT);
    for (std::size_t i = 0; i < 10; ++i) {
        r_test_model_part.CreateNewNode(i,0.0,0.0,0.0);
    }

    // Set fake values in the historical database
    array_1d<double,3> aux_vect;
    Matrix aux_mat = ZeroMatrix(3,3);
    for (std::size_t i = 0; i < 3; ++i) {
        aux_vect[i] = 1.0;
        for (std::size_t j = 0; j < 3; ++j) {
            aux_mat(i,j) = 1.0;
        }
    }
    for (auto& r_node : r_test_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(PRESSURE) = 1.0;
        r_node.FastGetSolutionStepValue(TEMPERATURE) = 1.0;
        noalias(r_node.FastGetSolutionStepValue(VELOCITY)) = aux_vect;
        noalias(r_node.FastGetSolutionStepValue(DISPLACEMENT)) = aux_vect;
        r_node.FastGetSolutionStepValue(DEFORMATION_GRADIENT) = aux_mat;
    }

    // Set some values to zero in the non-historical database
    VariableUtils::SetHistoricalVariablesToZero(r_test_model_part.Nodes(), PRESSURE, TEMPERATURE, VELOCITY, DISPLACEMENT, DEFORMATION_GRADIENT);

    // Values are properly allocated
    const double tolerance = 1.0e-12;
    const auto ZeroArray = array_1d<double, 3>(3, 0.0);
    for (const auto& r_node : r_test_model_part.Nodes()) {
        KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(PRESSURE), 0.0, tolerance);
        KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(TEMPERATURE), 0.0, tolerance);
        KRATOS_EXPECT_VECTOR_NEAR(r_node.FastGetSolutionStepValue(VELOCITY), ZeroArray, tolerance);
        KRATOS_EXPECT_VECTOR_NEAR(r_node.FastGetSolutionStepValue(DISPLACEMENT), ZeroArray, tolerance);
        KRATOS_EXPECT_MATRIX_NEAR(r_node.FastGetSolutionStepValue(DEFORMATION_GRADIENT), ZeroMatrix(0,0), tolerance);
    }
}

KRATOS_TEST_CASE_IN_SUITE(VariableUtilsSetNonHistoricalVariablesToZeroNodes, KratosCoreFastSuite)
{
    // Set auxilary nodal structure
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    for (std::size_t i = 0; i < 10; ++i) {
        r_test_model_part.CreateNewNode(i,0.0,0.0,0.0);
    }

    // Set some values to zero in the non-historical database
    VariableUtils::SetNonHistoricalVariablesToZero(r_test_model_part.Nodes(), PRESSURE, TEMPERATURE, VELOCITY, DISPLACEMENT, DEFORMATION_GRADIENT);

    // Values are properly allocated
    const double tolerance = 1.0e-12;
    const auto ZeroArray = array_1d<double, 3>(3, 0.0);
    for (const auto& r_node : r_test_model_part.Nodes()) {
        KRATOS_EXPECT_NEAR(r_node.GetValue(PRESSURE), 0.0, tolerance);
        KRATOS_EXPECT_NEAR(r_node.GetValue(TEMPERATURE), 0.0, tolerance);
        KRATOS_EXPECT_VECTOR_NEAR(r_node.GetValue(VELOCITY), ZeroArray, tolerance);
        KRATOS_EXPECT_VECTOR_NEAR(r_node.GetValue(DISPLACEMENT), ZeroArray, tolerance);
        KRATOS_EXPECT_MATRIX_NEAR(r_node.GetValue(DEFORMATION_GRADIENT), ZeroMatrix(0,0), tolerance);
    }
}

KRATOS_TEST_CASE_IN_SUITE(VariableUtilsSetNonHistoricalVariablesToZeroElements, KratosCoreFastSuite)
{
    // Set auxilary elemental structure
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    for (std::size_t i = 0; i < 6; ++i) {
        r_test_model_part.CreateNewNode(i,0.0,0.0,0.0);
    }
    auto p_prop = Kratos::make_shared<Properties>(0);
    r_test_model_part.CreateNewElement("Element2D2N",1,{{1,2}},p_prop);
    r_test_model_part.CreateNewElement("Element2D2N",2,{{2,3}},p_prop);
    r_test_model_part.CreateNewElement("Element2D2N",3,{{3,4}},p_prop);
    r_test_model_part.CreateNewElement("Element2D2N",4,{{4,5}},p_prop);

    // Set some values to zero in the non-historical database
    VariableUtils::SetNonHistoricalVariablesToZero(r_test_model_part.Elements(), PRESSURE, TEMPERATURE, VELOCITY, DISPLACEMENT, DEFORMATION_GRADIENT);

    // Values are properly allocated
    const double tolerance = 1.0e-12;
    const auto ZeroArray = array_1d<double, 3>(3, 0.0);
    for (const auto& r_element : r_test_model_part.Elements()) {
        KRATOS_EXPECT_NEAR(r_element.GetValue(PRESSURE), 0.0, tolerance);
        KRATOS_EXPECT_NEAR(r_element.GetValue(TEMPERATURE), 0.0, tolerance);
        KRATOS_EXPECT_VECTOR_NEAR(r_element.GetValue(VELOCITY), ZeroArray, tolerance);
        KRATOS_EXPECT_VECTOR_NEAR(r_element.GetValue(DISPLACEMENT), ZeroArray, tolerance);
        KRATOS_EXPECT_MATRIX_NEAR(r_element.GetValue(DEFORMATION_GRADIENT), ZeroMatrix(0,0), tolerance);
    }
}

KRATOS_TEST_CASE_IN_SUITE(VariableUtilsSetNonHistoricalVariablesToZeroConditions, KratosCoreFastSuite)
{
    // Set auxilary elemental structure
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    for (std::size_t i = 0; i < 6; ++i) {
        r_test_model_part.CreateNewNode(i,0.0,0.0,0.0);
    }
    auto p_prop = Kratos::make_shared<Properties>(0);
    r_test_model_part.CreateNewCondition("LineCondition2D2N",1,{{1,2}},p_prop);
    r_test_model_part.CreateNewCondition("LineCondition2D2N",2,{{2,3}},p_prop);
    r_test_model_part.CreateNewCondition("LineCondition2D2N",3,{{3,4}},p_prop);
    r_test_model_part.CreateNewCondition("LineCondition2D2N",4,{{4,5}},p_prop);

    // Set some values to zero in the non-historical database
    VariableUtils::SetNonHistoricalVariablesToZero(r_test_model_part.Conditions(), PRESSURE, TEMPERATURE, VELOCITY, DISPLACEMENT, DEFORMATION_GRADIENT);

    // Values are properly allocated
    const double tolerance = 1.0e-12;
    const auto ZeroArray = array_1d<double, 3>(3, 0.0);
    for (const auto& r_condition : r_test_model_part.Conditions()) {
        KRATOS_EXPECT_NEAR(r_condition.GetValue(PRESSURE), 0.0, tolerance);
        KRATOS_EXPECT_NEAR(r_condition.GetValue(TEMPERATURE), 0.0, tolerance);
        KRATOS_EXPECT_VECTOR_NEAR(r_condition.GetValue(VELOCITY), ZeroArray, tolerance);
        KRATOS_EXPECT_VECTOR_NEAR(r_condition.GetValue(DISPLACEMENT), ZeroArray, tolerance);
        KRATOS_EXPECT_MATRIX_NEAR(r_condition.GetValue(DEFORMATION_GRADIENT), ZeroMatrix(0,0), tolerance);
    }
}

KRATOS_TEST_CASE_IN_SUITE(VariableUtilsEraseNonHistoricalVariableConditions, KratosCoreFastSuite)
{
    // Set auxilary elemental structure
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    for (std::size_t i = 0; i < 6; ++i) {
        r_test_model_part.CreateNewNode(i,0.0,0.0,0.0);
    }
    auto p_prop = Kratos::make_shared<Properties>(0);
    r_test_model_part.CreateNewCondition("LineCondition2D2N",1,{{1,2}},p_prop);
    r_test_model_part.CreateNewCondition("LineCondition2D2N",2,{{2,3}},p_prop);
    r_test_model_part.CreateNewCondition("LineCondition2D2N",3,{{3,4}},p_prop);
    r_test_model_part.CreateNewCondition("LineCondition2D2N",4,{{4,5}},p_prop);

    // Set some values to zero in the non-historical database
    VariableUtils::SetNonHistoricalVariablesToZero(r_test_model_part.Conditions(), PRESSURE, TEMPERATURE);

    // Erase some values in the non-historical database
    VariableUtils().EraseNonHistoricalVariable(TEMPERATURE, r_test_model_part.Conditions());

    // Values are properly allocated
    for (const auto& r_condition : r_test_model_part.Conditions()) {
        KRATOS_EXPECT_TRUE(r_condition.Has(PRESSURE));
        KRATOS_EXPECT_FALSE(r_condition.Has(TEMPERATURE));
    }
}

KRATOS_TEST_CASE_IN_SUITE(VariableUtilsClearNonHistoricalDataConditions, KratosCoreFastSuite)
{
    // Set auxilary elemental structure
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    for (std::size_t i = 0; i < 6; ++i) {
        r_test_model_part.CreateNewNode(i,0.0,0.0,0.0);
    }
    auto p_prop = Kratos::make_shared<Properties>(0);
    r_test_model_part.CreateNewCondition("LineCondition2D2N",1,{{1,2}},p_prop);
    r_test_model_part.CreateNewCondition("LineCondition2D2N",2,{{2,3}},p_prop);
    r_test_model_part.CreateNewCondition("LineCondition2D2N",3,{{3,4}},p_prop);
    r_test_model_part.CreateNewCondition("LineCondition2D2N",4,{{4,5}},p_prop);

    // Set some values to zero in the non-historical database
    VariableUtils::SetNonHistoricalVariablesToZero(r_test_model_part.Conditions(), PRESSURE, TEMPERATURE);

    // Erase some values in the non-historical database
    VariableUtils().ClearNonHistoricalData(r_test_model_part.Conditions());

    // Values are properly allocated
    for (const auto& r_condition : r_test_model_part.Conditions()) {
        KRATOS_EXPECT_FALSE(r_condition.Has(PRESSURE));
        KRATOS_EXPECT_FALSE(r_condition.Has(TEMPERATURE));
    }
}

KRATOS_TEST_CASE_IN_SUITE(VariableUtilsGetPositionsVector, KratosCoreFastSuite)
{
    // Set auxilary elemental structure
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    for (std::size_t i = 0; i < 4; ++i) {
        r_test_model_part.CreateNewNode(i, 1,2,3);
    }

    const Vector init_pos_1 = VariableUtils().GetInitialPositionsVector(std::as_const(r_test_model_part.Nodes()), 3);
    const std::vector<double> init_pos_2 = VariableUtils().GetInitialPositionsVector<std::vector<double>>(std::as_const(r_test_model_part.Nodes()), 3);

    const Vector curr_pos_1 = VariableUtils().GetCurrentPositionsVector(std::as_const(r_test_model_part.Nodes()), 3);
    const std::vector<double> curr_pos_2 = VariableUtils().GetCurrentPositionsVector<std::vector<double>>(std::as_const(r_test_model_part.Nodes()), 3);

    const std::vector<double> exp_pos{1,2,3,1,2,3,1,2,3,1,2,3};

    KRATOS_EXPECT_VECTOR_EQ(exp_pos, init_pos_1);
    KRATOS_EXPECT_VECTOR_EQ(exp_pos, init_pos_2);

    KRATOS_EXPECT_VECTOR_EQ(exp_pos, curr_pos_1);
    KRATOS_EXPECT_VECTOR_EQ(exp_pos, curr_pos_2);

}

} // namespace Kratos::Testing
