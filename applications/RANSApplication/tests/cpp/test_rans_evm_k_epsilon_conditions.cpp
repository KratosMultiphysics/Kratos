//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//

// System includes
#include <functional>
#include <iomanip>
#include <string>
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/kratos_components.h"
#include "includes/variables.h"
#include "testing/testing.h"

// Application includes
#include "custom_conditions/evm_k_epsilon/rans_evm_k_epsilon_epsilon_wall.h"
#include "custom_conditions/evm_k_epsilon/rans_evm_k_epsilon_vms_monolithic_wall.h"
#include "rans_application_variables.h"

namespace Kratos
{
namespace Testing
{
namespace
{
void CreateRansEvmKEpsilonUnitTestModelPart(const std::string& rConditionName,
                                    const std::vector<std::string>& rDofVariableNamesList,
                                    ModelPart& rModelPart)
{
    const auto& r_proto = KratosComponents<Condition>::Get(rConditionName);
    auto node_ids = std::vector<ModelPart::IndexType>{};
    Matrix coords;
    r_proto.GetGeometry().PointsLocalCoordinates(coords);
    if (coords.size2() == 1)
    {
        for (std::size_t i = 0; i < coords.size1(); ++i)
            rModelPart.CreateNewNode(i + 1, coords(i, 0), 0.0, 0.0);
    }
    else
    {
        for (std::size_t i = 0; i < coords.size1(); ++i)
            rModelPart.CreateNewNode(i + 1, coords(i, 0), coords(i, 1), 0.0);
    }
    for (auto& r_node : rModelPart.Nodes())
    {
        for (std::string dof_variable_name : rDofVariableNamesList)
        {
            if (KratosComponents<Variable<double>>::Has(dof_variable_name))
            {
                r_node
                    .AddDof(KratosComponents<Variable<double>>::Get(dof_variable_name))
                    .SetEquationId(r_node.Id());
            }
            else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(dof_variable_name))
            {
                const auto& r_variable_x =
                    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Get(
                        dof_variable_name + "_X");
                const auto& r_variable_y =
                    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Get(
                        dof_variable_name + "_Y");
                r_node.AddDof(r_variable_x).SetEquationId(r_node.Id() * 5);
                r_node.AddDof(r_variable_y).SetEquationId(r_node.Id() * 5);

                if (coords.size2() != 2)
                {
                    const auto& r_variable_z =
                        KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Get(
                            dof_variable_name + "_Z");
                    r_node.AddDof(r_variable_z).SetEquationId(r_node.Id() * 10);
                }
            }
            else
            {
                KRATOS_ERROR
                    << "DOF variable name = " << dof_variable_name
                    << " not found in double, 3d-double variable lists.\n";
            }
        }
        node_ids.push_back(r_node.Id());
    }
    auto p_prop = rModelPart.CreateNewProperties(1);
    auto p_condition = rModelPart.CreateNewCondition(rConditionName, 1, node_ids, p_prop);
    rModelPart.SetBufferSize(1);
    p_condition->Check(rModelPart.GetProcessInfo());
}

void RansEvmKEpsilonEpsilonWall2D2N_SetUp(ModelPart& rModelPart)
{
    // rModelPart.AddNodalSolutionStepVariable(DISTANCE);
    rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE);
    rModelPart.AddNodalSolutionStepVariable(RANS_Y_PLUS);
    CreateRansEvmKEpsilonUnitTestModelPart("RansEvmKEpsilonEpsilonWall2D2N",
                                   {"TURBULENT_ENERGY_DISSIPATION_RATE"}, rModelPart);
}

void RansEvmKEpsilonVmsMonolithicWall2D2N_SetUp(ModelPart& rModelPart)
{
    rModelPart.AddNodalSolutionStepVariable(DISTANCE);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
    rModelPart.AddNodalSolutionStepVariable(RANS_Y_PLUS);
    rModelPart.AddNodalSolutionStepVariable(DENSITY);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE);
    rModelPart.AddNodalSolutionStepVariable(VISCOSITY);
    CreateRansEvmKEpsilonUnitTestModelPart("RansEvmKEpsilonVmsMonolithicWall2D2N",
                                   {"VELOCITY", "PRESSURE"}, rModelPart);
}

void RansEvmKEpsilonEpsilonWall2D2N_AssignTestData(ModelPart& rModelPart)
{
    rModelPart.GetProcessInfo()[TURBULENCE_RANS_C_MU] = 0.09;
    rModelPart.GetProcessInfo()[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA] = 0.98;
    auto& node1 = rModelPart.GetNode(1);
    node1.FastGetSolutionStepValue(KINEMATIC_VISCOSITY) = 0.24;
    node1.FastGetSolutionStepValue(RANS_Y_PLUS) = 11.24;
    node1.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 25.90;
    node1.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) = 15.90;
    node1.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = 1.29;
    auto& node2 = rModelPart.GetNode(2);
    node2.FastGetSolutionStepValue(KINEMATIC_VISCOSITY) = 0.28;
    node2.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 5.00;
    node2.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = 1.04;
    node2.FastGetSolutionStepValue(RANS_Y_PLUS) = 11.24;
    node2.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) = 11.90;
}

void RansEvmKEpsilonVmsMonolithicWall2D2N_AssignTestData(ModelPart& rModelPart)
{
    rModelPart.GetProcessInfo()[TURBULENCE_RANS_C_MU] = 0.09;
    rModelPart.GetProcessInfo()[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA] = 1.3;
    auto& node1 = rModelPart.GetNode(1);
    node1.FastGetSolutionStepValue(KINEMATIC_VISCOSITY) = 0.24;
    node1.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) = 12.90;
    node1.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 25.90;
    node1.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = 1.29;
    node1.FastGetSolutionStepValue(RANS_Y_PLUS) = 8.34;
    node1.FastGetSolutionStepValue(VELOCITY_X) = -1.85;
    node1.FastGetSolutionStepValue(VELOCITY_Y) = -2.15;
    node1.FastGetSolutionStepValue(PRESSURE) = -5.15;
    node1.FastGetSolutionStepValue(DISTANCE) = 2.15;
    node1.FastGetSolutionStepValue(VISCOSITY) = 5e-2;
    node1.FastGetSolutionStepValue(DENSITY) = 2.3;
    auto& node2 = rModelPart.GetNode(2);
    node2.FastGetSolutionStepValue(KINEMATIC_VISCOSITY) = 0.28;
    node2.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) = 5.00;
    node2.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 5.00;
    node2.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = 1.04;
    node2.FastGetSolutionStepValue(RANS_Y_PLUS) = 25.2;
    node2.FastGetSolutionStepValue(VELOCITY_X) = -0.78;
    node2.FastGetSolutionStepValue(VELOCITY_Y) = -0.06;
    node2.FastGetSolutionStepValue(PRESSURE) = -4.15;
    node2.FastGetSolutionStepValue(DISTANCE) = 8.15;
    node2.FastGetSolutionStepValue(VISCOSITY) = 5e-2;
    node1.FastGetSolutionStepValue(DENSITY) = 4.3;
}

void RansEvmKEpsilonVmsMonolithicWall2D2N_EvaluateTest(
    std::function<void(Condition&, const ProcessInfo&, const bool, const bool)> TestEvaluationMethod)
{
    // Setup:
    Model model;
    auto& model_part = model.CreateModelPart("test");
    RansEvmKEpsilonVmsMonolithicWall2D2N_SetUp(model_part);
    RansEvmKEpsilonVmsMonolithicWall2D2N_AssignTestData(model_part);
    auto& r_condition = model_part.Conditions().front();
    auto& r_process_info = model_part.GetProcessInfo();

    auto permutated_evaluation_method =
        [TestEvaluationMethod](Condition& rCondition, ProcessInfo& rProcessInfo,
                               const bool IsSlip, const bool IsCoSolvingProcessActive) {
            rCondition.Set(SLIP, IsSlip);
            rCondition.GetGeometry()[0].Set(SLIP, IsSlip);
            rCondition.GetGeometry()[1].Set(SLIP, IsSlip);
            rProcessInfo[IS_CO_SOLVING_PROCESS_ACTIVE] = IsCoSolvingProcessActive;
            TestEvaluationMethod(rCondition, rProcessInfo, IsSlip, IsCoSolvingProcessActive);
        };

    permutated_evaluation_method(r_condition, r_process_info, false, false);
    permutated_evaluation_method(r_condition, r_process_info, false, true);
    permutated_evaluation_method(r_condition, r_process_info, true, false);
    permutated_evaluation_method(r_condition, r_process_info, true, true);
}

void RansEvmKEpsilonEpsilonWall2D2N_EvaluateTest(
    std::function<void(Condition&, ProcessInfo&, const bool)> TestEvaluationMethod)
{
    // Setup:
    Model model;
    auto& model_part = model.CreateModelPart("test");
    RansEvmKEpsilonEpsilonWall2D2N_SetUp(model_part);
    RansEvmKEpsilonEpsilonWall2D2N_AssignTestData(model_part);
    auto& r_condition = model_part.Conditions().front();
    auto& r_process_info = model_part.GetProcessInfo();

    auto permutated_evaluation_method =
        [TestEvaluationMethod](Condition& rCondition, ProcessInfo& rProcessInfo,
                               const bool IsSlip) {
            rCondition.Set(SLIP, IsSlip);
            rCondition.GetGeometry()[0].Set(SLIP, IsSlip);
            rCondition.GetGeometry()[1].Set(SLIP, IsSlip);
            TestEvaluationMethod(rCondition, rProcessInfo, IsSlip);
        };

    permutated_evaluation_method(r_condition, r_process_info, false);
    permutated_evaluation_method(r_condition, r_process_info, true);
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonWall2D2N_EquationIdVector, KratosRansFastSuite)
{
    auto evaluation_method = [](Condition& rCondition,
                                ProcessInfo& rProcessInfo, const bool IsSlip) {
        auto eqn_ids = std::vector<std::size_t>{};
        rCondition.EquationIdVector(eqn_ids, rProcessInfo);
        KRATOS_CHECK_EQUAL(eqn_ids.size(), rCondition.GetGeometry().PointsNumber());
        for (std::size_t i = 0; i < eqn_ids.size(); ++i)
        {
            KRATOS_CHECK_EQUAL(eqn_ids[i], i + 1);
        }
    };

    RansEvmKEpsilonEpsilonWall2D2N_EvaluateTest(evaluation_method);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonWall2D2N_GetDofList, KratosRansFastSuite)
{
    auto evaluation_method = [](Condition& rCondition,
                                ProcessInfo& rProcessInfo, const bool IsSlip) {
        auto dofs = Element::DofsVectorType{};
        rCondition.GetDofList(dofs, rProcessInfo);
        KRATOS_CHECK_EQUAL(dofs.size(), rCondition.GetGeometry().PointsNumber());
        for (std::size_t i = 0; i < dofs.size(); ++i)
        {
            KRATOS_CHECK_EQUAL(dofs[i]->GetVariable(), TURBULENT_ENERGY_DISSIPATION_RATE);
            KRATOS_CHECK_EQUAL(dofs[i]->EquationId(), i + 1);
        }
    };

    RansEvmKEpsilonEpsilonWall2D2N_EvaluateTest(evaluation_method);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonWall2D2N_CalculateLocalSystem, KratosRansFastSuite)
{
    auto evaluation_method = [](Condition& rCondition,
                                ProcessInfo& rProcessInfo, const bool IsSlip) {
        Matrix LHS;
        Vector RHS;
        rCondition.CalculateLocalSystem(LHS, RHS, rProcessInfo);
        KRATOS_CHECK_VECTOR_EQUAL(RHS, ZeroVector(2));
        KRATOS_CHECK_MATRIX_EQUAL(LHS, ZeroMatrix(2, 2));
    };

    RansEvmKEpsilonEpsilonWall2D2N_EvaluateTest(evaluation_method);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonWall2D2N_CalculateRightHandSide, KratosRansFastSuite)
{
    auto evaluation_method = [](Condition& rCondition,
                                ProcessInfo& rProcessInfo, const bool IsSlip) {
        Vector RHS;
        rCondition.CalculateRightHandSide(RHS, rProcessInfo);
        KRATOS_CHECK_VECTOR_EQUAL(RHS, ZeroVector(2));
    };

    RansEvmKEpsilonEpsilonWall2D2N_EvaluateTest(evaluation_method);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonWall2D2N_CalculateLocalVelocityContribution,
                          KratosRansFastSuite)
{
    auto evaluation_method = [](Condition& rCondition,
                                ProcessInfo& rProcessInfo, const bool IsSlip) {
        Matrix LHS;
        Vector RHS;
        rCondition.CalculateLocalVelocityContribution(LHS, RHS, rProcessInfo);

        if (IsSlip)
        {
            Matrix ref_LHS(2, 2);
            ref_LHS(0, 0) = -8.8842866931915732e-01;
            ref_LHS(0, 1) = -3.5618553499123140e-01;
            ref_LHS(1, 0) = -3.5618553499123140e-01;
            ref_LHS(1, 1) = -5.3631347064576806e-01;
            Vector ref_RHS(2);
            ref_RHS[0] = 1.8364623708570253e+01;
            ref_RHS[1] = 1.2045480307045217e+01;
            KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
            KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
        }
        else
        {
            KRATOS_CHECK_VECTOR_EQUAL(RHS, ZeroVector(2));
            KRATOS_CHECK_MATRIX_EQUAL(LHS, ZeroMatrix(2, 2));
        }
    };

    RansEvmKEpsilonEpsilonWall2D2N_EvaluateTest(evaluation_method);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonWall2D2N_CalculateMassMatrix, KratosRansFastSuite)
{
    auto evaluation_method = [](Condition& rCondition,
                                ProcessInfo& rProcessInfo, const bool IsSlip) {
        Matrix M;
        rCondition.CalculateMassMatrix(M, rProcessInfo);
        KRATOS_CHECK_EQUAL(M.size1(), 0);
        KRATOS_CHECK_EQUAL(M.size2(), 0);
    };

    RansEvmKEpsilonEpsilonWall2D2N_EvaluateTest(evaluation_method);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonEpsilonWall2D2N_CalculateDampingMatrix, KratosRansFastSuite)
{
    auto evaluation_method = [](Condition& rCondition,
                                ProcessInfo& rProcessInfo, const bool IsSlip) {
        Matrix LHS;
        rCondition.CalculateDampingMatrix(LHS, rProcessInfo);

        if (IsSlip)
        {
            Matrix ref_LHS(2, 2);
            ref_LHS(0, 0) = -8.8842866931915732e-01;
            ref_LHS(0, 1) = -3.5618553499123140e-01;
            ref_LHS(1, 0) = -3.5618553499123140e-01;
            ref_LHS(1, 1) = -5.3631347064576806e-01;
            KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
        }
        else
        {
            KRATOS_CHECK_MATRIX_EQUAL(LHS, ZeroMatrix(2, 2));
        }
    };

    RansEvmKEpsilonEpsilonWall2D2N_EvaluateTest(evaluation_method);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVmsMonolithicWall2D2N_CalculateLocalSystem, KratosRansFastSuite)
{
    auto evaluation_method = [](Condition& rCondition, const ProcessInfo& rProcessInfo,
                                const bool IsSlip, const bool IsCoSolvingProcessActive) {
        Matrix LHS;
        Vector RHS;
        rCondition.CalculateLocalSystem(LHS, RHS, rProcessInfo);
        KRATOS_CHECK_VECTOR_EQUAL(RHS, ZeroVector(6));
        KRATOS_CHECK_MATRIX_EQUAL(LHS, ZeroMatrix(6, 6));
    };

    RansEvmKEpsilonVmsMonolithicWall2D2N_EvaluateTest(evaluation_method);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVmsMonolithicWall2D2N_CalculateLocalVelocityContribution,
                          KratosRansFastSuite)
{
    auto evaluation_method = [](Condition& rCondition, const ProcessInfo& rProcessInfo,
                                const bool IsSlip, const bool IsCoSolvingProcessActive) {
        Matrix LHS;
        Vector RHS;
        rCondition.CalculateLocalVelocityContribution(LHS, RHS, rProcessInfo);

        if (!IsSlip)
        {
            KRATOS_CHECK_VECTOR_EQUAL(RHS, ZeroVector(6));
            KRATOS_CHECK_MATRIX_EQUAL(LHS, ZeroMatrix(6, 6));
        }
        else if (IsSlip && !IsCoSolvingProcessActive)
        {
            Matrix ref_LHS(6, 6);
            ref_LHS.clear();
            ref_LHS(0, 0) = 9.9785637332056212e-02;
            ref_LHS(1, 1) = 9.9785637332056212e-02;
            Vector ref_RHS(6);
            ref_RHS.clear();
            ref_RHS[0] = 1.8460342906430399e-01;
            ref_RHS[1] = 2.1453912026392086e-01;

            KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
            KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
        }
        else if (IsSlip && IsCoSolvingProcessActive)
        {
            Matrix ref_LHS(6, 6);
            ref_LHS.clear();
            ref_LHS(0, 0) = 5.8700089206809922e+00;
            ref_LHS(0, 3) = 1.9261013323368170e+00;
            ref_LHS(1, 1) = 5.8700089206809922e+00;
            ref_LHS(1, 4) = 1.9261013323368170e+00;
            ref_LHS(3, 0) = 1.9261013323368170e+00;
            ref_LHS(3, 3) = 1.8343964086662754e+00;
            ref_LHS(4, 1) = 1.9261013323368170e+00;
            ref_LHS(4, 4) = 1.8343964086662754e+00;
            Vector ref_RHS(6);
            ref_RHS.clear();
            ref_RHS[0] = 1.2361875542482554e+01;
            ref_RHS[1] = 1.2736085259404343e+01;
            ref_RHS[3] = 4.9941166635828065e+00;
            ref_RHS[4] = 4.2511816490441330e+00;

            KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
            KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
        }
    };

    RansEvmKEpsilonVmsMonolithicWall2D2N_EvaluateTest(evaluation_method);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVmsMonolithicWall2D2N_CalculateMassMatrix, KratosRansFastSuite)
{
    auto evaluation_method = [](Condition& rCondition, const ProcessInfo& rProcessInfo,
                                const bool IsSlip, const bool IsCoSolvingProcessActive) {
        Matrix M;
        rCondition.CalculateMassMatrix(M, rProcessInfo);
        KRATOS_CHECK_EQUAL(M.size1(), 0);
        KRATOS_CHECK_EQUAL(M.size2(), 0);
    };

    RansEvmKEpsilonVmsMonolithicWall2D2N_EvaluateTest(evaluation_method);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKEpsilonVmsMonolithicWall2D2N_CalculateDampingMatrix, KratosRansFastSuite)
{
    auto evaluation_method = [](Condition& rCondition, const ProcessInfo& rProcessInfo,
                                const bool IsSlip, const bool IsCoSolvingProcessActive) {
        Matrix LHS;
        rCondition.CalculateDampingMatrix(LHS, rProcessInfo);

        if (!IsSlip)
        {
            KRATOS_CHECK_MATRIX_EQUAL(LHS, ZeroMatrix(6, 6));
        }
        else if (IsSlip && !IsCoSolvingProcessActive)
        {
            Matrix ref_LHS(6, 6);
            ref_LHS.clear();
            ref_LHS(0, 0) = 9.9785637332056212e-02;
            ref_LHS(1, 1) = 9.9785637332056212e-02;
            KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
        }
        else if (IsSlip && IsCoSolvingProcessActive)
        {
            Matrix ref_LHS(6, 6);
            ref_LHS.clear();
            ref_LHS(0, 0) = 5.8700089206809922e+00;
            ref_LHS(0, 3) = 1.9261013323368170e+00;
            ref_LHS(1, 1) = 5.8700089206809922e+00;
            ref_LHS(1, 4) = 1.9261013323368170e+00;
            ref_LHS(3, 0) = 1.9261013323368170e+00;
            ref_LHS(3, 3) = 1.8343964086662754e+00;
            ref_LHS(4, 1) = 1.9261013323368170e+00;
            ref_LHS(4, 4) = 1.8343964086662754e+00;
            KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
        }
    };

    RansEvmKEpsilonVmsMonolithicWall2D2N_EvaluateTest(evaluation_method);
}

} // namespace Testing
} // namespace Kratos.
