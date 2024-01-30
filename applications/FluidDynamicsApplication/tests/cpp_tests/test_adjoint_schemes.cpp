//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

// System includes
#include <functional>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/cfd_variables.h"
#include "includes/expect.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "spaces/ublas_space.h"
#include "testing/testing.h"
#include "utilities/normal_calculation_utils.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_response_functions/velocity_pressure_norm_square_response_function.h"
#include "custom_strategies/schemes/residualbased_predictorcorrector_velocity_bossak_scheme_turbulent.h"
#include "custom_strategies/schemes/residualbased_simple_steady_scheme.h"
#include "custom_strategies/schemes/simple_steady_adjoint_scheme.h"
#include "custom_strategies/schemes/velocity_bossak_adjoint_scheme.h"
#include "custom_strategies/schemes/simple_steady_sensitivity_builder_scheme.h"
#include "custom_strategies/schemes/velocity_bossak_sensitivity_builder_scheme.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
namespace Testing
{

void SetModelPartOldValues(
    ModelPart& rModelPart)
{
    // this correspods to t = t - 1 time step in forward time
    const double p_1 = 2.3;
    const double p_2 = 3.5;
    const double p_3 = 6.2;

    const array_1d<double, 3> v_1{2.3, 4.2, 3.6};
    const array_1d<double, 3> v_2{3.1, 2.4, 6.7};
    const array_1d<double, 3> v_3{2.5, 1.8, 6.3};
    const array_1d<double, 3> a_1{4.9, 3.5, 7.9};
    const array_1d<double, 3> a_2{1.7, 3.0, 3.6};
    const array_1d<double, 3> a_3{2.4, 2.0, 4.7};

    rModelPart.GetNode(1).FastGetSolutionStepValue(VELOCITY) = v_1;
    rModelPart.GetNode(1).FastGetSolutionStepValue(PRESSURE) = p_1;
    rModelPart.GetNode(1).FastGetSolutionStepValue(ACCELERATION) = a_1;
    rModelPart.GetNode(2).FastGetSolutionStepValue(VELOCITY) = v_2;
    rModelPart.GetNode(2).FastGetSolutionStepValue(PRESSURE) = p_2;
    rModelPart.GetNode(2).FastGetSolutionStepValue(ACCELERATION) = a_2;
    rModelPart.GetNode(3).FastGetSolutionStepValue(VELOCITY) = v_3;
    rModelPart.GetNode(3).FastGetSolutionStepValue(PRESSURE) = p_3;
    rModelPart.GetNode(3).FastGetSolutionStepValue(ACCELERATION) = a_3;
}

void SetModelPartCurrentValues(
    ModelPart& rModelPart)
{
    // this correspods to t = t time step in forward time
    // this is where we calculate testing derivatives
    const double p_1 = 1.5;
    const double p_2 = 2.7;
    const double p_3 = 1.9;

    const array_1d<double, 3> v_1{1.2, 3.1, 2.4};
    const array_1d<double, 3> v_2{2.2, 3.5, 3.8};
    const array_1d<double, 3> v_3{1.8, 2.3, 2.5};
    const array_1d<double, 3> a_1{3.2, 5.7, 8.4};
    const array_1d<double, 3> a_2{2.4, 3.9, 7.8};
    const array_1d<double, 3> a_3{3.9, 2.6, 9.5};

    rModelPart.GetNode(1).FastGetSolutionStepValue(VELOCITY) = v_1;
    rModelPart.GetNode(1).FastGetSolutionStepValue(PRESSURE) = p_1;
    rModelPart.GetNode(1).FastGetSolutionStepValue(ACCELERATION) = a_1;
    rModelPart.GetNode(2).FastGetSolutionStepValue(VELOCITY) = v_2;
    rModelPart.GetNode(2).FastGetSolutionStepValue(PRESSURE) = p_2;
    rModelPart.GetNode(2).FastGetSolutionStepValue(ACCELERATION) = a_2;
    rModelPart.GetNode(3).FastGetSolutionStepValue(VELOCITY) = v_3;
    rModelPart.GetNode(3).FastGetSolutionStepValue(PRESSURE) = p_3;
    rModelPart.GetNode(3).FastGetSolutionStepValue(ACCELERATION) = a_3;
}

void SetModelPartNextalues(
    ModelPart& rModelPart)
{
    // this correspods to t = t + 1 time step in forward time
    // this is required for adjoints since adjoints are calculated backwards in time,
    // hence these values are loaded to fille the buffer of adjoint model part
    const double p_1 = 7.5;
    const double p_2 = 9.7;
    const double p_3 = 2.9;

    const array_1d<double, 3> v_1{3.4, 7.1, 2.4};
    const array_1d<double, 3> v_2{2.7, 4.7, 4.8};
    const array_1d<double, 3> v_3{2.4, 6.3, 3.5};
    const array_1d<double, 3> a_1{3.9, 7.9, 4.8};
    const array_1d<double, 3> a_2{6.3, 4.5, 5.7};
    const array_1d<double, 3> a_3{4.6, 3.5, 3.2};

    rModelPart.GetNode(1).FastGetSolutionStepValue(VELOCITY) = v_1;
    rModelPart.GetNode(1).FastGetSolutionStepValue(PRESSURE) = p_1;
    rModelPart.GetNode(1).FastGetSolutionStepValue(ACCELERATION) = a_1;
    rModelPart.GetNode(2).FastGetSolutionStepValue(VELOCITY) = v_2;
    rModelPart.GetNode(2).FastGetSolutionStepValue(PRESSURE) = p_2;
    rModelPart.GetNode(2).FastGetSolutionStepValue(ACCELERATION) = a_2;
    rModelPart.GetNode(3).FastGetSolutionStepValue(VELOCITY) = v_3;
    rModelPart.GetNode(3).FastGetSolutionStepValue(PRESSURE) = p_3;
    rModelPart.GetNode(3).FastGetSolutionStepValue(ACCELERATION) = a_3;

    // set adjoint solution values
    const array_1d<double, 3> a_v_1_1{3.2, 1.1, 4.4};
    const array_1d<double, 3> a_v_2_1{4.8, 8.3, 5.1};
    const array_1d<double, 3> a_v_3_1{2.7, 7.5, 6.2};
    const array_1d<double, 3> a_v_1_2{1.3, 2.1, 4.4};
    const array_1d<double, 3> a_v_2_2{4.2, 5.3, 4.1};
    const array_1d<double, 3> a_v_3_2{6.2, 7.7, 3.2};
    const array_1d<double, 3> a_v_1_3{3.9, 2.9, 1.5};
    const array_1d<double, 3> a_v_2_3{4.8, 4.5, 3.3};
    const array_1d<double, 3> a_v_3_3{2.4, 3.1, 2.2};

    rModelPart.GetNode(1).FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1) = a_v_1_1;
    rModelPart.GetNode(1).FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_2) = a_v_2_1;
    rModelPart.GetNode(1).FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3) = a_v_3_1;
    rModelPart.GetNode(2).FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1) = a_v_1_2;
    rModelPart.GetNode(2).FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_2) = a_v_2_2;
    rModelPart.GetNode(2).FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3) = a_v_3_2;
    rModelPart.GetNode(3).FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1) = a_v_1_3;
    rModelPart.GetNode(3).FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_2) = a_v_2_3;
    rModelPart.GetNode(3).FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3) = a_v_3_3;
}

void AddPrimalDofs(ModelPart& rModelPart)
{
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(VELOCITY_X);
        r_node.AddDof(VELOCITY_Y);
        r_node.AddDof(PRESSURE);
    }
}

void SetPrimalModelPartValues(
    ModelPart& rModelPart)
{
    AddPrimalDofs(rModelPart);
    double previous_time = 1.2;
    double current_time = 1.6;

    rModelPart.CloneTimeStep(previous_time);
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    rModelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.3);

    for (auto& r_node : rModelPart.Nodes()) {
        r_node.FastGetSolutionStepValue(DENSITY) = 1.2;
        r_node.FastGetSolutionStepValue(VISCOSITY) = 2.5;
        r_node.SetValue(Y_WALL, 2.0);
    }

    SetModelPartOldValues(rModelPart);
    rModelPart.CloneTimeStep(current_time);
    SetModelPartCurrentValues(rModelPart);
}

void AddAdjointDofs(ModelPart& rModelPart)
{
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(ADJOINT_FLUID_VECTOR_1_X);
        r_node.AddDof(ADJOINT_FLUID_VECTOR_1_Y);
        r_node.AddDof(ADJOINT_FLUID_SCALAR_1);
    }
}

void UpdateAdjointModelPartRelaxedAcceleration(
    ModelPart& rAdjointModelPart,
    const ModelPart& rPrimalModelPart,
    const double AlphaBossak)
{
    // compute relaxed acceleration
    // TODO: remove this when VMS elements no longer depend on Acceleration for adjoints.
    //       Currently, primal_output_process of h5 is used to calculate relaxed acceleration
    //       and to overwrite Acceleration with the relaxed acceleration. This is planned
    //       to be removed, and replaced with RELAXED_ACCELERATION variable so ACCELERATION
    //       and RELAXED_ACCELERATION will be available.
    const auto& r_primal_nodes = rPrimalModelPart.Nodes();
    auto& r_adjoint_nodes = rAdjointModelPart.Nodes();
    for (unsigned int i = 0; i < r_primal_nodes.size(); ++i) {
        const auto& r_primal_node = *(r_primal_nodes.begin() + i);
        auto& r_adjoint_node = *(r_adjoint_nodes.begin() + i);

        const auto& tn_minus_1 = r_primal_node.FastGetSolutionStepValue(ACCELERATION, 1);
        const auto& tn_0 = r_primal_node.FastGetSolutionStepValue(ACCELERATION);
        const auto& tn_plus_1 = r_adjoint_node.FastGetSolutionStepValue(ACCELERATION, 1);

        r_adjoint_node.FastGetSolutionStepValue(ACCELERATION) =
            (1 - AlphaBossak) * tn_0 + AlphaBossak * tn_minus_1;
        r_adjoint_node.FastGetSolutionStepValue(ACCELERATION, 1) =
            (1 - AlphaBossak) * tn_plus_1 + AlphaBossak * tn_0;
    }
}

void SetAdjointModelPartValues(
    ModelPart& rAdjointPrimalModelPart,
    const ModelPart& rPrimalModelPart,
    const double AlphaBossak = -0.3)
{
    AddAdjointDofs(rAdjointPrimalModelPart);

    double current_time = 1.6;
    double next_time = 2.0;

    rAdjointPrimalModelPart.CloneTimeStep(next_time);
    rAdjointPrimalModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    rAdjointPrimalModelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.3);

    for (auto& r_node : rAdjointPrimalModelPart.Nodes()) {
        r_node.FastGetSolutionStepValue(DENSITY) = 1.2;
        r_node.FastGetSolutionStepValue(VISCOSITY) = 2.5;
        r_node.SetValue(Y_WALL, 2.0);
    }

    SetModelPartNextalues(rAdjointPrimalModelPart);
    rAdjointPrimalModelPart.CloneTimeStep(current_time);
    SetModelPartCurrentValues(rAdjointPrimalModelPart);
    UpdateAdjointModelPartRelaxedAcceleration(rAdjointPrimalModelPart, rPrimalModelPart, AlphaBossak);
}

void GenerateTestModelPart(
    ModelPart& rModelPart,
    const std::string& rElementName,
    const std::string& rConditionName)
{
    // Set buffer size
    rModelPart.SetBufferSize(2);

    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(DENSITY);
    rModelPart.AddNodalSolutionStepVariable(REACTION);
    rModelPart.AddNodalSolutionStepVariable(PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(BODY_FORCE);
    rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(NORMAL);
    rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_1);
    rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_2);
    rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_3);
    rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_SCALAR_1);
    rModelPart.AddNodalSolutionStepVariable(AUX_ADJOINT_FLUID_VECTOR_1);

    auto p_elem_prop = rModelPart.CreateNewProperties(1);

    // Element creation
    auto p_node_1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_node_3 = rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);

    std::vector<ModelPart::IndexType> cond_nodes{1, 2};
    std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3};

    p_node_1->Set(SLIP, true);
    p_node_2->Set(SLIP, true);

    p_node_1->SetValue(Y_WALL, 2.0);
    p_node_2->SetValue(Y_WALL, 3.0);

    rModelPart.CreateNewElement(rElementName, 1, elem_nodes, p_elem_prop);
    auto p_cond = rModelPart.CreateNewCondition(rConditionName, 1, cond_nodes, p_elem_prop);
    p_cond->Set(SLIP, true);

    auto& sub_model_part = rModelPart.CreateSubModelPart("Structure");
    sub_model_part.AddNodes({1, 2, 3});
    sub_model_part.AddElements({1});
}

template<class TSparseSpaceType, class TDenseSpaceType>
void EvaluateSystemContributions(
    Scheme<TSparseSpaceType, TDenseSpaceType>& rScheme,
    Matrix& rLHS,
    Vector& rRHS,
    ModelPart& rModelPart,
    const bool UpdateSolution = false)
{
    typename Scheme<TSparseSpaceType, TDenseSpaceType>::TSystemMatrixType A;
    typename Scheme<TSparseSpaceType, TDenseSpaceType>::TSystemVectorType b;
    typename Scheme<TSparseSpaceType, TDenseSpaceType>::TSystemVectorType x;
    typename Scheme<TSparseSpaceType, TDenseSpaceType>::DofsArrayType dx;

    std::vector<std::size_t> ids;

    NormalCalculationUtils().CalculateOnSimplex(rModelPart, 2);

    rScheme.Initialize(rModelPart);
    rScheme.InitializeElements(rModelPart);
    rScheme.InitializeConditions(rModelPart);
    rScheme.InitializeSolutionStep(rModelPart, A, x, b);
    rScheme.InitializeNonLinIteration(rModelPart, A, x, b);

    rScheme.CalculateSystemContributions(rModelPart.GetElement(1), rLHS, rRHS,
                                         ids, rModelPart.GetProcessInfo());

    Matrix aux_lhs;
    Vector aux_rhs;
    rScheme.CalculateSystemContributions(rModelPart.GetCondition(1), aux_lhs, aux_rhs,
                                         ids, rModelPart.GetProcessInfo());

    // do the assembly
    for (unsigned int i = 0; i < aux_lhs.size1(); ++i) {
        for (unsigned int j = 0; j < aux_lhs.size2(); ++j) {
            rLHS(i, j) += aux_lhs(i, j);
        }
        rRHS[i] += aux_rhs[i];
    }

    if (UpdateSolution) {
        rScheme.Update(rModelPart, dx, A, x, b);
    }
    rScheme.FinalizeNonLinIteration(rModelPart, A, x, b);
    rScheme.FinalizeSolutionStep(rModelPart, A, x, b);
}

Vector EvaluateSteadyResidual(
    ModelPart& rModelPart)
{
    using SparseSpaceType =
        UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

    ResidualBasedSimpleSteadyScheme<SparseSpaceType, LocalSpaceType> primal_scheme(
        1.0, 1.0, 2);

    Matrix lhs;
    Vector rhs;
    EvaluateSystemContributions(primal_scheme, lhs, rhs, rModelPart);

    return rhs;
}

Vector EvaluateBossakResidual(
    ModelPart& rModelPart)
{
    using SparseSpaceType =
        UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

    ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<SparseSpaceType, LocalSpaceType> primal_scheme(
        -0.3, 0.0, 2);

    Matrix lhs;
    Vector rhs;
    EvaluateSystemContributions(primal_scheme, lhs, rhs, rModelPart);

    return rhs;
}

template<class TSensitivityBuilderSchemeType, class TResponseType>
void EvaluateResidualShapeSensitivities(
    Matrix& rOutput,
    ModelPart& rModelPart,
    TSensitivityBuilderSchemeType& rSensitivityBuilderScheme,
    TResponseType& rResponse)
{
    Vector adjoint_values;
    GlobalPointersVector<ModelPart::NodeType> gp_vector;
    NormalCalculationUtils().CalculateOnSimplex(rModelPart, 2);
    NormalCalculationUtils().CalculateNormalShapeDerivativesOnSimplex(
        rModelPart.Conditions(), 2);
    rSensitivityBuilderScheme.Initialize(rModelPart, rModelPart, rResponse);
    rSensitivityBuilderScheme.InitializeSolutionStep(rModelPart, rModelPart, rResponse);

    rSensitivityBuilderScheme.CalculateResidualSensitivityMatrix(
        rModelPart.GetElement(1), adjoint_values, rOutput, gp_vector,
        SHAPE_SENSITIVITY, rModelPart.GetProcessInfo());

    Matrix aux_matrix;
    rSensitivityBuilderScheme.CalculateResidualSensitivityMatrix(
        rModelPart.GetCondition(1), adjoint_values, aux_matrix, gp_vector,
        SHAPE_SENSITIVITY, rModelPart.GetProcessInfo());

    // do the assembly
    for (unsigned int i = 0; i < aux_matrix.size1(); ++i) {
        for (unsigned int j = 0; j < aux_matrix.size2(); ++j) {
            rOutput(i, j) += aux_matrix(i, j);
        }
    }
}

void UpdateSensitivities(
    Matrix& rOutput,
    const Vector& rResidual,
    const Vector& rResidual_0,
    const double Delta,
    const unsigned int RowIndex)
{
    row(rOutput, RowIndex) = (rResidual - rResidual_0) / Delta;
}

void UpdateSensitivities(
    Vector& rOutput,
    const double rValue,
    const double rValue_0,
    const double Delta,
    const unsigned int RowIndex)
{
    rOutput[RowIndex] = (rValue - rValue_0) / Delta;
}

void ResizeOutput(
    Matrix& rOutput,
    const unsigned int Size)
{
    rOutput.resize(Size, 9);
}

void ResizeOutput(
    Vector& rOutput,
    const unsigned int Size)
{
    rOutput.resize(Size);
}

template<class TOutputType, class TResponseFunctor>
void CalculateFiniteDifferenceSensitivities(
    TOutputType& rOutput,
    ModelPart& rModelPart,
    const std::vector<const Variable<double>*>& rVariables,
    const TResponseFunctor& rResidualFunction,
    const double Delta,
    const unsigned int DerivativeStep = 0)
{
    const auto& rhs_0 = rResidualFunction(rModelPart);

    const int number_of_nodes = rModelPart.NumberOfNodes();
    const int number_of_derivatives = rVariables.size();
    ResizeOutput(rOutput, number_of_nodes * number_of_derivatives);
    unsigned int local_index = 0;
    for (auto& r_node : rModelPart.Nodes()) {
        for (const auto p_variable : rVariables) {
            if (*p_variable == SHAPE_SENSITIVITY_X) {
                r_node.X() += Delta;
            } else if (*p_variable == SHAPE_SENSITIVITY_Y) {
                r_node.Y() += Delta;
            } else if (*p_variable == SHAPE_SENSITIVITY_Z) {
                r_node.Z() += Delta;
            } else {
                r_node.FastGetSolutionStepValue(*p_variable, DerivativeStep) += Delta;
            }

            const auto& rhs = rResidualFunction(rModelPart);
            UpdateSensitivities(rOutput, rhs, rhs_0, Delta, local_index++);

            if (*p_variable == SHAPE_SENSITIVITY_X) {
                r_node.X() -= Delta;
            } else if (*p_variable == SHAPE_SENSITIVITY_Y) {
                r_node.Y() -= Delta;
            } else if (*p_variable == SHAPE_SENSITIVITY_Z) {
                r_node.Z() -= Delta;
            } else {
                r_node.FastGetSolutionStepValue(*p_variable, DerivativeStep) -= Delta;
            }
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(SimpleSteadySensitivityBuilderScheme, FluidDynamicsApplicationFastSuite)
{
    Model model;

    // Create a primal test element inside a modelpart
    auto& primal_model_part = model.CreateModelPart("Primal");
    GenerateTestModelPart(primal_model_part, "VMS2D3N", "MonolithicWallCondition2D2N");
    SetPrimalModelPartValues(primal_model_part);

    // Create a adjoint test element inside a modelpart
    auto& adjoint_model_part = model.CreateModelPart("Adjoint");
    GenerateTestModelPart(adjoint_model_part, "VMSAdjointElement2D", "AdjointMonolithicWallCondition2D2N");
    SetAdjointModelPartValues(adjoint_model_part, primal_model_part);

    VariableUtils().SetHistoricalVariableToZero(ACCELERATION, primal_model_part.Nodes());
    VariableUtils().SetHistoricalVariableToZero(ACCELERATION, adjoint_model_part.Nodes());

    // compute adjoint sensitivities
    VelocityPressureNormSquareResponseFunction response(
        Parameters(R"(
        {
            "main_model_part_name": "Adjoint",
            "norm_model_part_name": "Adjoint"
        })"),
        model);

    SimpleSteadySensitivityBuilderScheme sensitivity_builder_scheme(2, 3);
    Matrix analytical_sensitivities;
    EvaluateResidualShapeSensitivities(analytical_sensitivities, adjoint_model_part,
                                       sensitivity_builder_scheme, response);

    // compute finite difference sensitivities
    const double delta = 1e-8;
    Matrix fd_sensitivities;
    CalculateFiniteDifferenceSensitivities(fd_sensitivities, primal_model_part,
                                           {&SHAPE_SENSITIVITY_X, &SHAPE_SENSITIVITY_Y},
                                           EvaluateSteadyResidual, delta);

    KRATOS_EXPECT_MATRIX_NEAR(fd_sensitivities, analytical_sensitivities, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(SimpleSteadyAdjointScheme, FluidDynamicsApplicationFastSuite)
{
    using SparseSpaceType =
        UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

    Model model;

    // Create a primal test element inside a modelpart
    auto& primal_model_part = model.CreateModelPart("Primal");
    GenerateTestModelPart(primal_model_part, "VMS2D3N", "MonolithicWallCondition2D2N");
    SetPrimalModelPartValues(primal_model_part);

    // Create a adjoint test element inside a modelpart
    auto& adjoint_model_part = model.CreateModelPart("Adjoint");
    GenerateTestModelPart(adjoint_model_part, "VMSAdjointElement2D", "AdjointMonolithicWallCondition2D2N");
    SetAdjointModelPartValues(adjoint_model_part, primal_model_part);

    VariableUtils().SetHistoricalVariableToZero(ACCELERATION, primal_model_part.Nodes());
    VariableUtils().SetHistoricalVariableToZero(ACCELERATION, adjoint_model_part.Nodes());

    // compute adjoint sensitivities
    VelocityPressureNormSquareResponseFunction response(
        Parameters(R"(
        {
            "main_model_part_name": "Adjoint",
            "norm_model_part_name": "Adjoint"
        })"),
        model);
    SimpleSteadyAdjointScheme<SparseSpaceType, LocalSpaceType> adjoint_scheme(
        Kratos::make_shared<VelocityPressureNormSquareResponseFunction>(response), 2, 3);

    Matrix lhs;
    Vector rhs;
    EvaluateSystemContributions(adjoint_scheme, lhs, rhs, adjoint_model_part);

    // compute finite difference sensitivities
    const double delta = 1e-8;
    Matrix fd_sensitivities;
    CalculateFiniteDifferenceSensitivities(fd_sensitivities, primal_model_part,
                                           {&VELOCITY_X, &VELOCITY_Y, &PRESSURE},
                                           EvaluateSteadyResidual, delta);

    KRATOS_EXPECT_MATRIX_NEAR(fd_sensitivities, lhs, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(VelocityBossakSensitivityBuilderScheme, FluidDynamicsApplicationFastSuite)
{
    Model model;

    // Create a primal test element inside a modelpart
    auto& primal_model_part = model.CreateModelPart("Primal");
    GenerateTestModelPart(primal_model_part, "VMS2D3N", "MonolithicWallCondition2D2N");
    SetPrimalModelPartValues(primal_model_part);

    // Create a adjoint test element inside a modelpart
    auto& adjoint_model_part = model.CreateModelPart("Adjoint");
    GenerateTestModelPart(adjoint_model_part, "VMSAdjointElement2D", "AdjointMonolithicWallCondition2D2N");
    SetAdjointModelPartValues(adjoint_model_part, primal_model_part);

    // compute adjoint sensitivities
    VelocityPressureNormSquareResponseFunction response(
        Parameters(R"(
        {
            "main_model_part_name": "Adjoint",
            "norm_model_part_name": "Adjoint"
        })"),
        model);
    Matrix analytical_sensitivities;
    const double alpha_bossak = -0.3;
    VelocityBossakSensitivityBuilderScheme sensitivity_builder_scheme(alpha_bossak, 2, 3);
    EvaluateResidualShapeSensitivities(analytical_sensitivities, adjoint_model_part,
                                       sensitivity_builder_scheme, response);

    // compute finite difference sensitivities
    const double delta = 1e-8;
    Matrix fd_sensitivities;
    CalculateFiniteDifferenceSensitivities(fd_sensitivities, primal_model_part,
                                           {&SHAPE_SENSITIVITY_X, &SHAPE_SENSITIVITY_Y},
                                           EvaluateBossakResidual, delta);

    KRATOS_EXPECT_MATRIX_NEAR(fd_sensitivities, analytical_sensitivities, 1e-6);
}


KRATOS_TEST_CASE_IN_SUITE(VelocityBossakAdjointSchemeLHS, FluidDynamicsApplicationFastSuite)
{
    using SparseSpaceType =
        UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

    Model model;

    // Create a primal test element inside a modelpart
    auto& primal_model_part = model.CreateModelPart("Primal");
    GenerateTestModelPart(primal_model_part, "VMS2D3N", "MonolithicWallCondition2D2N");
    SetPrimalModelPartValues(primal_model_part);
    const double delta_time = primal_model_part.GetProcessInfo()[DELTA_TIME];

    // Create a adjoint test element inside a modelpart
    auto& adjoint_model_part = model.CreateModelPart("Adjoint");
    GenerateTestModelPart(adjoint_model_part, "VMSAdjointElement2D", "AdjointMonolithicWallCondition2D2N");
    SetAdjointModelPartValues(adjoint_model_part, primal_model_part);

    // compute adjoint sensitivities
    VelocityPressureNormSquareResponseFunction response(
        Parameters(R"(
        {
            "main_model_part_name": "Adjoint",
            "norm_model_part_name": "Adjoint"
        })"),
        model);
    VelocityBossakAdjointScheme<SparseSpaceType, LocalSpaceType> adjoint_scheme(
        Parameters(R"({
            "name"         : "adjoint_bossak",
            "scheme_type"  : "bossak",
            "alpha_bossak" : -0.3
        })"),
        Kratos::make_shared<VelocityPressureNormSquareResponseFunction>(response), 2, 3);

    Matrix lhs;
    Vector rhs;
    EvaluateSystemContributions(adjoint_scheme, lhs, rhs, adjoint_model_part);

    // compute finite difference sensitivities
    const double delta = 1e-7;
    Matrix fd_sensitivities_1;
    CalculateFiniteDifferenceSensitivities(fd_sensitivities_1, primal_model_part,
                                           {&VELOCITY_X, &VELOCITY_Y, &PRESSURE},
                                           EvaluateBossakResidual, delta);

    Matrix fd_sensitivities_2;
    // using DISPLACEMENT_X as a dummy variable. because, there is no variable for
    // PRESSURE time derivative.
    CalculateFiniteDifferenceSensitivities(fd_sensitivities_2, primal_model_part,
                                           {&ACCELERATION_X, &ACCELERATION_Y, &DISPLACEMENT_X},
                                           EvaluateBossakResidual, delta);

    const double alpha_bossak = -0.3;
    const double beta = 0.25 * (1.0 - alpha_bossak) * (1.0 - alpha_bossak);
    const double gamma = 0.5 - alpha_bossak;
    const double c6 = gamma / (beta * delta_time);
    const double c7 = 1.0 / (delta_time * delta_time * beta);

    fd_sensitivities_1 *= c6;
    fd_sensitivities_2 *= c7;

    noalias(fd_sensitivities_1) += fd_sensitivities_2;

    KRATOS_EXPECT_MATRIX_NEAR(fd_sensitivities_1, lhs, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(VelocityBossakAdjointSchemeRHS, FluidDynamicsApplicationFastSuite)
{
    using SparseSpaceType =
        UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

    Model model;

    Parameters primal_response_parameters(R"({
            "main_model_part_name" : "Primal",
            "norm_model_part_name" : "Primal.Structure"
        })");

    // Create a primal test element inside a modelpart
    auto& primal_model_part = model.CreateModelPart("Primal");
    GenerateTestModelPart(primal_model_part, "VMS2D3N", "MonolithicWallCondition2D2N");
    SetPrimalModelPartValues(primal_model_part);
    const double delta_time = primal_model_part.GetProcessInfo()[DELTA_TIME];
    VelocityPressureNormSquareResponseFunction primal_response(primal_response_parameters, model);
    primal_response.Initialize();

    // Create a adjoint test element inside a modelpart
    auto& adjoint_model_part = model.CreateModelPart("Adjoint");
    GenerateTestModelPart(adjoint_model_part, "VMSAdjointElement2D", "AdjointMonolithicWallCondition2D2N");
    SetAdjointModelPartValues(adjoint_model_part, primal_model_part);

    Parameters adjoint_response_parameters(R"({
            "main_model_part_name" : "Adjoint",
            "norm_model_part_name" : "Adjoint.Structure"
        })");

    // compute adjoint sensitivities
    VelocityPressureNormSquareResponseFunction adjoint_response(
        adjoint_response_parameters, model);
    adjoint_response.Initialize();

    VelocityBossakAdjointScheme<SparseSpaceType, LocalSpaceType> adjoint_scheme(
        Parameters(R"({
            "name"         : "adjoint_bossak",
            "scheme_type"  : "bossak",
            "alpha_bossak" : -0.3
        })"),
        Kratos::make_shared<VelocityPressureNormSquareResponseFunction>(adjoint_response), 2, 3);

    Matrix lhs;
    Vector rhs;
    EvaluateSystemContributions(adjoint_scheme, lhs, rhs, adjoint_model_part, true);

    const std::vector<const Variable<double>*> first_derivative_variables = {
        &VELOCITY_X, &VELOCITY_Y, &PRESSURE};
    const std::vector<const Variable<double>*> second_derivative_variables = {
        &ACCELERATION_X, &ACCELERATION_Y, &DISPLACEMENT_X};

    // compute time scheme constants
    const double alpha_bossak = -0.3;
    const double beta = 0.25 * (1.0 - alpha_bossak) * (1.0 - alpha_bossak);
    const double gamma = 0.5 - alpha_bossak;
    const double c0 = 1.0 - gamma / beta;
    const double c1 = 1.0 / (beta * delta_time);
    const double c2 = (1.0 - 0.5 * gamma / beta) * delta_time;
    const double c3 = (1.0 - 0.5 / beta);
    const double c4 = (beta - gamma * (gamma + 0.5)) / (beta * beta * delta_time);
    const double c5 = -1.0 * (gamma + 0.5) / (beta * beta * delta_time * delta_time);
    const double c6 = gamma / (beta * delta_time);
    const double c7 = 1.0 / (delta_time * delta_time * beta);

    const double delta = 1e-7;

    // compute objective finite difference sensitivities
    const auto& objective_Value_function = [&](ModelPart& rModelPart) -> double {
        return primal_response.CalculateValue(rModelPart);
    };
    Vector objective_first_derivatives;
    CalculateFiniteDifferenceSensitivities(
        objective_first_derivatives, primal_model_part,
        first_derivative_variables, objective_Value_function, delta);

    Vector objective_second_derivatives;
    CalculateFiniteDifferenceSensitivities(
        objective_second_derivatives, primal_model_part,
        second_derivative_variables, objective_Value_function, delta);

    // compute finite difference sensitivities
    Matrix residual_first_derivatives;
    CalculateFiniteDifferenceSensitivities(
        residual_first_derivatives, primal_model_part,
        first_derivative_variables, EvaluateBossakResidual, delta);

    Matrix residual_second_derivatives;
    // using DISPLACEMENT_X as a dummy variable. because, there is no variable
    // for PRESSURE time derivative.
    CalculateFiniteDifferenceSensitivities(
        residual_second_derivatives, primal_model_part,
        second_derivative_variables, EvaluateBossakResidual, delta);

    // compute old adjoint values
    Vector lambda_1(9, 0.0), lambda_2(9, 0.0), lambda_3(9, 0.0), lambda_1_old(9, 0.0), lambda_2_old(9, 0.0), lambda_3_old(9, 0.0);
    for (unsigned int a = 0; a < 3; ++a) {
        const auto& vec_1 = adjoint_model_part.GetNode(a + 1).FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1);
        const auto& vec_2 = adjoint_model_part.GetNode(a + 1).FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_2);
        const auto& vec_3 = adjoint_model_part.GetNode(a + 1).FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3);

        const auto& vec_old_1 = adjoint_model_part.GetNode(a + 1).FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1, 1);
        const auto& vec_old_2 = adjoint_model_part.GetNode(a + 1).FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_2, 1);
        const auto& vec_old_3 = adjoint_model_part.GetNode(a + 1).FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3, 1);
        for (unsigned int i = 0; i < 2; ++i) {
            lambda_1[a * 3 + i] = vec_1[i];
            lambda_2[a * 3 + i] = vec_2[i];
            lambda_3[a * 3 + i] = vec_3[i];

            lambda_1_old[a * 3 + i] = vec_old_1[i];
            lambda_2_old[a * 3 + i] = vec_old_2[i];
            lambda_3_old[a * 3 + i] = vec_old_3[i];
        }
        lambda_1[a * 3 + 2] = adjoint_model_part.GetNode(a + 1).FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1);
        lambda_1_old[a * 3 + 2] = adjoint_model_part.GetNode(a + 1).FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1, 1);
    }

    // compute and check lambda_2 using finite difference approach
    Vector fd_lambda_2 = -objective_first_derivatives -
                         prod(residual_first_derivatives, lambda_1) +
                         lambda_2_old * c0 - lambda_3_old * c1;

    // making the followings zero because, these values are not stored in ADJOINT_FLUID_VECTOR_2 variable, and those contributions are
    // directly taken to RHS from  objective_first_derivatives and prod(residual_first_derivatives, lambda_1)
    fd_lambda_2[2] = fd_lambda_2[5] = fd_lambda_2[8] = 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(fd_lambda_2, lambda_2, 1e-6);

    // compute and check lambda_3 using finite difference approach
    Vector fd_lambda_3 = -objective_second_derivatives -
                         prod(residual_second_derivatives, lambda_1) +
                         lambda_2_old * c2 + lambda_3_old * c3;
    // making the followings zero because, these values are not stored in ADJOINT_FLUID_VECTOR_2 variable, and those contributions are
    // directly taken to RHS from  objective_first_derivatives and prod(residual_first_derivatives, lambda_1)
    fd_lambda_3[2] = fd_lambda_3[5] = fd_lambda_3[8] = 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(fd_lambda_3, lambda_3, 1e-6);

    // compute and check rhs using finite difference sensitivities
    Vector fd_rhs = -objective_first_derivatives * c6 - objective_second_derivatives * c7 +
                    lambda_2_old * c4 + lambda_3_old * c5 - prod(lhs, lambda_1);
    KRATOS_EXPECT_VECTOR_NEAR(fd_rhs, rhs, 1e-6);
}

} // namespace Testing
} // namespace Kratos.
