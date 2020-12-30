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

// External includes

// Project includes
#include "containers/model.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "testing/testing.h"
#include "utilities/normal_calculation_utils.h"

// Application includes
#include "custom_response_functions/drag_response_function.h"
#include "custom_strategies/schemes/residualbased_simple_steady_scheme.h"
#include "custom_strategies/schemes/simple_steady_adjoint_scheme.h"
#include "custom_strategies/schemes/simple_steady_sensitivity_builder_scheme.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
namespace Testing
{
void SetModelPartValues(ModelPart& rModelPart)
{
    double delta_time = 0.1;
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    rModelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.3);
    rModelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);

    const double p_1 = 1.5;
    const double p_2 = 2.7;
    const double p_3 = 1.9;
    array_1d<double, 3> v_1 = ZeroVector(3);
    array_1d<double, 3> v_2 = ZeroVector(3);
    array_1d<double, 3> v_3 = ZeroVector(3);
    v_1[0] = 1.3;
    v_2[0] = 2.0;
    v_3[0] = 3.0;
    v_3[1] = 0.5;
    rModelPart.GetNode(1).FastGetSolutionStepValue(VELOCITY) = v_1;
    rModelPart.GetNode(1).FastGetSolutionStepValue(PRESSURE) = p_1;
    rModelPart.GetNode(2).FastGetSolutionStepValue(VELOCITY) = v_2;
    rModelPart.GetNode(2).FastGetSolutionStepValue(PRESSURE) = p_2;
    rModelPart.GetNode(3).FastGetSolutionStepValue(VELOCITY) = v_3;
    rModelPart.GetNode(3).FastGetSolutionStepValue(PRESSURE) = p_3;

    for (auto& r_node : rModelPart.Nodes()) {
        r_node.FastGetSolutionStepValue(DENSITY) = 1.2;
        r_node.FastGetSolutionStepValue(VISCOSITY) = 2.5;
    }
}

void GenerateTestModelPart(ModelPart& rModelPart, const std::string& rElementName)
{
    // Set buffer size
    rModelPart.SetBufferSize(1);

    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(DENSITY);
    rModelPart.AddNodalSolutionStepVariable(REACTION);
    rModelPart.AddNodalSolutionStepVariable(PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(BODY_FORCE);
    rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(NORMAL);
    rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_1);
    rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_SCALAR_1);

    auto p_elem_prop = rModelPart.CreateNewProperties(1);

    // Element creation
    auto p_node_1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_node_3 = rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);

    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(VELOCITY_X);
        r_node.AddDof(VELOCITY_Y);
        r_node.AddDof(PRESSURE);
        r_node.AddDof(ADJOINT_FLUID_VECTOR_1_X);
        r_node.AddDof(ADJOINT_FLUID_VECTOR_1_Y);
        r_node.AddDof(ADJOINT_FLUID_SCALAR_1);
    }

    std::vector<ModelPart::IndexType> cond_nodes{1, 2};
    std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3};

    p_node_1->Set(SLIP, true);
    p_node_2->Set(SLIP, true);

    rModelPart.CreateNewElement(rElementName, 1, elem_nodes, p_elem_prop);
    auto p_cond = rModelPart.CreateNewCondition("LineCondition2D2N", 1, cond_nodes, p_elem_prop);
    p_cond->Set(SLIP, true);

    SetModelPartValues(rModelPart);
}

void EvaluateResidual(
    Vector& rResidual,
    ModelPart& rModelPart)
{
    using SparseSpaceType =
        UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

    typename SparseSpaceType::MatrixType A;
    typename SparseSpaceType::VectorType b;
    typename SparseSpaceType::VectorType x;

    ResidualBasedSimpleSteadyScheme<SparseSpaceType, LocalSpaceType> primal_scheme(
        1.0, 1.0, 2);

    Matrix lhs;
    std::vector<std::size_t> ids;

    NormalCalculationUtils().CalculateOnSimplex(rModelPart, 2);

    primal_scheme.Initialize(rModelPart);
    primal_scheme.InitializeElements(rModelPart);
    primal_scheme.InitializeConditions(rModelPart);
    primal_scheme.InitializeSolutionStep(rModelPart, A, x, b);
    primal_scheme.InitializeNonLinIteration(rModelPart, A, x, b);
    primal_scheme.CalculateSystemContributions(
        rModelPart.GetElement(1), lhs, rResidual, ids,
        rModelPart.GetProcessInfo());
}

KRATOS_TEST_CASE_IN_SUITE(SimpleSteadySensitivityBuilderScheme, FluidDynamicsApplicationFastSuite)
{
    Model model;

    // Create a primal test element inside a modelpart
    auto& primal_model_part = model.CreateModelPart("Primal");
    GenerateTestModelPart(primal_model_part, "VMS2D3N");

    // Create a adjoint test element inside a modelpart
    auto& adjoint_model_part = model.CreateModelPart("Adjoint");
    GenerateTestModelPart(adjoint_model_part, "VMSAdjointElement2D");
    adjoint_model_part.GetProcessInfo()[DELTA_TIME] *= -1.0;

    Vector rhs_0;
    EvaluateResidual(rhs_0, primal_model_part);

    // compute adjoint sensitivities
    DragResponseFunction<2> response(Parameters(""
                                                "{}"
                                                ""),
                                     adjoint_model_part);
    Vector adjoint_values;
    Matrix analytical_sensitivities;
    GlobalPointersVector<ModelPart::NodeType> gp_vector;
    SimpleSteadySensitivityBuilderScheme<2> sensitivity_builder_scheme;
    NormalCalculationUtils().CalculateOnSimplex(adjoint_model_part, 2);
    sensitivity_builder_scheme.Initialize(adjoint_model_part, adjoint_model_part, response);
    sensitivity_builder_scheme.InitializeSolutionStep(
        adjoint_model_part, adjoint_model_part, response);
    sensitivity_builder_scheme.CalculateResidualSensitivityMatrix(
        adjoint_model_part.GetElement(1), adjoint_values, analytical_sensitivities,
        gp_vector, SHAPE_SENSITIVITY, adjoint_model_part.GetProcessInfo());

    // compute finite difference sensitivities
    Vector rhs;
    const double delta = 1e-8;
    Matrix fd_sensitivities = ZeroMatrix(6, 9);
    unsigned int local_index = 0;
    for (auto& r_node : primal_model_part.Nodes()) {
        r_node.X() += delta;
        EvaluateResidual(rhs, primal_model_part);
        row(fd_sensitivities, local_index++) = (rhs - rhs_0) / delta;
        r_node.X() -= delta;

        r_node.Y() += delta;
        EvaluateResidual(rhs, primal_model_part);
        row(fd_sensitivities, local_index++) = (rhs - rhs_0) / delta;
        r_node.Y() -= delta;
    }

    KRATOS_CHECK_MATRIX_NEAR(fd_sensitivities, analytical_sensitivities, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(SimpleSteadyAdjointScheme, FluidDynamicsApplicationFastSuite)
{
    using SparseSpaceType =
        UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

    Model model;

    // Create a primal test element inside a modelpart
    auto& primal_model_part = model.CreateModelPart("Primal");
    GenerateTestModelPart(primal_model_part, "VMS2D3N");

    // Create a adjoint test element inside a modelpart
    auto& adjoint_model_part = model.CreateModelPart("Adjoint");
    GenerateTestModelPart(adjoint_model_part, "VMSAdjointElement2D");
    adjoint_model_part.GetProcessInfo()[DELTA_TIME] *= -1.0;

    typename SparseSpaceType::MatrixType A;
    typename SparseSpaceType::VectorType b;
    typename SparseSpaceType::VectorType x;

    Vector rhs_0;
    EvaluateResidual(rhs_0, primal_model_part);

    // compute adjoint sensitivities
    DragResponseFunction<2> response(Parameters(""
                                                "{}"
                                                ""),
                                     adjoint_model_part);
    std::vector<std::size_t> ids;
    Matrix analytical_sensitivities;
    SimpleSteadyAdjointScheme<2, SparseSpaceType, LocalSpaceType> adjoint_scheme(
        Kratos::make_shared<DragResponseFunction<2>>(response));

    NormalCalculationUtils().CalculateOnSimplex(adjoint_model_part, 2);
    adjoint_scheme.InitializeElements(adjoint_model_part);
    adjoint_scheme.InitializeConditions(adjoint_model_part);
    adjoint_scheme.InitializeSolutionStep(adjoint_model_part, A, x, b);
    adjoint_scheme.InitializeNonLinIteration(adjoint_model_part, A, x, b);
    adjoint_scheme.CalculateLHSContribution(adjoint_model_part.GetElement(1),
                                            analytical_sensitivities, ids,
                                            adjoint_model_part.GetProcessInfo());

    // compute finite difference sensitivities
    Vector rhs;
    const double delta = 1e-8;
    Matrix fd_sensitivities = ZeroMatrix(9, 9);
    unsigned int local_index = 0;
    for (auto& r_node : primal_model_part.Nodes()) {
        auto& r_velocity = r_node.FastGetSolutionStepValue(VELOCITY);
        for (unsigned int k = 0; k < 2; ++k) {
            r_velocity[k] += delta;

            EvaluateResidual(rhs, primal_model_part);
            row(fd_sensitivities, local_index++) = (rhs - rhs_0) / delta;

            r_velocity[k] -= delta;
        }

        auto& r_pressure = r_node.FastGetSolutionStepValue(PRESSURE);
        r_pressure += delta;

        EvaluateResidual(rhs, primal_model_part);
        row(fd_sensitivities, local_index++) = (rhs - rhs_0) / delta;

        r_pressure -= delta;
    }

    KRATOS_CHECK_MATRIX_NEAR(fd_sensitivities, analytical_sensitivities, 1e-6);
}

} // namespace Testing
} // namespace Kratos.
