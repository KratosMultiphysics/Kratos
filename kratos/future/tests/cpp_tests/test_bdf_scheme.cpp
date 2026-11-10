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

// System includes
#include <limits>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "testing/testing.h"

#ifdef KRATOS_USE_FUTURE
#include "future/containers/define_linear_algebra_serial.h"
#include "future/solving_strategies/schemes/bdf_scheme.h"
#endif

namespace Kratos::Testing
{

namespace
{

/**
 * @brief A simple 1D dynamic test element
 * @details It provides a known stiffness (K), mass (M) and damping (D) matrix and reads the
 * displacement, velocity and acceleration from the nodal variables. This allows to verify the
 * dynamic (BDF) build of the scheme in a fully deterministic way:
 *  - LHS = K + c0^2*M + c0*D
 *  - RHS = (f - K*u) - M*a - D*v
 */
class TestDynamicElement1D : public Element
{
public:

    using IndexType = Element::IndexType;

    using DofsVectorType = Element::DofsVectorType;

    using EquationIdVectorType = Element::EquationIdVectorType;

    TestDynamicElement1D(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {}

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
        CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    }

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rLeftHandSideMatrix.size1() != n_nodes || rLeftHandSideMatrix.size2() != n_nodes) {
            rLeftHandSideMatrix.resize(n_nodes, n_nodes, false);
        }
        const double k = 1.0;
        rLeftHandSideMatrix(0, 0) = 2.0 * k;
        rLeftHandSideMatrix(0, 1) = -k;
        rLeftHandSideMatrix(1, 0) = -k;
        rLeftHandSideMatrix(1, 1) = 2.0 * k;
    }

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rRightHandSideVector.size() != n_nodes) {
            rRightHandSideVector.resize(n_nodes);
        }
        const double k = 1.0;
        const double f_ext = 1.0;
        const double u0 = r_geom[0].GetSolutionStepValue(DISPLACEMENT_X);
        const double u1 = r_geom[1].GetSolutionStepValue(DISPLACEMENT_X);
        rRightHandSideVector(0) = f_ext - (2.0 * k * u0 - k * u1);
        rRightHandSideVector(1) = f_ext - (-k * u0 + 2.0 * k * u1);
    }

    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rMassMatrix.size1() != n_nodes || rMassMatrix.size2() != n_nodes) {
            rMassMatrix.resize(n_nodes, n_nodes, false);
        }
        const double m = 1.0;
        rMassMatrix(0, 0) = m;
        rMassMatrix(0, 1) = 0.0;
        rMassMatrix(1, 0) = 0.0;
        rMassMatrix(1, 1) = m;
    }

    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rDampingMatrix.size1() != n_nodes || rDampingMatrix.size2() != n_nodes) {
            rDampingMatrix.resize(n_nodes, n_nodes, false);
        }
        const double c = 0.5;
        rDampingMatrix(0, 0) = c;
        rDampingMatrix(0, 1) = 0.0;
        rDampingMatrix(1, 0) = 0.0;
        rDampingMatrix(1, 1) = c;
    }

    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rValues.size() != n_nodes) {
            rValues.resize(n_nodes);
        }
        for (SizeType i = 0; i < n_nodes; ++i) {
            rValues[i] = r_geom[i].GetSolutionStepValue(VELOCITY_X, Step);
        }
    }

    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rValues.size() != n_nodes) {
            rValues.resize(n_nodes);
        }
        for (SizeType i = 0; i < n_nodes; ++i) {
            rValues[i] = r_geom[i].GetSolutionStepValue(ACCELERATION_X, Step);
        }
    }

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rResult.size() != n_nodes) {
            rResult.resize(n_nodes);
        }
        for (SizeType i = 0; i < n_nodes; ++i) {
            rResult[i] = r_geom[i].GetDof(DISPLACEMENT_X).EquationId();
        }
    }

    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rElementalDofList.size() != n_nodes) {
            rElementalDofList.resize(n_nodes);
        }
        for (SizeType i = 0; i < n_nodes; ++i) {
            rElementalDofList[i] = r_geom[i].pGetDof(DISPLACEMENT_X);
        }
    }
};

/**
 * @brief Sets up a 1D dynamic test model part
 * @details It creates two nodes connected by a single TestDynamicElement1D, adds the DISPLACEMENT
 * (as DOF) and the VELOCITY/ACCELERATION (as nodal state variables), and sets a time step of 0.5.
 */
void SetUpDynamicTestModelPart1D(ModelPart& rModelPart)
{
    const int domain_size = 2;
    const std::size_t buffer_size = 3;
    rModelPart.SetBufferSize(buffer_size);
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, domain_size);
    rModelPart.GetProcessInfo().SetValue(DELTA_TIME, 0.5);

    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(ACCELERATION);

    auto p_node_a = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_b = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_geom = Kratos::make_shared<Line2D2<Node>>(p_node_a, p_node_b);
    rModelPart.AddElement(Kratos::make_intrusive<TestDynamicElement1D>(1, p_geom));

    rModelPart.GetNodalSolutionStepVariablesList().AddDof(&DISPLACEMENT_X);
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
    }
}

} // End of anonymous namespace.

KRATOS_TEST_CASE_IN_SUITE(BDFSchemeCoefficientsOrder1, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    SetUpDynamicTestModelPart1D(r_test_model_part);

    // Create the scheme (order 1)
    Parameters scheme_settings = Parameters(R"({
        "integration_order" : 1
    })");
    using SchemeType = Future::BDFScheme<Future::SerialLinearAlgebraTraits>;
    auto p_scheme = Kratos::make_unique<SchemeType>(r_test_model_part, scheme_settings);

    // Compute the BDF coefficients for dt = 0.5
    Future::ImplicitStrategyData<Future::SerialLinearAlgebraTraits> strategy_data_container;
    p_scheme->Initialize(strategy_data_container);
    p_scheme->InitializeSolutionStep(strategy_data_container);

    // Order 1 with dt = 0.5 gives c0 = 1/dt = 2 and c1 = -1/dt = -2
    KRATOS_CHECK_EQUAL(p_scheme->GetOrder(), static_cast<std::size_t>(1));
    const auto& r_coefficients = p_scheme->GetBDFCoefficients();
    KRATOS_CHECK_EQUAL(r_coefficients.size(), 2);
    KRATOS_CHECK_NEAR(r_coefficients[0], 2.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_coefficients[1], -2.0, 1.0e-12);
#else
    true;
#endif
}

KRATOS_TEST_CASE_IN_SUITE(BDFSchemeBuildWithDynamicsOrder1, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    SetUpDynamicTestModelPart1D(r_test_model_part);

    // Set the previous step (step 1) values from which the current step will be extrapolated
    Node& r_node_1 = r_test_model_part.GetNode(1);
    Node& r_node_2 = r_test_model_part.GetNode(2);
    r_node_1.FastGetSolutionStepValue(DISPLACEMENT, 1)[0] = 0.1;
    r_node_1.FastGetSolutionStepValue(VELOCITY, 1)[0] = 0.3;
    r_node_1.FastGetSolutionStepValue(ACCELERATION, 1)[0] = 0.5;
    r_node_2.FastGetSolutionStepValue(DISPLACEMENT, 1)[0] = 0.2;
    r_node_2.FastGetSolutionStepValue(VELOCITY, 1)[0] = 0.4;
    r_node_2.FastGetSolutionStepValue(ACCELERATION, 1)[0] = 0.6;

    // Create the scheme (order 1)
    Parameters scheme_settings = Parameters(R"({
        "integration_order" : 1,
        "echo_level" : 0
    })");
    using SchemeType = Future::BDFScheme<Future::SerialLinearAlgebraTraits>;
    auto p_scheme = Kratos::make_unique<SchemeType>(r_test_model_part, scheme_settings);

    Future::ImplicitStrategyData<Future::SerialLinearAlgebraTraits> strategy_data_container;
    p_scheme->Initialize(strategy_data_container);
    p_scheme->InitializeSolutionStep(strategy_data_container);
    p_scheme->Predict(strategy_data_container);

    // Build the system
    const auto p_linear_system = strategy_data_container.pGetLinearSystem();
    auto& r_lhs = *(p_linear_system->pGetMatrix(Future::LinearSystemTags::SparseMatrixTag::LHS));
    auto& r_rhs = *(p_linear_system->pGetVector(Future::LinearSystemTags::DenseVectorTag::RHS));
    p_scheme->Build(r_lhs, r_rhs);

    // Expected values for order 1 (c0 = 2):
    // u0 = [0.3125, 0.475], v0 = [0.425, 0.55], a0 = [0.25, 0.3]
    // LHS = K + c0^2*M + c0*D = [[7, -1], [-1, 7]]
    // RHS = (f - K*u0) - M*a0 - D*v0 = [0.3875, -0.2125]
    const double tol = 1.0e-12;
    KRATOS_CHECK_EQUAL(r_lhs.size1(), 2);
    KRATOS_CHECK_NEAR(r_lhs(0, 0), 7.0, tol);
    KRATOS_CHECK_NEAR(r_lhs(0, 1), -1.0, tol);
    KRATOS_CHECK_NEAR(r_lhs(1, 0), -1.0, tol);
    KRATOS_CHECK_NEAR(r_lhs(1, 1), 7.0, tol);
    std::vector<double> expected_rhs = {0.3875, -0.2125};
    KRATOS_CHECK_VECTOR_NEAR(r_rhs, expected_rhs, tol);
#else
    true;
#endif
}

KRATOS_TEST_CASE_IN_SUITE(BDFSchemeBuildWithDynamicsOrder2, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    SetUpDynamicTestModelPart1D(r_test_model_part);

    // Set the current (step 0) values: the converged state of the previous time step (t = 0)
    Node& r_node_1 = r_test_model_part.GetNode(1);
    Node& r_node_2 = r_test_model_part.GetNode(2);
    r_node_1.FastGetSolutionStepValue(DISPLACEMENT, 0)[0] = 0.1;
    r_node_1.FastGetSolutionStepValue(VELOCITY, 0)[0] = 0.3;
    r_node_1.FastGetSolutionStepValue(ACCELERATION, 0)[0] = 0.5;
    r_node_2.FastGetSolutionStepValue(DISPLACEMENT, 0)[0] = 0.2;
    r_node_2.FastGetSolutionStepValue(VELOCITY, 0)[0] = 0.4;
    r_node_2.FastGetSolutionStepValue(ACCELERATION, 0)[0] = 0.6;

    // Mark the current step as a time step and advance to t = 0.5.
    // The previous time step (required by the 2nd order BDF) then holds TIME = 0 and DELTA_TIME = 0.5
    r_test_model_part.GetProcessInfo().SetValue(TIME, 0.0);
    r_test_model_part.GetProcessInfo().SetAsTimeStepInfo();
    r_test_model_part.GetProcessInfo().SetValue(DELTA_TIME, 0.5);
    r_test_model_part.CloneTimeStep(0.5);

    // Set the previous step (step 1) values from which the current step will be extrapolated
    r_node_1.FastGetSolutionStepValue(DISPLACEMENT, 1)[0] = 0.1;
    r_node_1.FastGetSolutionStepValue(VELOCITY, 1)[0] = 0.3;
    r_node_1.FastGetSolutionStepValue(ACCELERATION, 1)[0] = 0.5;
    r_node_2.FastGetSolutionStepValue(DISPLACEMENT, 1)[0] = 0.2;
    r_node_2.FastGetSolutionStepValue(VELOCITY, 1)[0] = 0.4;
    r_node_2.FastGetSolutionStepValue(ACCELERATION, 1)[0] = 0.6;

    // Create the scheme (order 2)
    Parameters scheme_settings = Parameters(R"({
        "integration_order" : 2,
        "echo_level" : 0
    })");
    using SchemeType = Future::BDFScheme<Future::SerialLinearAlgebraTraits>;
    auto p_scheme = Kratos::make_unique<SchemeType>(r_test_model_part, scheme_settings);

    Future::ImplicitStrategyData<Future::SerialLinearAlgebraTraits> strategy_data_container;
    p_scheme->Initialize(strategy_data_container);
    p_scheme->InitializeSolutionStep(strategy_data_container);
    p_scheme->Predict(strategy_data_container);

    // Build the system
    const auto p_linear_system = strategy_data_container.pGetLinearSystem();
    auto& r_lhs = *(p_linear_system->pGetMatrix(Future::LinearSystemTags::SparseMatrixTag::LHS));
    auto& r_rhs = *(p_linear_system->pGetVector(Future::LinearSystemTags::DenseVectorTag::RHS));
    p_scheme->Build(r_lhs, r_rhs);

    // Expected values for order 2 with constant dt = 0.5 (c0 = 3, c1 = -4, c2 = 1):
    // u0 = [0.3125, 0.475], v0 = [0.5375, 0.625], a0 = [0.4125, 0.275]
    // LHS = K + c0^2*M + c0*D = [[12.5, -1], [-1, 12.5]]
    // RHS = (f - K*u0) - M*a0 - D*v0 = [0.16875, -0.225]
    const double tol = 1.0e-12;
    KRATOS_CHECK_EQUAL(r_lhs.size1(), 2);
    KRATOS_CHECK_NEAR(r_lhs(0, 0), 12.5, tol);
    KRATOS_CHECK_NEAR(r_lhs(0, 1), -1.0, tol);
    KRATOS_CHECK_NEAR(r_lhs(1, 0), -1.0, tol);
    KRATOS_CHECK_NEAR(r_lhs(1, 1), 12.5, tol);
    std::vector<double> expected_rhs = {0.16875, -0.225};
    KRATOS_CHECK_VECTOR_NEAR(r_rhs, expected_rhs, tol);
#else
    true;
#endif
}

}  // namespace Kratos::Testing.
