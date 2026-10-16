//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Kratos Team
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/line_2d_2.h"
#include "future/includes/adjoint_derivatives_interface.h"
#include "includes/element.h"
#include "future/includes/linear_adjoint_derivatives.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "testing/testing.h"

namespace Kratos::Testing {

using namespace Kratos::Future;

namespace {

// Linear 2 node element with 1 DOF (DISPLACEMENT_X) per node, and known K and M. No damping.
class LinearTestElement
    : public Element
    , public LinearAdjointDerivatives<LinearTestElement>
{
public:
    LinearTestElement(const IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry) {}

    void GetDofList(DofsVectorType& rDofs, const ProcessInfo&) const override
    {
        rDofs = {GetGeometry()[0].pGetDof(DISPLACEMENT_X), GetGeometry()[1].pGetDof(DISPLACEMENT_X)};
    }

    void CalculateLeftHandSide(MatrixType& rLhs, const ProcessInfo&) override
    {
        rLhs.resize(2, 2, false);
        rLhs(0, 0) = 2.0; rLhs(0, 1) = -1.0;
        rLhs(1, 0) = -1.0; rLhs(1, 1) = 3.0;
    }

    void CalculateMassMatrix(MatrixType& rMass, const ProcessInfo&) override
    {
        rMass.resize(2, 2, false);
        rMass(0, 0) = 4.0; rMass(0, 1) = 1.0;
        rMass(1, 0) = 1.0; rMass(1, 1) = 5.0;
    }

    void CalculateDampingMatrix(MatrixType& rDamping, const ProcessInfo&) override
    {
        rDamping = ZeroMatrix(2, 2);
    }
};

// Coupled entity: owns u (DISPLACEMENT_X) and depends on p (PRESSURE). The implementor defines the order: p first, then u.
class CoupledTestElement
    : public Element
    , public AdjointDerivativesInterface
{
public:
    CoupledTestElement(const IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry) {}

    void GetAdjointDofList(
        const std::vector<DerivativeVariable>& rVariables,
        DofsVectorType& rAdjointDofs,
        const ProcessInfo&) const override
    {
        rAdjointDofs.clear();
        for (const auto& r_variable : rVariables) {
            if (r_variable.GetVariable().Key() == PRESSURE.Key()) {
                for (auto& r_node : GetGeometry()) {
                    rAdjointDofs.push_back(r_node.pGetDof(ADJOINT_PRESSURE));
                }
            }
        }
        for (auto& r_node : GetGeometry()) {
            rAdjointDofs.push_back(r_node.pGetDof(ADJOINT_DISPLACEMENT_X));
        }
    }

    void GetAdjointEquationIdVector(
        const std::vector<DerivativeVariable>& rVariables,
        EquationIdVectorType& rEquationIds,
        const ProcessInfo& rProcessInfo) const override
    {
        DofsVectorType dofs;
        GetAdjointDofList(rVariables, dofs, rProcessInfo);
        rEquationIds.clear();
        for (const auto& p_dof : dofs) {
            rEquationIds.push_back(p_dof->EquationId());
        }
    }
};

// Entity overriding nothing
class EmptyTestElement
    : public Element
    , public AdjointDerivativesInterface
{
public:
    EmptyTestElement(const IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry) {}
};

ModelPart& CreateModelPart(Model& rModel)
{
    auto& r_model_part = rModel.CreateModelPart("test");
    r_model_part.SetBufferSize(2);
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(ADJOINT_DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(ADJOINT_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(ADJOINT_VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ADJOINT_ACCELERATION);
    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    std::size_t id = 0;
    for (auto p_node : {p_node_1, p_node_2}) {
        p_node->AddDof(DISPLACEMENT_X).SetEquationId(id++);
        p_node->AddDof(PRESSURE).SetEquationId(id++);
        p_node->AddDof(ADJOINT_DISPLACEMENT_X).SetEquationId(10 + id);
        p_node->AddDof(ADJOINT_PRESSURE).SetEquationId(20 + id);
    }
    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(LinearAdjointDerivatives, KratosCoreFutureSuite)
{
    Model model;
    auto& r_model_part = CreateModelPart(model);
    auto p_geometry = Kratos::make_shared<Line2D2<Node>>(r_model_part.pGetNode(1), r_model_part.pGetNode(2));
    LinearTestElement element(1, p_geometry);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    // Dofs and ids follow the primal order
    std::vector<DerivativeVariable> variables{DerivativeVariable(DISPLACEMENT_X)};
    AdjointDerivativesInterface::DofsVectorType dofs;
    AdjointDerivativesInterface::EquationIdVectorType ids;
    element.GetAdjointDofList(variables, dofs, r_process_info);
    element.GetAdjointEquationIdVector(variables, ids, r_process_info);
    KRATOS_EXPECT_EQ(dofs.size(), 2u);
    KRATOS_EXPECT_EQ(dofs[0]->Id(), 1u);
    KRATOS_EXPECT_EQ(dofs[1]->Id(), 2u);
    KRATOS_EXPECT_TRUE(dofs[0]->GetVariable() == ADJOINT_DISPLACEMENT_X);
    KRATOS_EXPECT_EQ(ids.size(), 2u);
    KRATOS_EXPECT_EQ(ids[0], dofs[0]->EquationId());
    KRATOS_EXPECT_EQ(ids[1], dofs[1]->EquationId());

    // Values are read from nodal data
    r_model_part.pGetNode(1)->FastGetSolutionStepValue(ADJOINT_DISPLACEMENT_X) = 1.0;
    r_model_part.pGetNode(2)->FastGetSolutionStepValue(ADJOINT_DISPLACEMENT_X) = 2.0;
    r_model_part.pGetNode(1)->FastGetSolutionStepValue(ADJOINT_VELOCITY_X) = 3.0;
    r_model_part.pGetNode(2)->FastGetSolutionStepValue(ADJOINT_VELOCITY_X) = 4.0;
    r_model_part.pGetNode(1)->FastGetSolutionStepValue(ADJOINT_ACCELERATION_X) = 5.0;
    r_model_part.pGetNode(2)->FastGetSolutionStepValue(ADJOINT_ACCELERATION_X) = 6.0;
    Vector values;
    element.GetAdjointValues(variables, values, r_process_info);
    KRATOS_EXPECT_VECTOR_EQ(values, (std::vector<double>{1.0, 2.0}));
    element.GetAdjointFirstDerivativeValues(variables, values, r_process_info);
    KRATOS_EXPECT_VECTOR_EQ(values, (std::vector<double>{3.0, 4.0}));
    element.GetAdjointSecondDerivativeValues(variables, values, r_process_info);
    KRATOS_EXPECT_VECTOR_EQ(values, (std::vector<double>{5.0, 6.0}));

    // Values of the previous step
    r_model_part.pGetNode(1)->FastGetSolutionStepValue(ADJOINT_DISPLACEMENT_X, 1) = 7.0;
    r_model_part.pGetNode(2)->FastGetSolutionStepValue(ADJOINT_DISPLACEMENT_X, 1) = 8.0;
    element.GetAdjointValues(variables, values, r_process_info, 1);
    KRATOS_EXPECT_VECTOR_EQ(values, (std::vector<double>{7.0, 8.0}));

    // Derivatives
    std::vector<DerivativeVariable> request{
        DerivativeVariable(DISPLACEMENT_X),    // u_x
        DerivativeVariable(DISPLACEMENT_Y),    // u_y (no DOF)
        DerivativeVariable(ACCELERATION_X),    // a_x
        DerivativeVariable(DISPLACEMENT_X),    // u_x again
        DerivativeVariable(PRESSURE)           // no DOF
    };
    Matrix stiffness, mass, damping, force;
    element.CalculateResidualTermDerivatives(
        {{ResidualTerm::Stiffness, stiffness}, {ResidualTerm::Inertia, mass}, {ResidualTerm::Damping, damping}, {ResidualTerm::Force, force}},
        request, r_process_info);

    Matrix K(2, 2), M(2, 2);
    K(0, 0) = 2.0; K(0, 1) = -1.0; K(1, 0) = -1.0; K(1, 1) = 3.0;
    M(0, 0) = 4.0; M(0, 1) = 1.0; M(1, 0) = 1.0; M(1, 1) = 5.0;

    // Columns are concatenated over the 5 requested variables, 2 nodes each
    auto block = [&](const Matrix& rMatrix, const std::size_t Variable) {
        return Matrix(project(rMatrix, range(0, 2), range(2 * Variable, 2 * Variable + 2)));
    };
    const Matrix zero = ZeroMatrix(2, 2);
    for (const Matrix* p_matrix : {&stiffness, &mass, &damping}) {
        KRATOS_EXPECT_EQ(p_matrix->size1(), 2u);
        KRATOS_EXPECT_EQ(p_matrix->size2(), 10u);
    }
    // d(K u)/d u_x = K, for the array variable component and for the component variable
    KRATOS_EXPECT_MATRIX_NEAR(block(stiffness, 0), K, 1e-14);
    KRATOS_EXPECT_MATRIX_NEAR(block(stiffness, 3), K, 1e-14);
    // d(M a)/d a_x = M
    KRATOS_EXPECT_MATRIX_NEAR(block(mass, 2), M, 1e-14);
    // Everything else is zero
    KRATOS_EXPECT_MATRIX_NEAR(block(stiffness, 1), zero, 1e-14); // component without DOF
    KRATOS_EXPECT_MATRIX_NEAR(block(stiffness, 2), zero, 1e-14); // K does not depend on a
    KRATOS_EXPECT_MATRIX_NEAR(block(stiffness, 4), zero, 1e-14); // no DOF
    KRATOS_EXPECT_MATRIX_NEAR(block(mass, 0), zero, 1e-14);      // M does not depend on u
    KRATOS_EXPECT_MATRIX_NEAR(block(damping, 0), zero, 1e-14);   // damping is zero
    KRATOS_EXPECT_EQ(force.size1(), 0u);

    // A linear entity gives the same derivatives for the residual of any step
    Matrix stiffness_previous;
    element.CalculateResidualTermDerivatives({{ResidualTerm::Stiffness, stiffness_previous}}, request, r_process_info, 1);
    KRATOS_EXPECT_MATRIX_NEAR(stiffness_previous, stiffness, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(AdjointDerivativesCoupledOrder, KratosCoreFutureSuite)
{
    Model model;
    auto& r_model_part = CreateModelPart(model);
    auto p_geometry = Kratos::make_shared<Line2D2<Node>>(r_model_part.pGetNode(1), r_model_part.pGetNode(2));
    CoupledTestElement element(1, p_geometry);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    // The order is the one chosen by the implementor: p first, then u. Ids follow the same order.
    const std::vector<DerivativeVariable> variables{DerivativeVariable(DISPLACEMENT_X), DerivativeVariable(PRESSURE)};
    AdjointDerivativesInterface::DofsVectorType dofs;
    AdjointDerivativesInterface::EquationIdVectorType ids;
    element.GetAdjointDofList(variables, dofs, r_process_info);
    element.GetAdjointEquationIdVector(variables, ids, r_process_info);
    KRATOS_EXPECT_EQ(dofs.size(), 4u);
    KRATOS_EXPECT_TRUE(dofs[0]->GetVariable() == ADJOINT_PRESSURE);
    KRATOS_EXPECT_TRUE(dofs[1]->GetVariable() == ADJOINT_PRESSURE);
    KRATOS_EXPECT_TRUE(dofs[2]->GetVariable() == ADJOINT_DISPLACEMENT_X);
    KRATOS_EXPECT_TRUE(dofs[3]->GetVariable() == ADJOINT_DISPLACEMENT_X);
    for (std::size_t i = 0; i < dofs.size(); ++i) {
        KRATOS_EXPECT_EQ(ids[i], dofs[i]->EquationId());
    }

    // Without the pressure variable only the own DOFs are listed
    element.GetAdjointDofList({DerivativeVariable(DISPLACEMENT_X)}, dofs, r_process_info);
    KRATOS_EXPECT_EQ(dofs.size(), 2u);
}

KRATOS_TEST_CASE_IN_SUITE(AdjointDerivativesInterfaceNotImplemented, KratosCoreFutureSuite)
{
    Model model;
    auto& r_model_part = CreateModelPart(model);
    auto p_geometry = Kratos::make_shared<Line2D2<Node>>(r_model_part.pGetNode(1), r_model_part.pGetNode(2));
    EmptyTestElement element(1, p_geometry);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    AdjointDerivativesInterface::DofsVectorType dofs;
    AdjointDerivativesInterface::EquationIdVectorType ids;
    Matrix output;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element.GetAdjointDofList({}, dofs, r_process_info), "GetAdjointDofList is not implemented.");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element.GetAdjointEquationIdVector({}, ids, r_process_info), "GetAdjointEquationIdVector is not implemented.");
    Vector values;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element.GetAdjointValues({}, values, r_process_info), "GetAdjointValues is not implemented.");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element.GetAdjointFirstDerivativeValues({}, values, r_process_info), "GetAdjointFirstDerivativeValues is not implemented.");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element.GetAdjointSecondDerivativeValues({}, values, r_process_info), "GetAdjointSecondDerivativeValues is not implemented.");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element.CalculateResidualTermDerivatives({{ResidualTerm::Inertia, output}}, {}, r_process_info), "CalculateResidualTermDerivatives is not implemented.");
}

} // namespace Kratos::Testing
