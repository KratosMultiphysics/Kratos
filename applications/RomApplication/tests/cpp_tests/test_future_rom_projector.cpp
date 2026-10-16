//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Raul Bravo
//

#ifdef KRATOS_USE_FUTURE

// System includes

/* External includes */

/* Project includes */
#include "containers/model.h"
#include "geometries/line_2d_2.h"
#include "includes/kratos_parameters.h"
#include "testing/testing.h"
#include "future/containers/define_linear_algebra_serial.h"
#include "future/solving_strategies/schemes/static_scheme.h"
#include "future/solving_strategies/strategies/implicit_strategy_data.h"

/* Application includes */
#include "future/rom_projector.h"

namespace Kratos::Testing {
namespace FutureRomProjectorTestingInternal {

using SchemeType = Future::StaticScheme<Future::SerialLinearAlgebraTraits>;
using StrategyDataType = Future::ImplicitStrategyData<Future::SerialLinearAlgebraTraits>;
using RomProjectorType = Future::RomProjector<Future::SerialLinearAlgebraTraits>;

// Same element as in the test of the ROMBuilderAndSolver
class DummyLaplacianElement final : public Element
{
public:
    DummyLaplacianElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& rCurrentProcessInfo) const override
    {
        rResult.clear();
        rResult.reserve(NumNodes);

        for (const auto& r_node : GetGeometry()) {
            rResult.push_back(r_node.GetDof(TEMPERATURE).EquationId());
        }
    }

    void GetDofList(DofsVectorType& rElementalDofList,
                    const ProcessInfo& rCurrentProcessInfo) const override
    {
        rElementalDofList.clear();
        rElementalDofList.reserve(NumNodes);

        for (const auto& r_node : GetGeometry()) {
            rElementalDofList.push_back(r_node.pGetDof(TEMPERATURE));
        }
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override
    {
        if (rLeftHandSideMatrix.size1() != NumNodes) {
            rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);
        }
        if (rRightHandSideVector.size() != NumNodes) {
            rRightHandSideVector.resize(NumNodes, false);
        }

        // Laplacian
        BoundedMatrix<double, NumNodes, NumNodes> A;
        A(0,0) = 1;
        A(1,0) = -1;
        A(0,1) = -1;
        A(1,1) = 1;

        // Uniform source term
        BoundedVector<double, NumNodes> b;
        b(0) = 0.5;
        b(1) = 0.5;

        // Prev step solution
        BoundedVector<double, NumNodes> x;
        x[0] = GetGeometry()[0].GetSolutionStepValue(TEMPERATURE);
        x[1] = GetGeometry()[1].GetSolutionStepValue(TEMPERATURE);

        noalias(rLeftHandSideMatrix) = A;
        noalias(rRightHandSideVector) = b - prod(A,x);
    }

private:
    static constexpr IndexType NumNodes = 2;
};

ModelPart& FillModel(Model& rModel)
{
    auto& r_model_part = rModel.CreateModelPart("main");
    r_model_part.CreateNewProperties(0);
    r_model_part.SetBufferSize(1);
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)->AddDof(TEMPERATURE);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0)->AddDof(TEMPERATURE);
    r_model_part.CreateNewNode(3, 2.0, 0.0, 0.0)->AddDof(TEMPERATURE);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(TEMPERATURE) = 300;
    }

    // Elements
    auto p_properties = r_model_part.pGetProperties(0);
    for (std::size_t i = 1; i <= 2; ++i) {
        auto p_geometry = Kratos::make_shared<Line2D2<Node>>(r_model_part.pGetNode(i), r_model_part.pGetNode(i+1));
        r_model_part.AddElement(Kratos::make_intrusive<DummyLaplacianElement>(i, p_geometry, p_properties));
    }

    // Dirichlet
    r_model_part.GetNode(1).Fix(TEMPERATURE);

    return r_model_part;
}

}

KRATOS_TEST_CASE_IN_SUITE(FutureRomProjectorGalerkin, RomApplicationFastSuite)
{
    using namespace FutureRomProjectorTestingInternal;

    Model model{};
    ModelPart& r_model_part = FillModel(model);

    Parameters scheme_settings(R"(
    {
        "build_settings" : {
            "name" : "block_builder"
        }
    }
    )");
    auto p_scheme = Kratos::make_shared<SchemeType>(r_model_part, scheme_settings);
    auto p_strategy_data = Kratos::make_shared<StrategyDataType>();
    p_scheme->Initialize(*p_strategy_data);
    p_scheme->InitializeSolutionStep(*p_strategy_data);
    p_scheme->InitializeNonLinIteration(*p_strategy_data);

    // Basis with the modes 1 and x
    const std::size_t n_dofs = 3;
    RomProjectorType::EigenDynamicMatrix phi(n_dofs, 2);
    RomProjectorType::EigenDynamicVector solution(n_dofs);
    for (const auto& r_node : r_model_part.Nodes()) {
        const std::size_t eq_id = r_node.GetDof(TEMPERATURE).EffectiveEquationId();
        phi(eq_id, 0) = 1.0;
        phi(eq_id, 1) = r_node.X();
        solution[eq_id] = r_node.FastGetSolutionStepValue(TEMPERATURE);
    }

    RomProjectorType rom_projector(p_scheme, p_strategy_data);
    rom_projector.BuildEffectiveSystem();
    rom_projector.Project(phi);
    const RomProjectorType::EigenDynamicVector dq = rom_projector.SolveReduced();

    // Same values as in the test of the ROMBuilderAndSolver
    KRATOS_EXPECT_EQ(p_strategy_data->pGetEffectiveDofSet()->size(), n_dofs);
    KRATOS_EXPECT_NEAR(dq(0), 1.0, 1e-8);
    KRATOS_EXPECT_NEAR(dq(1), 0.5, 1e-8);

    // Testing that the solution is set in the free DOFs only
    solution += phi * dq;
    rom_projector.SetSolution(solution);
    KRATOS_EXPECT_NEAR(r_model_part.GetNode(1).FastGetSolutionStepValue(TEMPERATURE), 300.0, 1e-8);
    KRATOS_EXPECT_NEAR(r_model_part.GetNode(2).FastGetSolutionStepValue(TEMPERATURE), 301.5, 1e-8);
    KRATOS_EXPECT_NEAR(r_model_part.GetNode(3).FastGetSolutionStepValue(TEMPERATURE), 302.0, 1e-8);
}

} // namespace Kratos::Testing

#endif // KRATOS_USE_FUTURE
