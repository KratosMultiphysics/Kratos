//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Sebastian Ares de Parga
//

// System includes

/* External includes */

/* Project includes */
#include "containers/model.h"
#include "geometries/triangle_2d_3.h"
#include "includes/kratos_parameters.h"
#include "testing/testing.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "linear_solvers/linear_solver.h"
#include "custom_strategies/rom_builder_and_solver.h"
#include "custom_strategies/petrov_galerkin_rom_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "utilities/math_utils.h"

namespace Kratos::Testing {
namespace PetrovGalerkinROMBuilderAndSolverTestingInternal {

using SparseSpaceType = UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType >;
using BuilderAndSolverType = BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType >;
using PetrovGalerkinROMBuilderAndSolverType = PetrovGalerkinROMBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>;

using SchemeType = Scheme< SparseSpaceType, LocalSpaceType >;
using ResidualBasedIncrementalUpdateStaticSchemeType =  ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType>;

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


ModelPart& FillModel(Model& model)
{
    auto & mpart = model.CreateModelPart("main");
    mpart.CreateNewProperties(0);
    mpart.SetBufferSize(1);
    mpart.AddNodalSolutionStepVariable(TEMPERATURE);

    mpart.CreateNewNode(1, 0.0, 0.0, 0.0)->AddDof(TEMPERATURE);
    mpart.CreateNewNode(2, 1.0, 0.0, 0.0)->AddDof(TEMPERATURE);
    mpart.CreateNewNode(3, 2.0, 0.0, 0.0)->AddDof(TEMPERATURE);

    mpart.GetNode(1).FastGetSolutionStepValue(TEMPERATURE) = 300;
    mpart.GetNode(2).FastGetSolutionStepValue(TEMPERATURE) = 300;
    mpart.GetNode(3).FastGetSolutionStepValue(TEMPERATURE) = 300;

    // Elements
    using LineType = Line2D2<Node>;
    auto pProperties = mpart.pGetProperties(0);
    for(std::size_t i=1; i<=2; ++i)
    {
        auto p_geometry = std::make_shared<LineType>(mpart.pGetNode(i), mpart.pGetNode(i+1));
        auto pelem = Kratos::make_intrusive<DummyLaplacianElement>(i, p_geometry, pProperties);
        mpart.AddElement(pelem);
    }

    // Basis
    const auto phi_1 = [](Node const&){ return 1.0; };
    const auto phi_2 = [](Node const& r_node){ return r_node.X(); };
    Matrix phi_3 = ZeroMatrix(1,3);
    // Setting the 3rd mode to be the residual
    phi_3(0,0) = 0.5;
    phi_3(0,1) = 1.0;
    phi_3(0,2) = 0.5;

    Matrix basis = ZeroMatrix(1, 2); // 1 dof per node, 2 basis
    Matrix left_basis = ZeroMatrix(1, 3); // 1 dof per node, 3 basis
    for(auto& r_node: mpart.Nodes())
    {
        basis(0,0) = phi_1(r_node);
        basis(0,1) = phi_2(r_node);
        r_node.SetValue(ROM_BASIS, basis);

        left_basis(0,0) = phi_1(r_node);
        left_basis(0,1) = phi_2(r_node);
        left_basis(0,2) = phi_3(r_node);
        r_node.SetValue(ROM_LEFT_BASIS, left_basis);
    }

    // Dirichlet
    mpart.GetNode(1).Fix(TEMPERATURE);

    return mpart;
};

SparseSpaceType::VectorType BuildAndSolve(
        ModelPart& rModelPart,
        SchemeType::Pointer pScheme,
        BuilderAndSolverType& rBuilderAndSolver)
{
    rBuilderAndSolver.SetDofSetIsInitializedFlag(false);

    SparseSpaceType::MatrixPointerType pA;
    SparseSpaceType::VectorPointerType pDx;
    SparseSpaceType::VectorPointerType pb;

    rBuilderAndSolver.SetUpDofSet(pScheme, rModelPart);
    rBuilderAndSolver.SetUpSystem(rModelPart);
    rBuilderAndSolver.ResizeAndInitializeVectors(pScheme, pA, pDx, pb, rModelPart);

    SparseSpaceType::MatrixType& rA  = *pA;
    SparseSpaceType::VectorType& rDx = *pDx;
    SparseSpaceType::VectorType& rb  = *pb;

    rBuilderAndSolver.InitializeSolutionStep(rModelPart, rA, rDx, rb);
    pScheme->InitializeSolutionStep(rModelPart, rA, rDx, rb);
    pScheme->InitializeNonLinIteration(rModelPart, rA, rDx, rb);

    rBuilderAndSolver.BuildAndSolve(pScheme, rModelPart, rA, rDx, rb);

    return rDx;
}

}

KRATOS_TEST_CASE_IN_SUITE(PetrovGalerkinROMBuilderAndSolver, RomApplicationFastSuite)
{
    using namespace PetrovGalerkinROMBuilderAndSolverTestingInternal;

    Model model{};
    ModelPart& model_part = FillModel(model);

    Parameters parameters(R"(
    {
        "name" : "rom_builder_and_solver",
        "nodal_unknowns" : ["TEMPERATURE"],
        "number_of_rom_dofs" : 2,
        "petrov_galerkin_number_of_rom_dofs" : 3
    }
    )");

    auto plinear_solver = std::make_shared<LinearSolverType>();
    SchemeType::Pointer p_scheme = std::make_shared<ResidualBasedIncrementalUpdateStaticSchemeType>();
    auto romBnS = PetrovGalerkinROMBuilderAndSolverType(plinear_solver, parameters);

    const auto dx = BuildAndSolve(model_part, p_scheme, romBnS);
    const auto& dq = model_part.GetValue(ROM_SOLUTION_INCREMENT);

    KRATOS_EXPECT_NEAR(model_part.ElementsBegin()->GetValue(HROM_WEIGHT), 1, 1e-8);
    KRATOS_EXPECT_EQ(romBnS.GetEquationSystemSize(), 3);

    KRATOS_EXPECT_NEAR(dq(0), 1.0 , 1e-8);
    KRATOS_EXPECT_NEAR(dq(1), 0.5 , 1e-8);

    // Testing free dofs
    KRATOS_EXPECT_EQ(dx.size(), 3);
    KRATOS_EXPECT_NEAR(dx(1), 1.5 , 1e-8);
    KRATOS_EXPECT_NEAR(dx(2), 2.0 , 1e-8);
}

}