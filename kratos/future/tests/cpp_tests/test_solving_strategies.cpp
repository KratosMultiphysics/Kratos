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
//

// System includes
#include <limits>

// External includes

// Project includes
#include "constraints/linear_master_slave_constraint.h"
#include "containers/model.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "testing/testing.h"

#ifdef KRATOS_USE_FUTURE
#include "future/linear_solvers/amgcl_solver.h"
#include "future/linear_solvers/skyline_lu_factorization_solver.h"
#include "future/solving_strategies/schemes/static_scheme.h"
#include "future/solving_strategies/strategies/linear_strategy.h"
#endif

namespace Kratos::Testing
{

namespace
{

class TestLaplacianElement1D : public Element
{
public:

    using IndexType = Element::IndexType;

    using DofsArrayType = Element::DofsArrayType;

    using DofsVectorType = Element::DofsVectorType;

    using EquationIdVectorType = Element::EquationIdVectorType;

    TestLaplacianElement1D(
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
        MatrixType &rLeftHandSideMatrix,
        const ProcessInfo &rCurrentProcessInfo) override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rLeftHandSideMatrix.size1() != n_nodes || rLeftHandSideMatrix.size2() != n_nodes) {
            rLeftHandSideMatrix.resize(n_nodes, n_nodes, false);
        }
        const double aux = 1.0 / r_geom.Length();
        rLeftHandSideMatrix(0,0) = aux;
        rLeftHandSideMatrix(0,1) = -aux;
        rLeftHandSideMatrix(1,0) = -aux;
        rLeftHandSideMatrix(1,1) = aux;
    }

    void CalculateRightHandSide(
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rRightHandSideVector.size() != n_nodes) {
            rRightHandSideVector.resize(n_nodes);
        }
        const double f = 1.0;
        const double h = r_geom.Length();
        const double aux_f = f * h / 2.0;
        const double aux_phi =  (r_geom[0].GetSolutionStepValue(DISTANCE) - r_geom[1].GetSolutionStepValue(DISTANCE)) / h;
        rRightHandSideVector(0) = aux_f - aux_phi;
        rRightHandSideVector(1) = aux_f + aux_phi;
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
        for (unsigned int i = 0; i < n_nodes; ++i) {
            rResult[i] = r_geom[i].GetDof(DISTANCE).EquationId();
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
        for (unsigned int i = 0; i < n_nodes; i++) {
            rElementalDofList[i] = r_geom[i].pGetDof(DISTANCE);
        }
    }

};

class TestElasticityElement2D : public Element
{
public:

    using IndexType = Element::IndexType;

    using DofsArrayType = Element::DofsArrayType;

    using DofsVectorType = Element::DofsVectorType;

    using EquationIdVectorType = Element::EquationIdVectorType;

    TestElasticityElement2D(
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
        MatrixType &rLeftHandSideMatrix,
        const ProcessInfo &rCurrentProcessInfo) override
    {
        const auto& r_geom = GetGeometry();
        const std::size_t local_size = 2 * r_geom.PointsNumber();
        if (rLeftHandSideMatrix.size1() != local_size || rLeftHandSideMatrix.size2() != local_size) {
            rLeftHandSideMatrix.resize(local_size, local_size, false);
        }

        const double C = E / (1.0 - std::pow(nu,2));
        const double a = r_geom[1].X0() - r_geom[0].X0();
        const double b = r_geom[3].Y0() - r_geom[0].Y0();
        const double k1 = C *((b/(3.0*a))*(1+nu) + (a/(3.0*b))*(1-nu));
        const double k2 = C *((b/(6.0*a))*(1+nu) - (a/(6.0*b))*(1-nu));
        const double k3 = C *((b/(6.0*a))*(1+nu) + (a/(6.0*b))*(1-nu));
        const double k4 = C *((b/(3.0*a))*(1+nu) - (a/(3.0*b))*(1-nu));
        const double k5 = C *((b/(4.0*a))*(1-nu) + (a/(4.0*b))*(1-nu));
        const double k6 = C *((b/(4.0*a))*(1-nu) - (a/(4.0*b))*(1-nu));

        rLeftHandSideMatrix(0,0) = k1;
        rLeftHandSideMatrix(1,1) = k5;
        rLeftHandSideMatrix(2,2) = k1;
        rLeftHandSideMatrix(3,3) = k5;
        rLeftHandSideMatrix(4,4) = k1;
        rLeftHandSideMatrix(5,5) = k5;
        rLeftHandSideMatrix(6,6) = k1;
        rLeftHandSideMatrix(7,7) = k5;

        rLeftHandSideMatrix(0,2) = k2;
        rLeftHandSideMatrix(0,4) = -k3;
        rLeftHandSideMatrix(0,6) = -k4;
        rLeftHandSideMatrix(2,0) = rLeftHandSideMatrix(0,2);
        rLeftHandSideMatrix(4,0) = rLeftHandSideMatrix(0,4);
        rLeftHandSideMatrix(6,0) = rLeftHandSideMatrix(0,6);

        rLeftHandSideMatrix(1,3) = k6;
        rLeftHandSideMatrix(1,5) = -k5;
        rLeftHandSideMatrix(1,7) = -k6;
        rLeftHandSideMatrix(3,1) = rLeftHandSideMatrix(1,3);
        rLeftHandSideMatrix(5,1) = rLeftHandSideMatrix(1,5);
        rLeftHandSideMatrix(7,1) = rLeftHandSideMatrix(1,7);

        rLeftHandSideMatrix(2,4) = -k4;
        rLeftHandSideMatrix(2,6) = -k3;
        rLeftHandSideMatrix(4,2) = rLeftHandSideMatrix(2,4);
        rLeftHandSideMatrix(6,2) = rLeftHandSideMatrix(2,6);

        rLeftHandSideMatrix(3,5) = -k6;
        rLeftHandSideMatrix(3,7) = -k5;
        rLeftHandSideMatrix(5,3) = rLeftHandSideMatrix(3,5);
        rLeftHandSideMatrix(7,3) = rLeftHandSideMatrix(3,7);

        rLeftHandSideMatrix(4,6) = k2;
        rLeftHandSideMatrix(6,4) = rLeftHandSideMatrix(4,6);

        rLeftHandSideMatrix(5,7) = k6;
        rLeftHandSideMatrix(7,5) = rLeftHandSideMatrix(5,7);
    }

    void CalculateRightHandSide(
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override
    {
        const auto &r_geom = GetGeometry();
        const std::size_t local_size = 2 * r_geom.PointsNumber();
        if (rRightHandSideVector.size() != local_size) {
            rRightHandSideVector.resize(local_size);
        }
        const double f_x = 0.0;
        const double f_y = -1.0;
        const double a = r_geom[1].X0() - r_geom[0].X0();
        const double b = r_geom[3].Y0() - r_geom[0].Y0();
        const double aux = a * b / 4.0;

        const double C = E / (1.0 - std::pow(nu,2));
        const double k1 = C *((b/(3.0*a))*(1+nu) + (a/(3.0*b))*(1-nu));
        const double k2 = C *((b/(6.0*a))*(1+nu) - (a/(6.0*b))*(1-nu));
        const double k3 = C *((b/(6.0*a))*(1+nu) + (a/(6.0*b))*(1-nu));
        const double k4 = C *((b/(3.0*a))*(1+nu) - (a/(3.0*b))*(1-nu));
        const double k5 = C *((b/(4.0*a))*(1-nu) + (a/(4.0*b))*(1-nu));
        const double k6 = C *((b/(4.0*a))*(1-nu) - (a/(4.0*b))*(1-nu));

        const double u1 = r_geom[0].GetSolutionStepValue(DISPLACEMENT_X);
        const double v1 = r_geom[0].GetSolutionStepValue(DISPLACEMENT_Y);
        const double u2 = r_geom[1].GetSolutionStepValue(DISPLACEMENT_X);
        const double v2 = r_geom[1].GetSolutionStepValue(DISPLACEMENT_Y);
        const double u3 = r_geom[2].GetSolutionStepValue(DISPLACEMENT_X);
        const double v3 = r_geom[2].GetSolutionStepValue(DISPLACEMENT_Y);
        const double u4 = r_geom[3].GetSolutionStepValue(DISPLACEMENT_X);
        const double v4 = r_geom[3].GetSolutionStepValue(DISPLACEMENT_Y);

        rRightHandSideVector(0) = aux * f_x - (k1*u1 + k2*u2 - k3*u3 - k4*u4);
        rRightHandSideVector(1) = aux * f_y - (k5*v1 + k6*v2 - k5*v3 - k6*v4);
        rRightHandSideVector(2) = aux * f_x - (k2*u1 + k1*u2 - k4*u3 - k3*u4);
        rRightHandSideVector(3) = aux * f_y - (k6*v1 + k5*v2 - k6*v3 - k5*v4);
        rRightHandSideVector(4) = aux * f_x - (-k3*u1 - k4*u2 + k1*u3 + k2*u4);
        rRightHandSideVector(5) = aux * f_y - (-k5*v1 - k6*v2 + k5*v3 + k6*v4);
        rRightHandSideVector(6) = aux * f_x - (-k4*u1 - k3*u2 + k2*u3 + k1*u4);
        rRightHandSideVector(7) = aux * f_y - (-k6*v1 - k5*v2 + k6*v3 + k5*v4);
    }

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override
    {
        const auto& r_geom = GetGeometry();
        const std::size_t n_nodes = r_geom.PointsNumber();
        const std::size_t local_size = 2 * n_nodes;
        if (rResult.size() != local_size) {
            rResult.resize(local_size);
        }
        for (unsigned int i = 0; i < n_nodes; ++i) {
            rResult[2*i] = r_geom[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[2*i + 1] = r_geom[i].GetDof(DISPLACEMENT_Y).EquationId();
        }
    }

    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        const std::size_t local_size = 2 * n_nodes;
        if (rElementalDofList.size() != local_size) {
            rElementalDofList.resize(local_size);
        }
        for (unsigned int i = 0; i < n_nodes; i++) {
            rElementalDofList[2*i] = r_geom[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[2*i + 1] = r_geom[i].pGetDof(DISPLACEMENT_Y);
        }
    }

private:

    static constexpr double E = 1e6;

    static constexpr double nu = 0.3;

};

static void SetUpTestModelPart1D(
    const std::size_t NumElems,
    const double ElemSize,
    ModelPart& rModelPart)
{
    const int domain_size = 2;
    const std::size_t buffer_size = 3;
    rModelPart.SetBufferSize(buffer_size);
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, domain_size);

    rModelPart.AddNodalSolutionStepVariable(DISTANCE);

    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    for (std::size_t i_elem = 0; i_elem < NumElems; ++i_elem) {
        auto p_node_a = rModelPart.pGetNode(i_elem + 1);
        auto p_node_b = rModelPart.CreateNewNode(i_elem + 2, (i_elem + 1)*ElemSize, 0.0, 0.0);
        auto p_geom = Kratos::make_shared<Line2D2<Node>>(p_node_a, p_node_b);
        rModelPart.AddElement(Kratos::make_intrusive<TestLaplacianElement1D>(i_elem + 1, p_geom));
    }

    rModelPart.GetNodalSolutionStepVariablesList().AddDof(&DISTANCE);
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(DISTANCE);
    }

    for (auto& r_elem : rModelPart.Elements()) {
        r_elem.Set(ACTIVE, true);
    }
}

static void SetUpTestModelPart2D(
    const std::size_t NumElemsX,
    const std::size_t NumElemsY,
    const double ElemSizeX,
    const double ElemSizeY,
    ModelPart& rModelPart)
{
    const int domain_size = 2;
    const std::size_t buffer_size = 3;
    rModelPart.SetBufferSize(buffer_size);
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, domain_size);

    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);

    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    for (std::size_t i = 0; i < NumElemsY; ++i) {
        for (std::size_t j = 0; j < NumElemsX; ++j) {
            auto p_node_a = rModelPart.pGetNode(i * (NumElemsX + 1) + j + 1);
            auto p_node_b = rModelPart.CreateNewNode(i * (NumElemsX + 1) + (j + 1) + 1, p_node_a->X0() + ElemSizeX, p_node_a->Y0(), 0.0);
            auto p_node_c = rModelPart.CreateNewNode((i + 1) * (NumElemsX + 1) + (j + 1) + 1, p_node_a->X0() + ElemSizeX, p_node_a->Y0()+ ElemSizeY, 0.0);
            auto p_node_d = j == 0 ?
                rModelPart.CreateNewNode((i + 1) * (NumElemsX + 1) + j + 1, p_node_a->X0(), p_node_a->Y0() + ElemSizeY, 0.0) :
                rModelPart.pGetNode((i + 1) * (NumElemsX + 1) + j + 1);
            auto p_geom = Kratos::make_shared<Quadrilateral2D4<Node>>(p_node_a, p_node_b, p_node_c, p_node_d);
            rModelPart.AddElement(Kratos::make_intrusive<TestElasticityElement2D>(i * NumElemsX + j + 1, p_geom));
        }
    }

    rModelPart.GetNodalSolutionStepVariablesList().AddDof(&DISPLACEMENT_X);
    rModelPart.GetNodalSolutionStepVariablesList().AddDof(&DISPLACEMENT_Y);
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
    }

    for (auto& r_elem : rModelPart.Elements()) {
        r_elem.Set(ACTIVE, true);
    }
}

}

KRATOS_TEST_CASE_IN_SUITE(StaticSchemeBuild1D, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 2;
    const double elem_size = 1.0;
    SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "block_builder"
        }
    })");
    using SchemeType = Future::StaticScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    auto p_scheme = Kratos::make_unique<SchemeType>(r_test_model_part, scheme_settings);

    // Create the DOF set and set the global ids
    // Note that in a standard case this happens at the strategy level
    auto p_dof_set = Kratos::make_shared<ModelPart::DofsArrayType>();
    auto p_eff_dof_set = Kratos::make_shared<ModelPart::DofsArrayType>();

    // Set up the matrix graph and arrays
    // Note that in a standard case this happens at the strategy level
    Future::LinearSystemContainer<CsrMatrix<>, SystemVector<>> linear_system_container;

    // Call the initialize solution step (note that this sets all the arrays above)
    p_scheme->InitializeSolutionStep(p_dof_set, p_eff_dof_set, linear_system_container);

    // Call the build
    auto p_lhs = linear_system_container.pLhs;
    auto p_rhs = linear_system_container.pRhs;
    p_scheme->Build(*p_lhs, *p_rhs);

    // Check resultant matrices
    const double tol = 1.0e-12;
    std::vector<double> expected_rhs = {0.5,1.0,0.5};
    BoundedMatrix<double,3,3> expected_lhs;
    expected_lhs(0,0) = 1.0; expected_lhs(0,1) = -1.0; expected_lhs(0,2) = 0.0;
    expected_lhs(1,0) = -1.0; expected_lhs(1,1) = 2.0; expected_lhs(1,2) = -1.0;
    expected_lhs(2,0) = 0.0; expected_lhs(2,1) = -1.0; expected_lhs(2,2) = 1.0;
    KRATOS_CHECK_VECTOR_NEAR((*p_rhs), expected_rhs, tol); // Note that as there are not non-zero entries in the sparse vector we can use the standard macro
    for (unsigned int i = 0; i < p_lhs->size1(); ++i) {
        for (unsigned int j = 0; j < p_lhs->size2(); ++j) {
            const double expected_val = expected_lhs(i,j);
            if (std::abs(expected_val) > tol) {
                KRATOS_CHECK_NEAR(p_lhs->operator()(i,j), expected_val, tol); // Note that we check if the expected value is non-zero as this is a CSR matrix
            }
        }
    }
#else
    true;
#endif
}

KRATOS_TEST_CASE_IN_SUITE(StaticSchemeBuild2D, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems_x = 2500;
    const std::size_t num_elems_y = 2500;
    const double elem_size_x = 1.0;
    const double elem_size_y = 1.0;

    BuiltinTimer timer_mesh;
    SetUpTestModelPart2D(num_elems_x, num_elems_y, elem_size_x, elem_size_y, r_test_model_part);
    std::cout << "Mesh generation time: " << timer_mesh << std::endl;
    std::cout << "Number of elements: " << num_elems_x * num_elems_y << std::endl;

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "echo_level" : 1,
        "build_settings" : {
            "name" : "block_builder"
        }
    })");
    using SchemeType = Future::StaticScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    auto p_scheme = Kratos::make_unique<SchemeType>(r_test_model_part, scheme_settings);

    // Create the DOF set and set the global ids
    // Note that in a standard case this happens at the strategy level
    auto p_dof_set = Kratos::make_shared<ModelPart::DofsArrayType>();
    auto p_eff_dof_set = Kratos::make_shared<ModelPart::DofsArrayType>();

    // Set up the matrix graph and arrays
    // Note that in a standard case this happens at the strategy level
    Future::LinearSystemContainer<CsrMatrix<>, SystemVector<>> linear_system_container;

    // Call the initialize solution step (note that this sets all the arrays above)
    p_scheme->InitializeSolutionStep(p_dof_set, p_eff_dof_set, linear_system_container);

    // Call the build
    auto p_lhs = linear_system_container.pLhs;
    auto p_rhs = linear_system_container.pRhs;

    sleep(30);

    BuiltinTimer timer_build;
    for(unsigned int i=0; i<20; ++i)
        p_scheme->Build(*p_lhs, *p_rhs);
    std::cout << "Build time: " << timer_build << std::endl;


    p_lhs->SetValue(0.0);
    p_rhs->SetValue(0.0);

    sleep(30);

    BuiltinTimer timer_build_safe;
    for(unsigned int i=0; i<20; ++i)
    p_scheme->BuildWithSafeAssemble(*p_lhs, *p_rhs);
    std::cout << "Build w/ safe assemble time: " << timer_build_safe << std::endl;

    p_lhs->SetValue(0.0);
    p_rhs->SetValue(0.0);

    sleep(30);

    BuiltinTimer timer_build_thread_local;
    for(unsigned int i=0; i<20; ++i)
    p_scheme->BuildWithThreadLocal(*p_lhs, *p_rhs);
    std::cout << "Build w/ thread local: " << timer_build_thread_local << std::endl;

    p_lhs->SetValue(0.0);
    p_rhs->SetValue(0.0);

    sleep(30);

    BuiltinTimer timer_build_local_allocation;
    for(unsigned int i=0; i<20; ++i)
        p_scheme->BuildWithLocalAllocation(*p_lhs, *p_rhs);
    std::cout << "Build w/ local allocation: " << timer_build_local_allocation << std::endl;

    // // Check resultant matrices
    // const double tol = 1.0e-12;
    // std::vector<double> expected_rhs = {0.5,1.0,0.5};
    // BoundedMatrix<double,3,3> expected_lhs;
    // expected_lhs(0,0) = 1.0; expected_lhs(0,1) = -1.0; expected_lhs(0,2) = 0.0;
    // expected_lhs(1,0) = -1.0; expected_lhs(1,1) = 2.0; expected_lhs(1,2) = -1.0;
    // expected_lhs(2,0) = 0.0; expected_lhs(2,1) = -1.0; expected_lhs(2,2) = 1.0;
    // KRATOS_CHECK_VECTOR_NEAR((*p_rhs), expected_rhs, tol); // Note that as there are not non-zero entries in the sparse vector we can use the standard macro
    // for (unsigned int i = 0; i < p_lhs->size1(); ++i) {
    //     for (unsigned int j = 0; j < p_lhs->size2(); ++j) {
    //         const double expected_val = expected_lhs(i,j);
    //         if (std::abs(expected_val) > tol) {
    //             KRATOS_CHECK_NEAR(p_lhs->operator()(i,j), expected_val, tol); // Note that we check if the expected value is non-zero as this is a CSR matrix
    //         }
    //     }
    // }
#else
    true;
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategy, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 2;
    const double elem_size = 1.0;
    SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "block_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<CsrMatrix<>, SystemVector<>>;
    using LinearSolverType = Future::LinearSolver<CsrMatrix<>, SystemVector<>>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs
    auto p_node_1 = r_test_model_part.pGetNode(1);
    p_node_1->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_1->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->InitializeSolutionStep();
    p_strategy->Predict();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();
    p_strategy->Clear();

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 2.5, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 3.0, 1.0e-12);
#else
    true;
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyWithJumpConstraint, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 2;
    const double elem_size = 1.0;
    SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Mesh:  x ---- x ---- x
    // DOFs: u_1 -- u_2 -- u_3
    // Create a constraint with jump to impose u_3 = u_1 + 2
    // Note that as u_1 is fixed to 1.0 this actually imposes the Laplacian solution
    // This checks the MPC with fixed DOF but ensuring a non-overconstrained and compatible system of equations
    const double jump = 2.0;
    const double weight = 1.0;
    auto& r_master_node = r_test_model_part.GetNode(1);
    auto& r_slave_node = r_test_model_part.GetNode(3);
    auto p_const_1 = r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, r_master_node, DISTANCE, r_slave_node, DISTANCE, weight, jump);
    p_const_1->Set(ACTIVE, true);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "block_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<CsrMatrix<>, SystemVector<>>;
    using LinearSolverType = Future::LinearSolver<CsrMatrix<>, SystemVector<>>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs
    auto p_node_1 = r_test_model_part.pGetNode(1);
    p_node_1->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_1->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->InitializeSolutionStep();
    p_strategy->Predict();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();
    p_strategy->Clear();

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 2.5, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 3.0, 1.0e-12);
#else
    true;
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyWithPeriodicityConstraint, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 3;
    const double elem_size = 1.0;
    SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Mesh:  x ---- x ---- x ---- x
    // DOFs: u_1 -- u_2 -- u_3 -- u_4
    // Create a periodicity constraint with jump to impose u_4 = u_1 + 1
    const double jump = 1.0;
    const double weight = 1.0;
    auto& r_master_node = r_test_model_part.GetNode(1);
    auto& r_slave_node = r_test_model_part.GetNode(4);
    auto p_const_1 = r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, r_master_node, DISTANCE, r_slave_node, DISTANCE, weight, jump);
    p_const_1->Set(ACTIVE, true);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "block_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<CsrMatrix<>, SystemVector<>>;
    using LinearSolverType = Future::LinearSolver<CsrMatrix<>, SystemVector<>>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs
    auto p_node_1 = r_test_model_part.pGetNode(1);
    p_node_1->FastGetSolutionStepValue(DISTANCE, 0) = 0.0;
    p_node_1->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->InitializeSolutionStep();
    p_strategy->Predict();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();
    p_strategy->Clear();

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 0.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 4.0/3.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 5.0/3.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(4).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
#else
    true;
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyWithMultipleDofsConstraints, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 3;
    const double elem_size = 1.0;
    SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Mesh:  x ---- x ---- x ---- x
    // DOFs: u_1 -- u_2 -- u_3 -- u_4
    // Create a constraint involving multiple master DOFs such that u_4 = 0.5*u_1 + 0.5*u_2 + 1.0
    std::vector<typename Dof<double>::Pointer> slave_dofs(1);
    std::vector<typename Dof<double>::Pointer> master_dofs(2);
    slave_dofs[0] = r_test_model_part.GetNode(4).pGetDof(DISTANCE);
    master_dofs[0] = r_test_model_part.GetNode(1).pGetDof(DISTANCE);
    master_dofs[1] = r_test_model_part.GetNode(2).pGetDof(DISTANCE);

    Vector jump_vector(1);
    jump_vector[0] = 1.0;

    Matrix relation_matrix(1, 2);
    relation_matrix(0,0) = 0.5;
    relation_matrix(0,1) = 0.5;

    auto p_const_1 = r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, master_dofs, slave_dofs, relation_matrix, jump_vector);
    p_const_1->Set(ACTIVE, true);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "block_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<CsrMatrix<>, SystemVector<>>;
    using LinearSolverType = Future::LinearSolver<CsrMatrix<>, SystemVector<>>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs
    auto p_node_1 = r_test_model_part.pGetNode(1);
    p_node_1->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_1->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->InitializeSolutionStep();
    p_strategy->Predict();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();
    p_strategy->Clear();

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 3.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 3.5, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(4).FastGetSolutionStepValue(DISTANCE), 3.0, 1.0e-12);
#else
    true;
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyWithTieConstraints, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 3;
    const double elem_size = 1.0;
    SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Create the nodes to apply the tie constraints
    auto p_node_5 = r_test_model_part.CreateNewNode(5, 0.0, 0.0, 0.0);
    auto p_node_6 = r_test_model_part.CreateNewNode(6, 0.0, 0.0, 0.0);
    p_node_5->AddDof(DISTANCE);
    p_node_6->AddDof(DISTANCE);

    // Mesh:  x    x ---- x ---- x ---- x    x
    // DOFs: u_5  u_1 -- u_2 -- u_3 -- u_4  u_6
    // Create a tie constraint at each side of the domain s.t. u_1 = u_5 and u_4 = u_6
    // Note that these nodes are not connected to the domain (i.e., are somehow "flying")
    const double jump = 0.0;
    const double weight = 1.0;

    auto& r_slave_node_a = r_test_model_part.GetNode(1);
    auto& r_master_node_a = r_test_model_part.GetNode(5);
    auto p_const_1 = r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, r_master_node_a, DISTANCE, r_slave_node_a, DISTANCE, weight, jump);
    p_const_1->Set(ACTIVE, true);

    auto& r_slave_node_b = r_test_model_part.GetNode(4);
    auto& r_master_node_b = r_test_model_part.GetNode(6);
    auto p_const_2 = r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 2, r_master_node_b, DISTANCE, r_slave_node_b, DISTANCE, weight, jump);
    p_const_2->Set(ACTIVE, true);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "block_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<CsrMatrix<>, SystemVector<>>;
    using LinearSolverType = Future::LinearSolver<CsrMatrix<>, SystemVector<>>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs to the tying ("flying") nodes
    p_node_5->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_5->Fix(DISTANCE);

    p_node_6->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_6->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->InitializeSolutionStep();
    p_strategy->Predict();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();
    p_strategy->Clear();

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 2.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 2.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(4).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
#else
    true;
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyWithRigidBodyMotionConstraint, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 3;
    const double elem_size = 1.0;
    SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Mesh:  x ---- x ---- x ---- x
    // DOFs: u_1 -- u_2 -- u_3 -- u_4
    // Create an average constraint to enforce the average solution to be zero (i.e., no rigid body motion)
    // We do this as u_4 = - u_1 - u_2 - u_3 as u_1 + u_2 + u_3 + u_4 = u_1 + u_2 + u_3 + (- u_1 - u_2 - u_3) = 0
    typename LinearMasterSlaveConstraint::DofPointerVectorType slave_dofs(1);
    slave_dofs[0] = r_test_model_part.GetNode(4).pGetDof(DISTANCE);
    typename LinearMasterSlaveConstraint::DofPointerVectorType master_dofs(3);
    master_dofs[0] = r_test_model_part.GetNode(1).pGetDof(DISTANCE);
    master_dofs[1] = r_test_model_part.GetNode(2).pGetDof(DISTANCE);
    master_dofs[2] = r_test_model_part.GetNode(3).pGetDof(DISTANCE);
    Matrix relation_matrix(1,3);
    relation_matrix(0,0) = -1.0;
    relation_matrix(0,1) = -1.0;
    relation_matrix(0,2) = -1.0;
    Vector constant_vector = ZeroVector(1);
    auto p_const_1 = r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, master_dofs, slave_dofs, relation_matrix, constant_vector);
    p_const_1->Set(ACTIVE, true);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "block_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<CsrMatrix<>, SystemVector<>>;
    using LinearSolverType = Future::LinearSolver<CsrMatrix<>, SystemVector<>>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->InitializeSolutionStep();
    p_strategy->Predict();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();
    p_strategy->Clear();

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), -0.125, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 0.125, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 0.125, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(4).FastGetSolutionStepValue(DISTANCE), -0.125, 1.0e-12);
#else
    true;
#endif
}

}  // namespace Kratos::Testing.

