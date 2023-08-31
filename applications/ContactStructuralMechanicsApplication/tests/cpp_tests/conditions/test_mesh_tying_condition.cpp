// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "contact_structural_mechanics_application_variables.h"
#include "custom_conditions/mesh_tying_mortar_condition.h"
#include "tests/test_utilities/test_laplacian_element.h"
// #include "input_output/vtk_output.h"

/* Geometries */
#include "geometries/line_2d_2.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"

/* Linear solvers */
#include "spaces/ublas_space.h"
#include "linear_solvers/reorderer.h"
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
//#include "linear_solvers/amgcl_solver.h"

/* Strategies */
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
//#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

/* The most basic scheme (static) */
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

/* The builder and solvers */
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

namespace Kratos::Testing 
{

/// Initial definitons
using GeometryType = Geometry<Node>;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

/// The direct solver
using ReordererType = Reorderer<SparseSpaceType,  LocalSpaceType>;
using DirectSolverType = DirectSolver<SparseSpaceType,  LocalSpaceType, ReordererType>;
using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;
using SkylineLUFactorizationSolverType = SkylineLUFactorizationSolver<SparseSpaceType,  LocalSpaceType, ReordererType>;
// using AMGCLSolverType = AMGCLSolver<SparseSpaceType, LocalSpaceType, ReordererType>;

/// The builder and solver type
using BuilderAndSolverType = BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>;
using ResidualBasedBlockBuilderAndSolverType = ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>;

/// The time scheme
using SchemeType = Scheme<SparseSpaceType, LocalSpaceType>;
using ResidualBasedIncrementalUpdateStaticSchemeType = ResidualBasedIncrementalUpdateStaticScheme<SparseSpaceType, LocalSpaceType>;

// Strategies
using ResidualBasedLinearStrategyType = ResidualBasedLinearStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;
//using ResidualBasedNewtonRaphsonStrategyType = ResidualBasedNewtonRaphsonStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;

// /**
//  * @brief Generates a function comment for the given function body.
//  * @param rModelPart the model part to perform the function on
//  */
// void IOMeshTyingDebug(ModelPart& rModelPart)
// {
//     Parameters io_parameters = Parameters(R"(
//     {
//         "file_format"                        : "ascii",
//         "nodal_solution_step_data_variables" : ["TEMPERATURE"]
//     })" );
//     VtkOutput(rModelPart, io_parameters).PrintOutput();
// }

/**
 * @brief Set boundary conditions for the given model part.
 * @param rModelPart the model part for which to set the boundary conditions
 */
void SetBC(ModelPart& rModelPart)
{
    // Set initial solution
    for (auto& r_node : rModelPart.Nodes()) {
        if (r_node.X() < 0.01) {
            r_node.FastGetSolutionStepValue(TEMPERATURE) = 1.0;
            r_node.FastGetSolutionStepValue(TEMPERATURE, 1) = 1.0;
            r_node.FastGetSolutionStepValue(TEMPERATURE, 2) = 1.0;
        } else {
            r_node.FastGetSolutionStepValue(TEMPERATURE) = 0.0;
            r_node.FastGetSolutionStepValue(TEMPERATURE, 1) = 0.0;
            r_node.FastGetSolutionStepValue(TEMPERATURE, 2) = 0.0;
        }
    }

    // Fix dofs
    for (auto& r_node : rModelPart.Nodes()) {
        if (r_node.X() < 0.01 || r_node.X() > 0.99) {
            r_node.Fix(TEMPERATURE);
        }
    }
}

/**
 * @brief Check the solutions of the given ModelPart.
 * @param rModelPart The ModelPart to check the solutions of.
 */
void CheckSolution(ModelPart& rModelPart)
{
    // Check solutions
    const double tolerance = 1.0e-6;
    for (auto& r_node : rModelPart.Nodes()) {
        if (r_node.X() < 0.01) {
            KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(TEMPERATURE), 1.0, tolerance);
        } else if (r_node.X() > 0.99) {
            KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(TEMPERATURE), 0.0, tolerance);
        } else {
            KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(TEMPERATURE), 0.5, tolerance);
        }
    }
}

// /**
//  * @brief Generates a the simplest reference model part.
//  * @param rModelPart the model part to generate the reference model for
//  */
// void GenerateReferenceSimplestModelPart(ModelPart& rModelPart)
// {
//     // Adding variable 
//     rModelPart.AddNodalSolutionStepVariable(NORMAL);
//     rModelPart.AddNodalSolutionStepVariable(TEMPERATURE);
//     rModelPart.AddNodalSolutionStepVariable(REACTION_FLUX);

//     // Creating properties
//     auto p_prop = rModelPart.CreateNewProperties(1, 0);

//     // Set conductivity
//     p_prop->SetValue(CONDUCTIVITY, 1.0);

//     // Creating nodes 
//     auto pnode1 = rModelPart.CreateNewNode(1, 0.0, 0.5, 0.0);
//     auto pnode2 = rModelPart.CreateNewNode(2, 0.5, 0.0, 0.0);
//     auto pnode3 = rModelPart.CreateNewNode(3, 0.5, 1.0, 0.0);
//     auto pnode4 = rModelPart.CreateNewNode(4, 1.0, 0.5, 0.0);

//     /// Add dof
//     for (auto& r_node : rModelPart.Nodes()) {
//         r_node.AddDof(TEMPERATURE, REACTION_FLUX);
//     }

//     // Creating elements
//     GeometryType::Pointer p_elem_geom1 = Kratos::make_shared<Triangle2D3<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode1, pnode2, pnode3})});
//     rModelPart.AddElement(Kratos::make_intrusive<TestLaplacianElement>( 1, p_elem_geom1, p_prop));
//     GeometryType::Pointer p_elem_geom2 = Kratos::make_shared<Triangle2D3<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode2, pnode4, pnode3})});
//     rModelPart.AddElement(Kratos::make_intrusive<TestLaplacianElement>( 2, p_elem_geom2, p_prop));

//     // Set BC
//     SetBC(rModelPart);
// }

// /**
//  * @brief Generates a reference model part.
//  * @param rModelPart the model part to generate the reference model for
//  */
// void GenerateReferenceModelPart(ModelPart& rModelPart)
// {
//     // Adding variable
//     rModelPart.AddNodalSolutionStepVariable(NORMAL);
//     rModelPart.AddNodalSolutionStepVariable(TEMPERATURE);
//     rModelPart.AddNodalSolutionStepVariable(REACTION_FLUX);

//     // Creating properties
//     auto p_prop = rModelPart.CreateNewProperties(1, 0);

//     // Set conductivity
//     p_prop->SetValue(CONDUCTIVITY, 1.0);

//     // Creating nodes 
//     auto pnode1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
//     auto pnode2 = rModelPart.CreateNewNode(2, 0.5, 0.0, 0.0);
//     auto pnode3 = rModelPart.CreateNewNode(3, 0.5, 0.5, 0.0);
//     auto pnode4 = rModelPart.CreateNewNode(4, 0.0, 0.5, 0.0);
//     auto pnode5 = rModelPart.CreateNewNode(5, 0.5, 1.0, 0.0);
//     auto pnode6 = rModelPart.CreateNewNode(6, 0.0, 1.0, 0.0);
//     auto pnode7 = rModelPart.CreateNewNode(7, 1.0, 0.0, 0.0);
//     auto pnode8 = rModelPart.CreateNewNode(8, 1.0, 0.5, 0.0);
//     auto pnode9 = rModelPart.CreateNewNode(9, 1.0, 1.0, 0.0);

//     /// Add dof
//     for (auto& r_node : rModelPart.Nodes()) {
//         r_node.AddDof(TEMPERATURE, REACTION_FLUX);
//     }

//     // Creating elements
//     GeometryType::Pointer p_elem_geom1 = Kratos::make_shared<Quadrilateral2D4<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode1, pnode2, pnode3, pnode4})});
//     rModelPart.AddElement(Kratos::make_intrusive<TestLaplacianElement>( 1, p_elem_geom1, p_prop));
//     GeometryType::Pointer p_elem_geom2 = Kratos::make_shared<Quadrilateral2D4<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode4, pnode3, pnode5, pnode6})});
//     rModelPart.AddElement(Kratos::make_intrusive<TestLaplacianElement>( 2, p_elem_geom2, p_prop));
//     GeometryType::Pointer p_elem_geom3 = Kratos::make_shared<Quadrilateral2D4<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode2, pnode7, pnode8, pnode3})});
//     rModelPart.AddElement(Kratos::make_intrusive<TestLaplacianElement>( 3, p_elem_geom3, p_prop));
//     GeometryType::Pointer p_elem_geom4 = Kratos::make_shared<Quadrilateral2D4<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode8, pnode9, pnode5, pnode3})});
//     rModelPart.AddElement(Kratos::make_intrusive<TestLaplacianElement>( 4, p_elem_geom4, p_prop));

//     // Set BC
//     SetBC(rModelPart);
// }

/**
 * @brief Generates a simplest mesh tying model part.
 * @param rModelPart the model part to generate the mesh tying model for
 */
void GenerateMeshTyingSimplestModelPart(ModelPart& rModelPart)
{
    // Adding variable
    rModelPart.AddNodalSolutionStepVariable(NORMAL);
    rModelPart.AddNodalSolutionStepVariable(TEMPERATURE);
    rModelPart.AddNodalSolutionStepVariable(REACTION_FLUX);
    rModelPart.AddNodalSolutionStepVariable(SCALAR_LAGRANGE_MULTIPLIER);

    // Creating properties
    auto p_prop = rModelPart.CreateNewProperties(1, 0);

    // Set conductivity
    p_prop->SetValue(CONDUCTIVITY, 1.0);
    p_prop->SetValue(TYING_VARIABLE, "TEMPERATURE");

    // Creating nodes 
    // First we create first side nodes
    auto pnode1 = rModelPart.CreateNewNode(1, 0.0, 0.5, 0.0);
    auto pnode2 = rModelPart.CreateNewNode(2, 0.5, 0.0, 0.0);
    auto pnode3 = rModelPart.CreateNewNode(3, 0.5, 1.0, 0.0);

    // Now we create second side nodes
    auto pnode4 = rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);
    auto pnode5 = rModelPart.CreateNewNode(5, 0.5, 1.0, 0.0);
    auto pnode6 = rModelPart.CreateNewNode(6, 1.0, 0.5, 0.0);

    /// Add dof
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(TEMPERATURE, REACTION_FLUX);
        r_node.AddDof(SCALAR_LAGRANGE_MULTIPLIER);
    }

    // Creating elements
    GeometryType::Pointer p_elem_geom1 = Kratos::make_shared<Triangle2D3<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode1, pnode2, pnode3})});
    rModelPart.AddElement(Kratos::make_intrusive<TestLaplacianElement>( 1, p_elem_geom1, p_prop));
    GeometryType::Pointer p_elem_geom2 = Kratos::make_shared<Triangle2D3<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode5, pnode4, pnode6})});
    rModelPart.AddElement(Kratos::make_intrusive<TestLaplacianElement>( 2, p_elem_geom2, p_prop));

    // Conditions geometries
    GeometryType::Pointer p_cond_geom1 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode2, pnode3})});
    GeometryType::Pointer p_cond_geom2 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode5, pnode4})});

    // Calculate normal
    const array_1d<double, 3> slave_normal = p_cond_geom1->UnitNormal(p_cond_geom1->Center());
    for (auto& r_node : p_cond_geom1->Points()) {
        r_node.FastGetSolutionStepValue(NORMAL) = slave_normal;
    }
    const array_1d<double, 3> master_normal = p_cond_geom2->UnitNormal(p_cond_geom2->Center());
    for (auto& r_node : p_cond_geom2->Points()) {
        r_node.FastGetSolutionStepValue(NORMAL) = master_normal;
    }

    // Creating conditions
    auto& r_prototype = dynamic_cast<const PairedCondition&>(KratosComponents<Condition>::Get("MeshTyingMortarCondition2D2N"));
    Condition::Pointer p_cond_1 = r_prototype.Create(1, p_cond_geom1, p_prop, p_cond_geom2);
    p_cond_1->SetValue(NORMAL, slave_normal);
    rModelPart.AddCondition(p_cond_1);

    // Set BC
    SetBC(rModelPart);
}

/**
 * @brief Generates a mesh tying model part.
 * @param rModelPart the model part to generate the mesh tying model for
 */
void GenerateMeshTyingModelPart(ModelPart& rModelPart)
{
    // Adding variable
    rModelPart.AddNodalSolutionStepVariable(NORMAL);
    rModelPart.AddNodalSolutionStepVariable(TEMPERATURE);
    rModelPart.AddNodalSolutionStepVariable(REACTION_FLUX);
    rModelPart.AddNodalSolutionStepVariable(SCALAR_LAGRANGE_MULTIPLIER);

    // Creating properties
    auto p_prop = rModelPart.CreateNewProperties(1, 0);

    // Set conductivity
    p_prop->SetValue(CONDUCTIVITY, 1.0);
    p_prop->SetValue(TYING_VARIABLE, "TEMPERATURE");

    // Creating nodes 
    // First we create first side nodes
    auto pnode1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto pnode2 = rModelPart.CreateNewNode(2, 0.5, 0.0, 0.0);
    auto pnode3 = rModelPart.CreateNewNode(3, 0.5, 0.5, 0.0);
    auto pnode4 = rModelPart.CreateNewNode(4, 0.0, 0.5, 0.0);
    auto pnode5 = rModelPart.CreateNewNode(5, 0.5, 1.0, 0.0);
    auto pnode6 = rModelPart.CreateNewNode(6, 0.0, 1.0, 0.0);

    // Now we create second side nodes
    auto pnode7 = rModelPart.CreateNewNode(7, 0.5, 0.0, 0.0);
    auto pnode8 = rModelPart.CreateNewNode(8, 1.0, 0.0, 0.0);
    auto pnode9 = rModelPart.CreateNewNode(9, 1.0, 0.5, 0.0);
    auto pnode10 = rModelPart.CreateNewNode(10, 0.5, 0.5, 0.0);
    auto pnode11 = rModelPart.CreateNewNode(11, 1.0, 1.0, 0.0);
    auto pnode12 = rModelPart.CreateNewNode(12, 0.5, 1.0, 0.0);

    /// Add dof
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(TEMPERATURE, REACTION_FLUX);
        r_node.AddDof(SCALAR_LAGRANGE_MULTIPLIER);
    }

    // Creating elements
    GeometryType::Pointer p_elem_geom1 = Kratos::make_shared<Quadrilateral2D4<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode1, pnode2, pnode3, pnode4})});
    rModelPart.AddElement(Kratos::make_intrusive<TestLaplacianElement>( 1, p_elem_geom1, p_prop));
    GeometryType::Pointer p_elem_geom2 = Kratos::make_shared<Quadrilateral2D4<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode4, pnode3, pnode5, pnode6})});
    rModelPart.AddElement(Kratos::make_intrusive<TestLaplacianElement>( 2, p_elem_geom2, p_prop));

    GeometryType::Pointer p_elem_geom3 = Kratos::make_shared<Quadrilateral2D4<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode7, pnode8, pnode9, pnode10})});
    rModelPart.AddElement(Kratos::make_intrusive<TestLaplacianElement>( 3, p_elem_geom3, p_prop));
    GeometryType::Pointer p_elem_geom4 = Kratos::make_shared<Quadrilateral2D4<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode10, pnode9, pnode11, pnode12})});
    rModelPart.AddElement(Kratos::make_intrusive<TestLaplacianElement>( 4, p_elem_geom4, p_prop));

    // Conditions geometries
    GeometryType::Pointer p_cond_geom1 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode2, pnode3})});
    GeometryType::Pointer p_cond_geom2 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode3, pnode5})});
    GeometryType::Pointer p_cond_geom3 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode10, pnode7})});
    GeometryType::Pointer p_cond_geom4 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode12, pnode10})});

    // Calculate normal
    const array_1d<double, 3> slave_normal = p_cond_geom1->UnitNormal(p_cond_geom1->Center());
    for (auto& r_node : p_cond_geom1->Points()) {
        r_node.FastGetSolutionStepValue(NORMAL) = slave_normal;
    }
    for (auto& r_node : p_cond_geom2->Points()) {
        r_node.FastGetSolutionStepValue(NORMAL) = slave_normal;
    }
    const array_1d<double, 3> master_normal = p_cond_geom3->UnitNormal(p_cond_geom3->Center());
    for (auto& r_node : p_cond_geom3->Points()) {
        r_node.FastGetSolutionStepValue(NORMAL) = master_normal;
    }
    for (auto& r_node : p_cond_geom4->Points()) {
        r_node.FastGetSolutionStepValue(NORMAL) = master_normal;
    }

    // Creating conditions
    auto& r_prototype = dynamic_cast<const PairedCondition&>(KratosComponents<Condition>::Get("MeshTyingMortarCondition2D2N"));
    Condition::Pointer p_cond_1 = r_prototype.Create(1, p_cond_geom1, p_prop, p_cond_geom3);
    p_cond_1->SetValue(NORMAL, slave_normal);
    rModelPart.AddCondition(p_cond_1);
    Condition::Pointer p_cond_2 = r_prototype.Create(2, p_cond_geom2, p_prop, p_cond_geom4);
    p_cond_2->SetValue(NORMAL, slave_normal);
    rModelPart.AddCondition(p_cond_2);

    // Set BC
    SetBC(rModelPart);
}

/**
 * @brief Builds the system and solves it.
 * @param rModelPart The model part to be solved.
 */
static void SolveSystem(ModelPart& rModelPart)
{
    // Build and solve system
    SchemeType::Pointer p_scheme = SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
    LinearSolverType::Pointer p_solver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
    //LinearSolverType::Pointer p_solver = LinearSolverType::Pointer( new AMGCLSolverType() );
    Parameters parameters = Parameters(R"(
    {
        "diagonal_values_for_dirichlet_dofs" : "no_scaling",
        "silent_warnings"                    : false
    })" );
    BuilderAndSolverType::Pointer p_builder_and_solver = BuilderAndSolverType::Pointer( new ResidualBasedBlockBuilderAndSolverType(p_solver, parameters) );

    // Creating linear strategy
    auto linear_solver = ResidualBasedLinearStrategyType(rModelPart, p_scheme, p_builder_and_solver);
    
    // Solving the system
    linear_solver.Solve();
}

/** 
* Checks the correct work of the mesh tying condition (simplest case)
*/
KRATOS_TEST_CASE_IN_SUITE(MeshTyingCondition1, KratosContactStructuralMechanicsFastSuite)
{
    // Create model part
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 3);

    // Fill test model part
    //GenerateReferenceSimplestModelPart(r_model_part);
    GenerateMeshTyingSimplestModelPart(r_model_part);

    // Solve system
    SolveSystem(r_model_part);

    // // Check results
    // IOMeshTyingDebug(r_model_part);

    // Check solutions
    CheckSolution(r_model_part);
}

/** 
* Checks the correct work of the mesh tying condition
*/
KRATOS_TEST_CASE_IN_SUITE(MeshTyingCondition2, KratosContactStructuralMechanicsFastSuite)
{
    // Create model part
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 3);

    // Fill test model part
    // GenerateReferenceModelPart(r_model_part);
    GenerateMeshTyingModelPart(r_model_part);

    // Solve system
    SolveSystem(r_model_part);

    // // Check results
    // IOMeshTyingDebug(r_model_part);

    // Check solutions
    CheckSolution(r_model_part);
}

}  // namespace Kratos::Testing.
