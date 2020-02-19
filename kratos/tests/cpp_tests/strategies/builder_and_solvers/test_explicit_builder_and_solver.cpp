//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes
#include <limits>
#include <iomanip>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/define.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "geometries/line_2d_2.h"
#include "tests/cpp_tests/auxiliar_files_for_cpp_unnitest/test_bar_element.h"
#include "solving_strategies/builder_and_solvers/explicit_builder_and_solver.h"

namespace Kratos
{
namespace Testing
{

/// Tests
typedef Node<3> NodeType;
typedef Geometry<NodeType> GeometryType;
typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;

// The builder ans solver type
typedef ExplicitBuilderAndSolver< SparseSpaceType, LocalSpaceType > ExplicitBuilderAndSolverType;

/**
 * @brief It generates a truss structure with an expected solution
 */
static inline void GenerateTestModelPart(
    ModelPart& rModelPart,
    const bool FixDofs = false)
{
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(REACTION);

    // Create the test elements
    auto p_prop = rModelPart.CreateNewProperties(1, 0);
    auto p_node_1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_node_3 = rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    auto p_geom_1 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({p_node_1, p_node_2})});
    auto p_geom_2 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({p_node_2, p_node_3})});
    rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>(1, p_geom_1, p_prop));
    rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>(2, p_geom_2, p_prop));

    /// Add dof
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X, REACTION_X);
        r_node.AddDof(DISPLACEMENT_Y, REACTION_Y);
        r_node.AddDof(DISPLACEMENT_Z, REACTION_Z);
    }

    /// Initialize elements
    const auto& r_process_info = rModelPart.GetProcessInfo();
    for (auto& r_elem : rModelPart.Elements()) {
        r_elem.Initialize(r_process_info);
        r_elem.InitializeSolutionStep(r_process_info);
    }

    // Set initial solution
    for (auto& r_node : rModelPart.Nodes()) {
        (r_node.FastGetSolutionStepValue(DISPLACEMENT)) = ZeroVector(3);
        (r_node.FastGetSolutionStepValue(DISPLACEMENT, 1)) = ZeroVector(3);
        (r_node.FastGetSolutionStepValue(DISPLACEMENT, 2)) = ZeroVector(3);
    }

    // Fix dofs
    if (FixDofs) {
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.Fix(DISPLACEMENT_Y);
            r_node.Fix(DISPLACEMENT_Z);
        }
    }
}

/**
 * Checks if the explicit builder and solver performs correctly the assemble of the DOF set and the lumped mass vector
 */
KRATOS_TEST_CASE_IN_SUITE(ExplicitBlockBuilderAndSolverDofSet, KratosCoreFastSuite)
{
    // Generate the test model part
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("TestModelPart", 3);
    GenerateTestModelPart(r_model_part);

    // Test the explicit builder and solver
    auto p_builder_and_solver = Kratos::make_unique<ExplicitBuilderAndSolverType>();
    p_builder_and_solver->SetUpDofSet(r_model_part);
    p_builder_and_solver->SetUpDofSetEquationIds();
    const auto& r_dof_set = p_builder_and_solver->GetDofSet();
    // KRATOS_WATCH(r_dof_set)

    p_builder_and_solver->SetUpLumpedMassMatrixVector(r_model_part);
    const auto& r_lumped_mass_vector = p_builder_and_solver->GetLumpedMassMatrixVector();
    KRATOS_WATCH(r_lumped_mass_vector)

    // // The solution check
    // constexpr double tolerance = 1e-8;
    // KRATOS_CHECK(rA.size1() == 6);
    // KRATOS_CHECK(rA.size2() == 6);
    // KRATOS_CHECK_LESS_EQUAL(std::abs((rA(0,0) - 2069000000.000000000)/rA(0,0)), tolerance);
    // KRATOS_CHECK_LESS_EQUAL(std::abs((rA(1,1) - 1.000000000)/rA(1,1)), tolerance);
    // KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,2) - 4138000000.000000000)/rA(2,2)), tolerance);
    // KRATOS_CHECK_LESS_EQUAL(std::abs((rA(2,4) - -2069000000.000000000)/rA(2,4)), tolerance);
    // KRATOS_CHECK_LESS_EQUAL(std::abs((rA(3,3) - 1.000000000)/rA(3,3)), tolerance);
    // KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,2) - -2069000000.000000000)/rA(4,2)), tolerance);
    // KRATOS_CHECK_LESS_EQUAL(std::abs((rA(4,4) - 2069000000.000000000)/rA(4,4)), tolerance);
    // KRATOS_CHECK_LESS_EQUAL(std::abs((rA(5,5) - 1.000000000)/rA(5,5)), tolerance);

}

} // namespace Testing
} // namespace Kratos.
