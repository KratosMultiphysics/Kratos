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
#include "solving_strategies/builder_and_solvers/explicit_builder.h"

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
typedef ExplicitBuilder< SparseSpaceType, LocalSpaceType > ExplicitBuilderType;

/**
 * @brief It generates a truss structure with an expected solution
 */
static inline void GenerateTestExplicitBuilderModelPart(
    ModelPart& rModelPart,
    const bool FixDofs = false)
{
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(REACTION);

    // Create the test elements
    auto p_prop = rModelPart.CreateNewProperties(1, 0);
    p_prop->SetValue(DENSITY, 1.0);
    p_prop->SetValue(NODAL_AREA, 1.0);
    auto p_node_1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_node_3 = rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    auto p_geom_1 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({p_node_1, p_node_2})});
    auto p_geom_2 = Kratos::make_shared<Line2D2<NodeType>>(PointerVector<NodeType>{std::vector<NodeType::Pointer>({p_node_2, p_node_3})});
    rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>(1, p_geom_1, p_prop));
    rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>(2, p_geom_2, p_prop));

    /// Add dof
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X, REACTION_X);
        r_node.AddDof(DISPLACEMENT_Y, REACTION_Y);
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
        }
    }
}

/**
 * Checks if the explicit builder and solver performs correctly the assemble of the DOF set
 */
KRATOS_TEST_CASE_IN_SUITE(ExplicitBlockBuilderAndSolverInitialization, KratosCoreFastSuite)
{
    // Generate the test model part
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("TestModelPart", 3);
    GenerateTestExplicitBuilderModelPart(r_model_part);

    // Test the explicit builder and solver dof set set up
    auto p_builder_and_solver = Kratos::make_unique<ExplicitBuilderType>();
    p_builder_and_solver->Initialize(r_model_part);
    const auto& r_dof_set = p_builder_and_solver->GetDofSet();
    const auto& r_lumped_mass_vector = p_builder_and_solver->GetLumpedMassMatrixVector();

    // Check the DOF set
    KRATOS_CHECK_EQUAL(r_dof_set.size(), 6);
    for (unsigned int i_dof = 0; i_dof < r_dof_set.size(); ++i_dof) {
        const auto it_dof = r_dof_set.begin() + i_dof;
        if (i_dof % 2 == 0) {
            KRATOS_CHECK(it_dof->GetVariable() == DISPLACEMENT_X);
            KRATOS_CHECK(it_dof->GetReaction() == REACTION_X);
        } else {
            KRATOS_CHECK(it_dof->GetVariable() == DISPLACEMENT_Y);
            KRATOS_CHECK(it_dof->GetReaction() == REACTION_Y);
        }
    }

    // Check the lumped mass vector values
    const double tolerance = 1.0e-8;
    std::vector<double> expected_solution = {0.5, 0.0, 0.5, 0.5, 0.0, 0.5};
    KRATOS_CHECK_VECTOR_NEAR(r_lumped_mass_vector, expected_solution, tolerance);
}

} // namespace Testing
} // namespace Kratos.
