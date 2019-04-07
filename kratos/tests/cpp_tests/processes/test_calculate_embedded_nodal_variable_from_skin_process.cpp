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
//

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/checks.h"
#include "linear_solvers/linear_solver.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "processes/calculate_discontinuous_distance_to_skin_process.h"
#include "processes/calculate_embedded_nodal_variable_from_skin_process.h"
#include "spaces/ublas_space.h"

namespace Kratos {
namespace Testing {

    KRATOS_TEST_CASE_IN_SUITE(CalculateEmbeddedNodalVariableFromSkinProcessHorizontalPlane2D, KratosCoreFastSuite)
    {
        Model current_model;

        // Generate the element
        ModelPart &fluid_part = current_model.CreateModelPart("Surface");
        fluid_part.AddNodalSolutionStepVariable(DISTANCE);
        fluid_part.AddNodalSolutionStepVariable(TEMPERATURE);
        fluid_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        fluid_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        fluid_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        fluid_part.CreateNewNode(4, 1.0, 1.0, 0.0);
        Properties::Pointer p_properties_0(new Properties(0));
        fluid_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties_0);
        fluid_part.CreateNewElement("Element2D3N", 2, {2, 4, 3}, p_properties_0);

        // Generate the skin
        const double plane_height = 0.5;
        ModelPart &skin_part = current_model.CreateModelPart("Skin");
        skin_part.AddNodalSolutionStepVariable(TEMPERATURE);
        skin_part.CreateNewNode(1, -1.0, plane_height, 0.0);
        skin_part.CreateNewNode(2,  0.6, plane_height, 0.0);
        Properties::Pointer p_properties_1(new Properties(1));
        skin_part.CreateNewElement("Element2D2N", 1, {{1, 2}}, p_properties_1);
        skin_part.GetNode(1).FastGetSolutionStepValue(TEMPERATURE) = 1.0;
        skin_part.GetNode(1).FastGetSolutionStepValue(TEMPERATURE) = 2.0;

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(fluid_part, skin_part);
        disc_dist_proc.Execute();

        // Check values
        const auto &r_elem_dist = (fluid_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
        KRATOS_CHECK_NEAR(r_elem_dist[0], -0.5, 1e-12);
        KRATOS_CHECK_NEAR(r_elem_dist[1], -0.5, 1e-12);
        KRATOS_CHECK_NEAR(r_elem_dist[2],  0.5, 1e-12);

        // Compute the embedded nodal variable values
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
        typedef CalculateEmbeddedNodalVariableFromSkinProcess<2, double, SparseSpaceType, LocalSpaceType, LinearSolverType> EmbeddedNodalVariableProcess;
        EmbeddedNodalVariableProcess emb_nod_var_from_skin_proc(
            fluid_part,
            skin_part,
            TEMPERATURE,
            TEMPERATURE,
            "discontinuous");

        emb_nod_var_from_skin_proc.Execute();

        KRATOS_WATCH(fluid_part.GetNode(1).FastGetSolutionStepValue(DISTANCE))
        KRATOS_WATCH(fluid_part.GetNode(2).FastGetSolutionStepValue(DISTANCE))
        KRATOS_WATCH(fluid_part.GetNode(3).FastGetSolutionStepValue(DISTANCE))
        KRATOS_WATCH(fluid_part.GetNode(4).FastGetSolutionStepValue(DISTANCE))
        KRATOS_WATCH(fluid_part.GetNode(1).FastGetSolutionStepValue(TEMPERATURE))
        KRATOS_WATCH(fluid_part.GetNode(2).FastGetSolutionStepValue(TEMPERATURE))
        KRATOS_WATCH(fluid_part.GetNode(3).FastGetSolutionStepValue(TEMPERATURE))
        KRATOS_WATCH(fluid_part.GetNode(4).FastGetSolutionStepValue(TEMPERATURE))
    }

}  // namespace Testing.
}  // namespace Kratos.
