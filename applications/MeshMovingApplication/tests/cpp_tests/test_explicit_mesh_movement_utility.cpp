//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/checks.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "utilities/variable_utils.h"
#include "includes/mesh_moving_variables.h"

// Application includes
#include "custom_utilities/explicit_mesh_moving_utilities.h"

namespace Kratos {
namespace Testing {

    KRATOS_TEST_CASE_IN_SUITE(ExplicitMeshMovingUtilities2D, MeshMovingApplicationFastSuite)
    {
        Model current_model;

        // Generate the origin model part (done with the StructuredMeshGeneratorProcess)
        Node<3>::Pointer p_point_1 = Kratos::make_intrusive<Node<3>>(1, 0.0, 0.0, 0.0);
        Node<3>::Pointer p_point_2 = Kratos::make_intrusive<Node<3>>(2, 0.0, 1.0, 0.0);
        Node<3>::Pointer p_point_3 = Kratos::make_intrusive<Node<3>>(3, 1.0, 1.0, 0.0);
        Node<3>::Pointer p_point_4 = Kratos::make_intrusive<Node<3>>(4, 1.0, 0.0, 0.0);

        Quadrilateral2D4<Node<3> > geometry(p_point_1, p_point_2, p_point_3, p_point_4);

        Parameters mesher_parameters(R"(
        {
            "number_of_divisions": 3,
            "element_name": "Element2D3N"
        })");

        ModelPart& origin_model_part = current_model.CreateModelPart("OriginModelPart");
        origin_model_part.SetBufferSize(3);
        origin_model_part.AddNodalSolutionStepVariable(VELOCITY);
        origin_model_part.AddNodalSolutionStepVariable(PRESSURE);
        origin_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
        origin_model_part.AddNodalSolutionStepVariable(MESH_DISPLACEMENT);
        StructuredMeshGeneratorProcess(geometry, origin_model_part, mesher_parameters).Execute();

        // Fix the boundary mesh displacement
        const double max_x = 1.0;
        const double min_x = 0.0;
        const double max_y = 1.0;
        const double min_y = 0.0;
        const double coord_tol = 1e-3;
        for (auto it_node : origin_model_part.NodesArray()){
            if ((std::abs(it_node->X() - min_x) < coord_tol) || (std::abs(it_node->X() - max_x) < coord_tol))
                it_node->Fix(MESH_DISPLACEMENT_X);
            if ((std::abs(it_node->Y() - min_y) < coord_tol) || (std::abs(it_node->Y() - max_y) < coord_tol))
                it_node->Fix(MESH_DISPLACEMENT_Y);
        }

        // Create a fake time loop to fill the buffer
        const double n_steps = 3;
        const double delta_time = 0.1;
        for (unsigned int i_step = 0; i_step < n_steps; ++i_step){
            origin_model_part.CloneTimeStep(i_step * delta_time);
            double p_val = i_step * delta_time;
            array_1d<double,3> v_val = ZeroVector(3);
            for (auto it_node : origin_model_part.NodesArray()){
                it_node->GetSolutionStepValue(PRESSURE) = p_val;
                v_val(0) = i_step * delta_time * it_node->X();
                v_val(1) = i_step * delta_time * it_node->Y();
                it_node->GetSolutionStepValue(VELOCITY) = v_val;
            }
        }

        // Set the virtual model part
        ModelPart& virtual_model_part = current_model.CreateModelPart("VirtualModelPart");
        virtual_model_part.SetBufferSize(3);
        virtual_model_part.AddNodalSolutionStepVariable(VELOCITY);
        virtual_model_part.AddNodalSolutionStepVariable(PRESSURE);
        virtual_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
        virtual_model_part.AddNodalSolutionStepVariable(MESH_DISPLACEMENT);

        // Set the structure model part
        ModelPart& str_model_part =current_model.CreateModelPart("StructureModelPart");
        str_model_part.SetBufferSize(3);
        str_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        Properties::Pointer p_prop = Kratos::make_shared<Properties>(0);
        auto p_str_point_1 = str_model_part.CreateNewNode(100, 0.6, 0.2, 0.0);
        auto p_str_point_2 = str_model_part.CreateNewNode(200, 0.6, 0.4, 0.0);
        auto p_str_point_3 = str_model_part.CreateNewNode(300, 0.8, 0.4, 0.0);
        auto p_str_point_4 = str_model_part.CreateNewNode(400, 0.8, 0.2, 0.0);
        auto p_elem = str_model_part.CreateNewElement("Element2D4N", 1, {100,200,300,400}, p_prop);

        // Set the structure mesh movement
        array_1d<double,3> str_mov = ZeroVector(3);
        for (auto it_str_node : str_model_part.NodesArray()){
            str_mov(0) = -(it_str_node->X()) * 0.1;
            str_mov(1) = (it_str_node->Y()) * 0.2;
            it_str_node->FastGetSolutionStepValue(DISPLACEMENT) = str_mov;
        }

        // Set the explicit mesh moving utility
        const double search_radius = 1.0;
        ExplicitMeshMovingUtilities::Pointer p_mesh_moving =
            Kratos::make_shared<ExplicitMeshMovingUtilities>(virtual_model_part, str_model_part, search_radius);

        // Fill the virtual model part geometry
        p_mesh_moving->FillVirtualModelPart(origin_model_part);

        // Copy the origin model part data to the virtual one
        VariableUtils::Pointer p_var_utils = Kratos::make_shared<VariableUtils>();
        for (unsigned int i_step = 0; i_step < n_steps; ++i_step){
            p_var_utils->CopyModelPartNodalVar<Variable<double>>(PRESSURE, origin_model_part, virtual_model_part, i_step);
            p_var_utils->CopyModelPartNodalVar<Variable<array_1d<double,3>>>(VELOCITY, origin_model_part, virtual_model_part, i_step);
        }

        // Execute the explicit mesh movement operations
        const unsigned int buffer_size = 3;
        p_mesh_moving->ComputeExplicitMeshMovement(delta_time);
        p_mesh_moving->ProjectVirtualValues<2>(origin_model_part, buffer_size);
        p_mesh_moving->UndoMeshMovement();

        // Check the obtained results
        auto p_node_10 = origin_model_part.pGetNode(10);
        const auto &r_vel_n1 = p_node_10->GetSolutionStepValue(VELOCITY,1);
        const auto &r_mesh_vel = p_node_10->FastGetSolutionStepValue(MESH_VELOCITY);

        const double tol = 1e-6;
        KRATOS_CHECK_NEAR(r_vel_n1(0), 0.0710645, tol);
        KRATOS_CHECK_NEAR(r_vel_n1(1), 0.0294852, tol);
        KRATOS_CHECK_NEAR(r_mesh_vel(0), -0.439785, tol);
        KRATOS_CHECK_NEAR(r_mesh_vel(1), 0.384812, tol);

    }
}
}  // namespace Kratos.
