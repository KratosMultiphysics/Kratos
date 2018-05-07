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
#include "testing/testing.h"
#include "includes/checks.h"
#include "includes/gid_io.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "utilities/variable_utils.h"

// Application includes
#include "mesh_moving_application.h"
#include "custom_utilities/explicit_mesh_moving_utilities.h"

namespace Kratos {
namespace Testing {

    KRATOS_TEST_CASE_IN_SUITE(ExplicitMeshMovingUtilities2D, MeshMovingApplicationFastSuite)
    {

        // Generate the origin model part (done with the StructuredMeshGeneratorProcess)
        Node<3>::Pointer p_point_1 = Kratos::make_shared<Node<3>>(1, 0.0, 0.0, 0.0);
        Node<3>::Pointer p_point_2 = Kratos::make_shared<Node<3>>(2, 0.0, 1.0, 0.0);
        Node<3>::Pointer p_point_3 = Kratos::make_shared<Node<3>>(3, 1.0, 1.0, 0.0);
        Node<3>::Pointer p_point_4 = Kratos::make_shared<Node<3>>(4, 1.0, 0.0, 0.0);

        Quadrilateral2D4<Node<3> > geometry(p_point_1, p_point_2, p_point_3, p_point_4);

        Parameters mesher_parameters(R"(
        {
            "number_of_divisions": 2,
            "element_name": "Element2D3N"
        })");

        ModelPart origin_model_part("OriginModelPart");
        origin_model_part.SetBufferSize(3);
        origin_model_part.AddNodalSolutionStepVariable(VELOCITY);
        origin_model_part.AddNodalSolutionStepVariable(PRESSURE);
        origin_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
        origin_model_part.AddNodalSolutionStepVariable(MESH_DISPLACEMENT);
        StructuredMeshGeneratorProcess(geometry, origin_model_part, mesher_parameters).Execute();

        // Create a fake time loop to fill the buffer
        const double n_steps = 3;
        const double delta_time = 0.1;
        for (unsigned int i_step = 0; i_step < n_steps; ++i_step){
            origin_model_part.CloneTimeStep(i_step * delta_time);
        }

        // Set the virtual model part
        ModelPart virtual_model_part("VirtualModelPart");
        virtual_model_part.SetBufferSize(3);
        virtual_model_part.AddNodalSolutionStepVariable(VELOCITY);
        virtual_model_part.AddNodalSolutionStepVariable(PRESSURE);
        virtual_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
        virtual_model_part.AddNodalSolutionStepVariable(MESH_DISPLACEMENT);

        // Set the structure model part
        ModelPart str_model_part("StructureModelPart");
        str_model_part.SetBufferSize(3);
        str_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        Properties::Pointer p_prop = Kratos::make_shared<Properties>(0);
        auto p_str_point_1 = str_model_part.CreateNewNode(1, 0.6, 0.2, 0.0);
        auto p_str_point_2 = str_model_part.CreateNewNode(2, 0.6, 0.4, 0.0);
        auto p_str_point_3 = str_model_part.CreateNewNode(3, 0.8, 0.4, 0.0);
        auto p_str_point_4 = str_model_part.CreateNewNode(4, 0.8, 0.2, 0.0);
        auto p_elem = str_model_part.CreateNewElement("Element2D4N", 1, {1,2,3,4}, p_prop);

        // Set the structure mesh movement
        array_1d<double,3> str_mov = ZeroVector(3);
        str_mov(0) = -0.4;
        str_mov(1) = 0.4;
        for (auto it_str_node : str_model_part.NodesArray()){
            it_str_node->FastGetSolutionStepValue(DISPLACEMENT);
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
    }
}
}  // namespace Kratos.
