//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Núñez
//
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_processes/move_model_part_process.h"
#include "custom_processes/define_2d_wake_process.h"

namespace Kratos {
  namespace Testing {

    KRATOS_TEST_CASE_IN_SUITE(MoveModelModelPartProcess, CompressiblePotentialApplicationFastSuite)
    {

      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      // Nodes creation
      model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
      model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
      model_part.CreateNewNode(3, -1.0, 0.0, 0.0);

      // Process parameters
      Parameters moving_parameters = Parameters(R"(
        {
            "origin"                        : [5.0,5.0,0.0],
            "sizing_multiplier"             : 2.0

        })" );

      moving_parameters.AddEmptyValue("rotation_angle");
      moving_parameters["rotation_angle"].SetDouble(Globals::Pi/2);

      MoveModelPartProcess MoveModelPartProcess(model_part, moving_parameters);
      MoveModelPartProcess.Execute();

      Matrix reference = ZeroMatrix(3,3);
      reference(0,0) = 5.0; reference(0,1) = 5.0;
      reference(1,0) = 5.0; reference(1,1) = 7.0;
      reference(2,0) = 5.0; reference(2,1) = 3.0;

      for (std::size_t i_node = 0; i_node<reference.size1(); i_node++){
        for (std::size_t i_dim = 0; i_dim<reference.size2(); i_dim++){
          KRATOS_CHECK_NEAR(model_part.GetNode(i_node+1).Coordinates()[i_dim], reference(i_node,i_dim), 1e-6);
        }
      }
    }

    KRATOS_TEST_CASE_IN_SUITE(Define2DWakeProcessProcess, CompressiblePotentialApplicationFastSuite)
    {
      // Create model_part
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);
      model_part.AddNodalSolutionStepVariable(DISTANCE);

      // Set model_part properties
      BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
      free_stream_velocity(0) = 10.0;
      model_part.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;

      // Create nodes
      model_part.CreateNewNode(1, 2.0, 0.0, 0.0);
      model_part.CreateNewNode(2, 2.0, 2.0, 0.0);
      ModelPart::NodeType::Pointer pNode = model_part.CreateNewNode(3, 0.0, 1.0, 0.0);

      // Create element
      model_part.CreateNewProperties(0);
      Properties::Pointer pElemProp = model_part.pGetProperties(0);
      std::vector<ModelPart::IndexType> elemNodes{ 1, 2, 3 };
      ModelPart::ElementType::Pointer pElement = model_part.CreateNewElement("IncompressiblePotentialFlowElement2D3N", 1, elemNodes, pElemProp);

      // Create body sub_model_part
      ModelPart& body_model_part = model_part.CreateSubModelPart("body_model_part");
      body_model_part.AddNode(pNode);

      // Set Tolerance
      const double tolerance = 1e-9;

      // Construct the Define2DWakeProcess
      Define2DWakeProcess Define2DWakeProcess(body_model_part, tolerance);

      // Execute the Define2DWakeProcess
      Define2DWakeProcess.ExecuteInitialize();

      const int wake = pElement->GetValue(WAKE);
      KRATOS_CHECK_NEAR(wake, 1, 1e-6);
    }
  } // namespace Testing
}  // namespace Kratos.
