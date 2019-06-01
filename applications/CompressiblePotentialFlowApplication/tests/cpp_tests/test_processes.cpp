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
#include "compressible_potential_flow_application_variables.h"
#include "custom_processes/move_model_part_process.h"
#include "custom_processes/compute_embedded_lift_process.h"

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

    KRATOS_TEST_CASE_IN_SUITE(ComputeEmbeddedLiftProcess, CompressiblePotentialApplicationFastSuite)
    {

      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      // Set the element properties
      model_part.CreateNewProperties(0);
      Properties::Pointer pElemProp = model_part.pGetProperties(0);

      // Geometry creation
      model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
      model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
      model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
      std::vector<std::size_t> elemNodes{ 1, 2, 3 };
      model_part.CreateNewElement("EmbeddedIncompressiblePotentialFlowElement2D3N", 1, elemNodes, pElemProp);

      Element::Pointer pElement = model_part.pGetElement(1);
      pElement -> Set(TO_SPLIT);

      KRATOS_WATCH("HEY1")
      // Define the nodal values
      Vector potential(3);
      potential(0) = 1.0;
      potential(1) = 2.0;
      potential(2) = 3.0;

      Vector distances(3);
      potential(0) = -1.0;
      potential(1) = -1.0;
      potential(2) = 1.0;
      KRATOS_WATCH("HEY2")


      for (unsigned int i = 0; i < 3; i++)
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential(i);

      for (unsigned int i = 0; i < 3; i++)
        pElement->GetGeometry()[i].FastGetSolutionStepValue(GEOMETRY_DISTANCE) = distances(i);
      KRATOS_WATCH("HEY3")

      Vector resultant_force(3);
      ComputeEmbeddedLiftProcess ComputeEmbeddedLiftProcess(model_part, resultant_force);
      ComputeEmbeddedLiftProcess.Execute();
      KRATOS_WATCH("HEY4")

      // std::vector<double> reference({0.5, 0.0, -0.5});

      // for (unsigned int i = 0; i < 3; i++) {
      //   KRATOS_CHECK_NEAR(resultant_force(i), reference[i], 1e-6);
      // }
    }
  } // namespace Testing
}  // namespace Kratos.
