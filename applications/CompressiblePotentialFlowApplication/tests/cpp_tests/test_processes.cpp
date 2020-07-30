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
#include "custom_processes/define_2d_wake_process.h"
#include "custom_processes/apply_far_field_process.h"
#include "custom_processes/compute_nodal_value_process.h"

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

      model_part.AddNodalSolutionStepVariable(GEOMETRY_DISTANCE);
      model_part.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
      model_part.AddNodalSolutionStepVariable(AUXILIARY_VELOCITY_POTENTIAL);

      array_1d<double, 3> free_stream_velocity;
      free_stream_velocity[0] = 1.0; free_stream_velocity[1] = 0.0; free_stream_velocity[2] = 0.0;
      model_part.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;

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
      pElement -> Set(ACTIVE);

      // Define the nodal values
      std::array<double,3> potential;
      potential[0] = 1.0;
      potential[1] = 2.0;
      potential[2] = 3.0;

      std::array<double,3> distances;
      distances[0] = -1.0;
      distances[1] = -1.0;
      distances[2] = 1.0;


      for (unsigned int i = 0; i < 3; i++)
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];

      for (unsigned int i = 0; i < 3; i++)
        pElement->GetGeometry()[i].FastGetSolutionStepValue(GEOMETRY_DISTANCE) = distances[i];

      Vector resultant_force(3);
      ComputeEmbeddedLiftProcess<2,3>(model_part, resultant_force).Execute();

      std::array<double,3> reference{0.0, 0.5, 0.0};
      KRATOS_WATCH(resultant_force)

      for (unsigned int i = 0; i < 3; i++) {
        KRATOS_CHECK_NEAR(resultant_force(i), reference[i], 1e-6);
      }
    }

    KRATOS_TEST_CASE_IN_SUITE(Define2DWakeProcessProcess, CompressiblePotentialApplicationFastSuite)
    {
      // Create model_part
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

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

    KRATOS_TEST_CASE_IN_SUITE(ApplyFarFieldProcess, CompressiblePotentialApplicationFastSuite)
    {
      // Create model_part
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      // Variables addition
      model_part.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
      model_part.AddNodalSolutionStepVariable(AUXILIARY_VELOCITY_POTENTIAL);

      // Set model_part properties
      BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
      free_stream_velocity(0) = 10.0;
      model_part.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;

      // Create nodes
      model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
      model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
      model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
      model_part.CreateNewNode(4, 1.0, 1.0, 0.0);

      for (auto& r_node : model_part.Nodes()) {
        r_node.AddDof(VELOCITY_POTENTIAL);
      }

      model_part.CreateNewProperties(0);
      Properties::Pointer pProp = model_part.pGetProperties(0);

      std::vector<ModelPart::IndexType> cond1{1, 2};
      std::vector<ModelPart::IndexType> cond2{2, 4};
      std::vector<ModelPart::IndexType> cond3{4, 3};
      std::vector<ModelPart::IndexType> cond4{3, 1};

      model_part.CreateNewCondition("LineCondition2D2N", 1, cond1, pProp);
      model_part.CreateNewCondition("LineCondition2D2N", 2, cond2, pProp);
      model_part.CreateNewCondition("LineCondition2D2N", 3, cond3, pProp);
      model_part.CreateNewCondition("LineCondition2D2N", 4, cond4, pProp);


      // Set initial potential
      const double initial_potential = 1.0;
      const bool initialize_flow_field = true;
      const bool perturbation_field = false;

      // Construct the ApplyFarFieldProcess
      ApplyFarFieldProcess ApplyFarFieldProcess(model_part, initial_potential, initialize_flow_field, perturbation_field);

      // Execute the ApplyFarFieldProcess
      ApplyFarFieldProcess.Execute();

      for (auto& r_node : model_part.Nodes()) {
        if (r_node.Id() == 1 || r_node.Id() == 3) {
          KRATOS_CHECK(r_node.IsFixed(VELOCITY_POTENTIAL));
          KRATOS_CHECK_NEAR(r_node.FastGetSolutionStepValue(VELOCITY_POTENTIAL), initial_potential, 1e-6);
        }
        else {
          KRATOS_CHECK_IS_FALSE(r_node.IsFixed(VELOCITY_POTENTIAL));
          KRATOS_CHECK_NEAR(r_node.FastGetSolutionStepValue(VELOCITY_POTENTIAL), 11.0, 1e-6);
        }
      }
    }

    KRATOS_TEST_CASE_IN_SUITE(ComputeNodalValueProcess, CompressiblePotentialApplicationFastSuite)
    {
      // Create model_part
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);
      model_part.GetProcessInfo()[DOMAIN_SIZE] = 2;

      // Set model_part properties
      BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
      free_stream_velocity(0) = 10.0;
      model_part.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;

      // Variables addition
      model_part.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
      model_part.AddNodalSolutionStepVariable(AUXILIARY_VELOCITY_POTENTIAL);

      // Set the element properties
      model_part.CreateNewProperties(0);
      Properties::Pointer pElemProp = model_part.pGetProperties(0);

      // Geometry creation
      model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
      model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
      model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
      std::vector<ModelPart::IndexType> elemNodes_1{ 1, 2, 3 };
      model_part.CreateNewElement("IncompressiblePotentialFlowElement2D3N", 1, elemNodes_1, pElemProp);

      auto& r_element = model_part.GetElement(1);

      r_element.GetValue(WAKE) = 1;
      Vector wake_elemental_distances(3);
      wake_elemental_distances(0) = 1.0;
      wake_elemental_distances(1) = 1.0;
      wake_elemental_distances(2) = -1.0;
      r_element.GetValue(WAKE_ELEMENTAL_DISTANCES) = wake_elemental_distances;

      // Define the nodal values
      Vector potential(3);
      potential(0) = 1.0;
      potential(1) = 2.0;
      potential(2) = 3.0;

      for (unsigned int i = 0; i < 3; i++){
        if (wake_elemental_distances(i) > 0.0)
          r_element.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential(i);
        else
          r_element.GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential(i);
      }
      for (unsigned int i = 0; i < 3; i++){
        if (wake_elemental_distances(i) < 0.0)
          r_element.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential(i)+1;
        else
          r_element.GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential(i)+1;
      }
      // Construct the ComputeNodalPotentialFlowVelocityProcess
      const std::vector<std::string> variable_array = {"VELOCITY","PRESSURE_COEFFICIENT"};
      ComputeNodalValueProcess ComputeNodalValueProcess(model_part, variable_array);

      // Execute the ComputeNodalPotentialFlowVelocityProcess
      ComputeNodalValueProcess.Execute();
      for (auto& r_node : model_part.Nodes()) {
        auto nodal_area = r_node.GetValue(NODAL_AREA);
        KRATOS_CHECK_NEAR(nodal_area, 0.166667, 1e-6);
        auto nodal_velocity = r_node.GetValue(VELOCITY);
        KRATOS_CHECK_NEAR(nodal_velocity[0], 1, 1e-6);
        KRATOS_CHECK_NEAR(nodal_velocity[1], 1, 1e-6);
        KRATOS_CHECK_NEAR(nodal_velocity[2], 0, 1e-6);
        auto nodal_pressure = r_node.GetValue(PRESSURE_COEFFICIENT);
        KRATOS_CHECK_NEAR(nodal_pressure, 0.98, 1e-6);
      }
    }
  } // namespace Testing
}  // namespace Kratos.
