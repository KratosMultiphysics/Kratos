//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez
//
//

// Project includes
#include "tests/cpp_tests/compressible_potential_flow_fast_suite.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_elements/incompressible_potential_flow_element.h"
#include "custom_elements/embedded_incompressible_potential_flow_element.h"
#include "custom_utilities/potential_flow_utilities.h"
#include "tests/cpp_tests/test_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos {
  namespace Testing {

    typedef ModelPart::IndexType IndexType;
    typedef ModelPart::NodeIterator NodeIteratorType;

    void GenerateElement(ModelPart& rModelPart)
    {
      // Variables addition
      rModelPart.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
      rModelPart.AddNodalSolutionStepVariable(AUXILIARY_VELOCITY_POTENTIAL);

      // Set the element properties
      rModelPart.CreateNewProperties(0);
      Properties::Pointer pElemProp = rModelPart.pGetProperties(0);
      BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
      free_stream_velocity(0) = 10.0;

      rModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;
      rModelPart.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.0;

      // Geometry creation
      rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
      rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
      rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
      std::vector<ModelPart::IndexType> elemNodes{ 1, 2, 3 };
      rModelPart.CreateNewElement("IncompressiblePotentialFlowElement2D3N", 1, elemNodes, pElemProp);
    }

    void GenerateEmbeddedElement(ModelPart& rModelPart)
    {
      // Variables addition
      rModelPart.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
      rModelPart.AddNodalSolutionStepVariable(AUXILIARY_VELOCITY_POTENTIAL);
      rModelPart.AddNodalSolutionStepVariable(GEOMETRY_DISTANCE);


      // Set the element properties
      rModelPart.CreateNewProperties(0);
      Properties::Pointer pElemProp = rModelPart.pGetProperties(0);
      rModelPart.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.0;

      // Geometry creation
      rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
      rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
      rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
      std::vector<ModelPart::IndexType> elemNodes{ 1, 2, 3 };
      rModelPart.CreateNewElement("EmbeddedIncompressiblePotentialFlowElement2D3N", 1, elemNodes, pElemProp);
    }

    void AssignPotentialsToNormalElement(Element::Pointer pElement)
    {
      // Define the nodal values
      Vector potential(3);
      potential(0) = 1.0;
      potential(1) = 2.0;
      potential(2) = 3.0;

      for (unsigned int i = 0; i < 3; i++)
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential(i);
    }

    void AssignPotentialsToWakeElement(Element::Pointer pElement, const array_1d<double, 3>& rDistances)
    {
      // Define the nodal values
      Vector potential(3);
      potential(0) = 1.0;
      potential(1) = 2.0;
      potential(2) = 3.0;

      for (unsigned int i = 0; i < 3; i++){
        if (rDistances(i) > 0.0)
          pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential(i);
        else
          pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential(i);
      }
      for (unsigned int i = 0; i < 3; i++){
        if (rDistances(i) < 0.0)
          pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential(i)+5;
        else
          pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential(i)+5;
      }
    }

    BoundedVector<double,3> AssignDistances()
    {
      BoundedVector<double,3> distances;
      distances(0) = 1.0;
      distances(1) = -1.0;
      distances(2) = -1.0;
      return distances;
    }

    /** Checks the IncompressiblePotentialFlowElement element.
     * Checks the LHS and RHS computation.
     */
    KRATOS_TEST_CASE_IN_SUITE(IncompressiblePotentialFlowElementCalculateLocalSystem, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);

      AssignPotentialsToNormalElement(p_element);

      // Compute RHS and LHS
      Vector RHS = ZeroVector(3);
      Matrix LHS = ZeroMatrix(3, 3);

      const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
      p_element->CalculateLocalSystem(LHS, RHS, r_current_process_info);

      // Check the RHS values (the RHS is computed as the LHS x previous_solution,
      // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
      std::array<double, 3> reference{0.5, 0.0, -0.5};

      for (unsigned int i = 0; i < RHS.size(); i++) {
        KRATOS_EXPECT_NEAR(RHS(i), reference[i], 1e-6);
      }
    }

    KRATOS_TEST_CASE_IN_SUITE(IncompressiblePotentialFlowElementCalculateLocalSystemWake, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);

      BoundedVector<double,3> distances = AssignDistances();

      p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
      p_element->GetValue(WAKE) = true;

      AssignPotentialsToWakeElement(p_element, distances);

      // Compute RHS and LHS
      Vector RHS = ZeroVector(6);
      Matrix LHS = ZeroMatrix(6, 6);

      const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
      p_element->CalculateLocalSystem(LHS, RHS, r_current_process_info);

      // Check the RHS values (the RHS is computed as the LHS x previous_solution,
      // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
      std::array<double, 6> reference{0.5, 0.0, 0.0, 0.0, 0.0, -0.5};

      for (unsigned int i = 0; i < RHS.size(); i++) {
        KRATOS_EXPECT_NEAR(RHS(i), reference[i], 1e-6);
      }
    }

    KRATOS_TEST_CASE_IN_SUITE(EmbeddedIncompressiblePotentialFlowElementCalculateLocalSystem, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateEmbeddedElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);

      AssignPotentialsToNormalElement(p_element);

      // Define the distance values
      Vector level_set(3);
      level_set(0) = 1.0;
      level_set(1) = -1.0;
      level_set(2) = -1.0;

      for (unsigned int i = 0; i < 3; i++){
        p_element->GetGeometry()[i].FastGetSolutionStepValue(GEOMETRY_DISTANCE) = level_set(i);
      }

      // Compute RHS and LHS
      Vector RHS = ZeroVector(3);
      Matrix LHS = ZeroMatrix(3, 3);

      const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
      p_element->CalculateLocalSystem(LHS, RHS, r_current_process_info);

      // Check the RHS values (the RHS is computed as the LHS x previous_solution,
      // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
      std::array<double, 3> reference{0.125, 0.0, -0.125};

      for (unsigned int i = 0; i < RHS.size(); i++) {
        KRATOS_EXPECT_NEAR(RHS(i), reference[i], 1e-6);
      }
    }

    KRATOS_TEST_CASE_IN_SUITE(PingEmbeddedIncompressiblePotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateEmbeddedElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);
      const unsigned int number_of_nodes = p_element->GetGeometry().size();

      // Define the distance values
      std::array<double, 3> level_set{1.0, -1.0, -1.0};
      for (unsigned int i = 0; i < 3; i++) {
          p_element->GetGeometry()[i].FastGetSolutionStepValue(GEOMETRY_DISTANCE) = level_set[i];
      }

      const std::array<double, 3> potential{1.0, 20.0, 50.0};

      Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
      Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

      PotentialFlowTestUtilities::ComputeElementalSensitivities<3>(model_part, LHS_finite_diference, LHS_analytical, potential);

      KRATOS_EXPECT_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
    }

    KRATOS_TEST_CASE_IN_SUITE(PingEmbeddedIncompressiblePotentialFlowElementLHSPenalty, CompressiblePotentialApplicationFastSuite) {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateEmbeddedElement(model_part);
      model_part.GetProcessInfo()[PENALTY_COEFFICIENT] = 100.0;
      Element::Pointer p_element = model_part.pGetElement(1);
      const unsigned int number_of_nodes = p_element->GetGeometry().size();

      // Define the distance values
      std::array<double, 3> level_set{1.0, -1.0, -1.0};
      for (unsigned int i = 0; i < 3; i++) {
          p_element->GetGeometry()[i].FastGetSolutionStepValue(GEOMETRY_DISTANCE) = level_set[i];
      }

      const std::array<double, 3> potential{1.0, 20.0, 50.0};

      Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
      Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

      PotentialFlowTestUtilities::ComputeElementalSensitivities<3>(model_part, LHS_finite_diference, LHS_analytical, potential);

      KRATOS_EXPECT_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
    }


    /** Checks the IncompressiblePotentialFlowElement element.
 * Checks the EquationIdVector.
 */
    KRATOS_TEST_CASE_IN_SUITE(IncompressiblePotentialFlowElementEquationIdVector, CompressiblePotentialApplicationFastSuite)
    {

      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);

      for (unsigned int i = 0; i < 3; i++)
        p_element->GetGeometry()[i].AddDof(VELOCITY_POTENTIAL);

      Element::DofsVectorType ElementalDofList;
      const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
      p_element->GetDofList(ElementalDofList, r_current_process_info);

      for (int i = 0; i < 3; i++)
        ElementalDofList[i]->SetEquationId(i);

      Element::EquationIdVectorType EquationIdVector;
      p_element->EquationIdVector(EquationIdVector, r_current_process_info);

      // Check the EquationIdVector values
      for (unsigned int i = 0; i < EquationIdVector.size(); i++) {
        KRATOS_EXPECT_TRUE(EquationIdVector[i] == i);
      }
    }

    /** Checks the IncompressiblePotentialFlowElement element.
 * Checks the EquationIdVector for the Wake.
 */
    KRATOS_TEST_CASE_IN_SUITE(IncompressiblePotentialFlowElementEquationIdVectorWake, CompressiblePotentialApplicationFastSuite)
    {

      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);
      p_element->SetValue(WAKE, true);

      BoundedVector<double,3> distances = AssignDistances();
      p_element->SetValue(WAKE_ELEMENTAL_DISTANCES, distances);

      for (unsigned int i = 0; i < 3; i++) {
        p_element->GetGeometry()[i].AddDof(VELOCITY_POTENTIAL);
        p_element->GetGeometry()[i].AddDof(AUXILIARY_VELOCITY_POTENTIAL);
      }

      Element::DofsVectorType ElementalDofList;
      const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
      p_element->GetDofList(ElementalDofList, r_current_process_info);

      for (int i = 0; i < 6; i++)
        ElementalDofList[i]->SetEquationId(i);

      Element::EquationIdVectorType EquationIdVector;
      p_element->EquationIdVector(EquationIdVector, r_current_process_info);

      //Check the EquationIdVector values
      for (unsigned int i = 0; i < EquationIdVector.size(); i++) {
        KRATOS_EXPECT_TRUE(EquationIdVector[i] == i);
      }
    }

    // Checks the function GetWakeDistances
    KRATOS_TEST_CASE_IN_SUITE(GetWakeDistances, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);
      p_element->SetValue(WAKE, true);

      BoundedVector<double,3> distances = AssignDistances();
      p_element->SetValue(WAKE_ELEMENTAL_DISTANCES, distances);

      const auto returned_distances = PotentialFlowUtilities::GetWakeDistances<2, 3>(*p_element);

      std::array<double, 3> reference{1.0, -1.0, -1.0};

      for (unsigned int i = 0; i < returned_distances.size(); i++) {
        KRATOS_EXPECT_NEAR(returned_distances(i), reference[i], 1e-7);
      }
    }

    // Checks the function GetPotentialOnNormalElement
    KRATOS_TEST_CASE_IN_SUITE(GetPotentialOnNormalElement, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);

      AssignPotentialsToNormalElement(p_element);

      auto potentials = PotentialFlowUtilities::GetPotentialOnNormalElement<2,3>(*p_element);

      std::array<double, 3> reference{1.0, 2.0, 3.0};

      for (unsigned int i = 0; i < potentials.size(); i++) {
        KRATOS_EXPECT_NEAR(potentials(i), reference[i], 1e-7);
      }
    }

    // Checks the function GetPotentialOnUpperWakeElement
    KRATOS_TEST_CASE_IN_SUITE(GetPotentialOnUpperWakeElement, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);
      p_element->SetValue(WAKE, true);

      BoundedVector<double,3> distances = AssignDistances();

      AssignPotentialsToWakeElement(p_element, distances);

      array_1d<double, 3> potentials =
          PotentialFlowUtilities::GetPotentialOnUpperWakeElement<2, 3>(*p_element, distances);

      std::array<double, 3> reference{1.0, 2.0, 3.0};

      for (unsigned int i = 0; i < potentials.size(); i++) {
        KRATOS_EXPECT_NEAR(potentials(i), reference[i], 1e-7);
      }
    }

    // Checks the function GetPotentialOnLowerWakeElement
    KRATOS_TEST_CASE_IN_SUITE(GetPotentialOnLowerWakeElement, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);
      p_element->SetValue(WAKE, true);

      BoundedVector<double,3> distances = AssignDistances();

      AssignPotentialsToWakeElement(p_element, distances);

      array_1d<double, 3> potentials =
          PotentialFlowUtilities::GetPotentialOnLowerWakeElement<2, 3>(*p_element, distances);

      std::array<double, 3> reference{6.0, 7.0, 8.0};

      for (unsigned int i = 0; i < potentials.size(); i++) {
        KRATOS_EXPECT_NEAR(potentials(i), reference[i], 1e-7);
      }
    }

    // Checks the function GetPotentialOnWakeElement
    KRATOS_TEST_CASE_IN_SUITE(GetPotentialOnWakeElement, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);
      p_element->SetValue(WAKE, true);

      BoundedVector<double,3> distances = AssignDistances();

      AssignPotentialsToWakeElement(p_element, distances);

      BoundedVector<double, 2 * 3> potentials =
          PotentialFlowUtilities::GetPotentialOnWakeElement<2, 3>(*p_element, distances);

      std::array<double, 6> reference{1.0, 2.0, 3.0, 6.0, 7.0, 8.0};

      for (unsigned int i = 0; i < potentials.size(); i++) {
        KRATOS_EXPECT_NEAR(potentials(i), reference[i], 1e-7);
      }
    }

    // Checks the function ComputeVelocityNormalElement
    KRATOS_TEST_CASE_IN_SUITE(ComputeVelocityNormalElement, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);

      AssignPotentialsToNormalElement(p_element);

      auto velocity = PotentialFlowUtilities::ComputeVelocityNormalElement<2,3>(*p_element);

      std::array<double, 2> reference{1.0, 1.0};

      for (unsigned int i = 0; i < velocity.size(); i++) {
        KRATOS_EXPECT_NEAR(velocity(i), reference[i], 1e-7);
      }
    }

    // Checks the function ComputeVelocityUpperWakeElement
    KRATOS_TEST_CASE_IN_SUITE(ComputeVelocityUpperWakeElement, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);
      p_element->SetValue(WAKE, true);

      BoundedVector<double,3> distances = AssignDistances();
      p_element->SetValue(WAKE_ELEMENTAL_DISTANCES, distances);

      AssignPotentialsToWakeElement(p_element, distances);

      auto velocity = PotentialFlowUtilities::ComputeVelocityUpperWakeElement<2, 3>(*p_element);

      std::array<double, 2> reference{1.0, 1.0};

      for (unsigned int i = 0; i < velocity.size(); i++) {
        KRATOS_EXPECT_NEAR(velocity(i), reference[i], 1e-7);
      }
    }

    // Checks the function ComputeVelocityLowerWakeElement
    KRATOS_TEST_CASE_IN_SUITE(ComputeVelocityLowerWakeElement, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);
      p_element->SetValue(WAKE, true);

      BoundedVector<double,3> distances = AssignDistances();
      p_element->SetValue(WAKE_ELEMENTAL_DISTANCES, distances);

      AssignPotentialsToWakeElement(p_element, distances);

      auto velocity = PotentialFlowUtilities::ComputeVelocityLowerWakeElement<2, 3>(*p_element);

      std::array<double, 2> reference{1.0, 1.0};

      for (unsigned int i = 0; i < velocity.size(); i++) {
        KRATOS_EXPECT_NEAR(velocity(i), reference[i], 1e-7);
      }
    }

    // Checks the function ComputeVelocity
    KRATOS_TEST_CASE_IN_SUITE(ComputeVelocity, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);

      AssignPotentialsToNormalElement(p_element);

      auto velocity = PotentialFlowUtilities::ComputeVelocity<2,3>(*p_element);

      std::array<double, 2> reference{1.0, 1.0};

      for (unsigned int i = 0; i < velocity.size(); i++) {
        KRATOS_EXPECT_NEAR(velocity(i), reference[i], 1e-7);
      }
    }

    // Checks the function ComputeIncompressiblePressureCoefficient
    KRATOS_TEST_CASE_IN_SUITE(ComputeIncompressiblePressureCoefficient, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);

      AssignPotentialsToNormalElement(p_element);

      const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
      double pressure_coefficient = PotentialFlowUtilities::ComputeIncompressiblePressureCoefficient<2,3>(*p_element, r_current_process_info);

      KRATOS_EXPECT_NEAR(pressure_coefficient, 0.98, 1e-7);
    }

  } // namespace Testing
}  // namespace Kratos.
