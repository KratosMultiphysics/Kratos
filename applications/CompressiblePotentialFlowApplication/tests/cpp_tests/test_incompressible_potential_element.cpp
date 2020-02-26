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
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_elements/incompressible_potential_flow_element.h"
#include "custom_elements/embedded_incompressible_potential_flow_element.h"
#include "custom_utilities/potential_flow_utilities.h"

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
      Element::Pointer pElement = model_part.pGetElement(1);

      AssignPotentialsToNormalElement(pElement);

      // Compute RHS and LHS
      Vector RHS = ZeroVector(3);
      Matrix LHS = ZeroMatrix(3, 3);

      pElement->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

      // Check the RHS values (the RHS is computed as the LHS x previous_solution,
      // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
      std::array<double, 3> reference{0.5, 0.0, -0.5};

      for (unsigned int i = 0; i < RHS.size(); i++) {
        KRATOS_CHECK_NEAR(RHS(i), reference[i], 1e-6);
      }
    }

    KRATOS_TEST_CASE_IN_SUITE(IncompressiblePotentialFlowElementCalculateLocalSystemWake, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);

      BoundedVector<double,3> distances = AssignDistances();

      pElement->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
      pElement->GetValue(WAKE) = true;

      AssignPotentialsToWakeElement(pElement, distances);

      // Compute RHS and LHS
      Vector RHS = ZeroVector(6);
      Matrix LHS = ZeroMatrix(6, 6);

      pElement->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

      // Check the RHS values (the RHS is computed as the LHS x previous_solution,
      // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
      std::array<double, 6> reference{0.5, 0.0, 0.0, 0.0, 0.0, -0.5};

      for (unsigned int i = 0; i < RHS.size(); i++) {
        KRATOS_CHECK_NEAR(RHS(i), reference[i], 1e-6);
      }
    }

    KRATOS_TEST_CASE_IN_SUITE(EmbeddedIncompressiblePotentialFlowElementCalculateLocalSystem, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateEmbeddedElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);

      AssignPotentialsToNormalElement(pElement);

      // Define the distance values
      Vector level_set(3);
      level_set(0) = 1.0;
      level_set(1) = -1.0;
      level_set(2) = -1.0;

      for (unsigned int i = 0; i < 3; i++){
        pElement->GetGeometry()[i].FastGetSolutionStepValue(GEOMETRY_DISTANCE) = level_set(i);
      }

      // Compute RHS and LHS
      Vector RHS = ZeroVector(3);
      Matrix LHS = ZeroMatrix(3, 3);

      pElement->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

      // Check the RHS values (the RHS is computed as the LHS x previous_solution,
      // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
      std::array<double, 3> reference{0.125, 0.0, -0.125};

      for (unsigned int i = 0; i < RHS.size(); i++) {
        KRATOS_CHECK_NEAR(RHS(i), reference[i], 1e-6);
      }
    }


    /** Checks the IncompressiblePotentialFlowElement element.
 * Checks the EquationIdVector.
 */
    KRATOS_TEST_CASE_IN_SUITE(IncompressiblePotentialFlowElementEquationIdVector, CompressiblePotentialApplicationFastSuite)
    {

      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);

      for (unsigned int i = 0; i < 3; i++)
        pElement->GetGeometry()[i].AddDof(VELOCITY_POTENTIAL);

      Element::DofsVectorType ElementalDofList;
      pElement->GetDofList(ElementalDofList, model_part.GetProcessInfo());

      for (int i = 0; i < 3; i++)
        ElementalDofList[i]->SetEquationId(i);

      Element::EquationIdVectorType EquationIdVector;
      pElement->EquationIdVector(EquationIdVector, model_part.GetProcessInfo());

      // Check the EquationIdVector values
      for (unsigned int i = 0; i < EquationIdVector.size(); i++) {
        KRATOS_CHECK(EquationIdVector[i] == i);
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
      Element::Pointer pElement = model_part.pGetElement(1);
      pElement->SetValue(WAKE, true);

      BoundedVector<double,3> distances = AssignDistances();
      pElement->SetValue(WAKE_ELEMENTAL_DISTANCES, distances);

      for (unsigned int i = 0; i < 3; i++) {
        pElement->GetGeometry()[i].AddDof(VELOCITY_POTENTIAL);
        pElement->GetGeometry()[i].AddDof(AUXILIARY_VELOCITY_POTENTIAL);
      }

      Element::DofsVectorType ElementalDofList;
      pElement->GetDofList(ElementalDofList, model_part.GetProcessInfo());

      for (int i = 0; i < 6; i++)
        ElementalDofList[i]->SetEquationId(i);

      Element::EquationIdVectorType EquationIdVector;
      pElement->EquationIdVector(EquationIdVector, model_part.GetProcessInfo());

      //Check the EquationIdVector values
      for (unsigned int i = 0; i < EquationIdVector.size(); i++) {
        KRATOS_CHECK(EquationIdVector[i] == i);
      }
    }

    // Checks the function GetWakeDistances
    KRATOS_TEST_CASE_IN_SUITE(GetWakeDistances, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);
      pElement->SetValue(WAKE, true);

      BoundedVector<double,3> distances = AssignDistances();
      pElement->SetValue(WAKE_ELEMENTAL_DISTANCES, distances);

      const auto returned_distances = PotentialFlowUtilities::GetWakeDistances<2, 3>(*pElement);

      std::array<double, 3> reference{1.0, -1.0, -1.0};

      for (unsigned int i = 0; i < returned_distances.size(); i++) {
        KRATOS_CHECK_NEAR(returned_distances(i), reference[i], 1e-7);
      }
    }

    // Checks the function GetPotentialOnNormalElement
    KRATOS_TEST_CASE_IN_SUITE(GetPotentialOnNormalElement, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);

      AssignPotentialsToNormalElement(pElement);

      auto potentials = PotentialFlowUtilities::GetPotentialOnNormalElement<2,3>(*pElement);

      std::array<double, 3> reference{1.0, 2.0, 3.0};

      for (unsigned int i = 0; i < potentials.size(); i++) {
        KRATOS_CHECK_NEAR(potentials(i), reference[i], 1e-7);
      }
    }

    // Checks the function GetPotentialOnUpperWakeElement
    KRATOS_TEST_CASE_IN_SUITE(GetPotentialOnUpperWakeElement, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);
      pElement->SetValue(WAKE, true);

      BoundedVector<double,3> distances = AssignDistances();

      AssignPotentialsToWakeElement(pElement, distances);

      array_1d<double, 3> potentials =
          PotentialFlowUtilities::GetPotentialOnUpperWakeElement<2, 3>(*pElement, distances);

      std::array<double, 3> reference{1.0, 2.0, 3.0};

      for (unsigned int i = 0; i < potentials.size(); i++) {
        KRATOS_CHECK_NEAR(potentials(i), reference[i], 1e-7);
      }
    }

    // Checks the function GetPotentialOnLowerWakeElement
    KRATOS_TEST_CASE_IN_SUITE(GetPotentialOnLowerWakeElement, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);
      pElement->SetValue(WAKE, true);

      BoundedVector<double,3> distances = AssignDistances();

      AssignPotentialsToWakeElement(pElement, distances);

      array_1d<double, 3> potentials =
          PotentialFlowUtilities::GetPotentialOnLowerWakeElement<2, 3>(*pElement, distances);

      std::array<double, 3> reference{6.0, 7.0, 8.0};

      for (unsigned int i = 0; i < potentials.size(); i++) {
        KRATOS_CHECK_NEAR(potentials(i), reference[i], 1e-7);
      }
    }

    // Checks the function GetPotentialOnWakeElement
    KRATOS_TEST_CASE_IN_SUITE(GetPotentialOnWakeElement, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);
      pElement->SetValue(WAKE, true);

      BoundedVector<double,3> distances = AssignDistances();

      AssignPotentialsToWakeElement(pElement, distances);

      BoundedVector<double, 2 * 3> potentials =
          PotentialFlowUtilities::GetPotentialOnWakeElement<2, 3>(*pElement, distances);

      std::array<double, 6> reference{1.0, 2.0, 3.0, 6.0, 7.0, 8.0};

      for (unsigned int i = 0; i < potentials.size(); i++) {
        KRATOS_CHECK_NEAR(potentials(i), reference[i], 1e-7);
      }
    }

    // Checks the function ComputeVelocityNormalElement
    KRATOS_TEST_CASE_IN_SUITE(ComputeVelocityNormalElement, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);

      AssignPotentialsToNormalElement(pElement);

      auto velocity = PotentialFlowUtilities::ComputeVelocityNormalElement<2,3>(*pElement);

      std::array<double, 2> reference{1.0, 1.0};

      for (unsigned int i = 0; i < velocity.size(); i++) {
        KRATOS_CHECK_NEAR(velocity(i), reference[i], 1e-7);
      }
    }

    // Checks the function ComputeVelocityUpperWakeElement
    KRATOS_TEST_CASE_IN_SUITE(ComputeVelocityUpperWakeElement, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);
      pElement->SetValue(WAKE, true);

      BoundedVector<double,3> distances = AssignDistances();
      pElement->SetValue(WAKE_ELEMENTAL_DISTANCES, distances);

      AssignPotentialsToWakeElement(pElement, distances);

      auto velocity = PotentialFlowUtilities::ComputeVelocityUpperWakeElement<2, 3>(*pElement);

      std::array<double, 2> reference{1.0, 1.0};

      for (unsigned int i = 0; i < velocity.size(); i++) {
        KRATOS_CHECK_NEAR(velocity(i), reference[i], 1e-7);
      }
    }

    // Checks the function ComputeVelocityLowerWakeElement
    KRATOS_TEST_CASE_IN_SUITE(ComputeVelocityLowerWakeElement, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);
      pElement->SetValue(WAKE, true);

      BoundedVector<double,3> distances = AssignDistances();
      pElement->SetValue(WAKE_ELEMENTAL_DISTANCES, distances);

      AssignPotentialsToWakeElement(pElement, distances);

      auto velocity = PotentialFlowUtilities::ComputeVelocityLowerWakeElement<2, 3>(*pElement);

      std::array<double, 2> reference{1.0, 1.0};

      for (unsigned int i = 0; i < velocity.size(); i++) {
        KRATOS_CHECK_NEAR(velocity(i), reference[i], 1e-7);
      }
    }

    // Checks the function ComputeVelocity
    KRATOS_TEST_CASE_IN_SUITE(ComputeVelocity, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);

      AssignPotentialsToNormalElement(pElement);

      auto velocity = PotentialFlowUtilities::ComputeVelocity<2,3>(*pElement);

      std::array<double, 2> reference{1.0, 1.0};

      for (unsigned int i = 0; i < velocity.size(); i++) {
        KRATOS_CHECK_NEAR(velocity(i), reference[i], 1e-7);
      }
    }

    // Checks the function ComputeIncompressiblePressureCoefficient
    KRATOS_TEST_CASE_IN_SUITE(ComputeIncompressiblePressureCoefficient, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);

      AssignPotentialsToNormalElement(pElement);

      double pressure_coefficient = PotentialFlowUtilities::ComputeIncompressiblePressureCoefficient<2,3>(*pElement, model_part.GetProcessInfo());

      KRATOS_CHECK_NEAR(pressure_coefficient, 0.98, 1e-7);
    }

  } // namespace Testing
}  // namespace Kratos.
