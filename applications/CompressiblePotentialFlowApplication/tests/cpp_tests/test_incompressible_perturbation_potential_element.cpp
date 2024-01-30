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
#include "custom_elements/incompressible_perturbation_potential_flow_element.h"
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos {
  namespace Testing {

    typedef ModelPart::IndexType IndexType;
    typedef ModelPart::NodeIterator NodeIteratorType;

    void GenerateIncompressiblePerturbationElement(ModelPart& rModelPart)
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
      rModelPart.CreateNewElement("IncompressiblePerturbationPotentialFlowElement2D3N", 1, elemNodes, pElemProp);
    }

    void AssignPotentialsToNormalPerturbationElement(Element::Pointer pElement)
    {
      // Define the nodal values
      std::array<double, 3> potential{1.0, 2.0, 3.0};

      for (unsigned int i = 0; i < 3; i++)
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
    }

    void AssignPotentialsToWakePerturbationElement(Element::Pointer pElement, const array_1d<double, 3>& rDistances, const std::array<double, 6>& rPotential)
    {
      for (unsigned int i = 0; i < 3; i++){
        if (rDistances(i) > 0.0)
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = rPotential[i];
        else
            pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = rPotential[i];
    }
    for (unsigned int i = 0; i < 3; i++){
        if (rDistances(i) < 0.0)
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = rPotential[i+3];
        else
            pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = rPotential[i+3];
    }
    }

    BoundedVector<double,3> AssignDistancesToPerturbationElement()
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
    KRATOS_TEST_CASE_IN_SUITE(IncompressiblePerturbationPotentialFlowElementRHS, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateIncompressiblePerturbationElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);

      AssignPotentialsToNormalPerturbationElement(p_element);

      // Compute RHS
      Vector RHS = ZeroVector(3);

      const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
      p_element->CalculateRightHandSide(RHS, r_current_process_info);

      // Check the RHS values
      std::vector<double> reference{5.5, -5, -0.5};

      KRATOS_EXPECT_VECTOR_NEAR(RHS, reference, 1e-6);
    }

    /** Checks the IncompressiblePotentialFlowElement element.
     * Checks the LHS and RHS computation.
     */
    KRATOS_TEST_CASE_IN_SUITE(IncompressiblePerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateIncompressiblePerturbationElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);

      AssignPotentialsToNormalPerturbationElement(p_element);

      // Compute LHS
      Matrix LHS = ZeroMatrix(3,3);

      const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
      p_element->CalculateLeftHandSide(LHS, r_current_process_info);

      // Check the LHS values
      std::array<double, 9> reference{0.5, -0.5, 0, -0.5, 1, -0.5, 0.0, -0.5, 0.5};

      for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
          KRATOS_EXPECT_NEAR(LHS(i,j), reference[3*i+j], 1e-6);
        }
      }
    }

    KRATOS_TEST_CASE_IN_SUITE(WakeIncompressiblePerturbationPotentialFlowElementCalculateRHS, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateIncompressiblePerturbationElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);

      BoundedVector<double,3> distances = AssignDistancesToPerturbationElement();

      p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
      p_element->GetValue(WAKE) = true;

      const std::array<double, 6> potential{1.5348, 2.31532, 3.2874, 4.1642, 6.8254, 5.174};
      AssignPotentialsToWakePerturbationElement(p_element, distances, potential);

      // Compute RHS
      Vector RHS = ZeroVector(6);

      const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
      p_element->CalculateRightHandSide(RHS, r_current_process_info);

      // Check the RHS values
      std::vector<double> reference{5.39026,2.252080000000001,-1.31174,0.9403400000000008,-7.1563,0.8256999999999999};

      KRATOS_EXPECT_VECTOR_NEAR(RHS, reference, 1e-13);
    }

    KRATOS_TEST_CASE_IN_SUITE(WakeIncompressiblePerturbationPotentialFlowElementCalculateLHS, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateIncompressiblePerturbationElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);

      BoundedVector<double,3> distances = AssignDistancesToPerturbationElement();

      p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
      p_element->GetValue(WAKE) = true;

      const std::array<double, 6> potential{1.0, 2.0, 3.0, 6.0, 7.0, 8.0};
      AssignPotentialsToWakePerturbationElement(p_element, distances, potential);

      // Compute LHS
      Matrix LHS = ZeroMatrix(6,6);

      const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
      p_element->CalculateLeftHandSide(LHS, r_current_process_info);

      // Check the LHS values
      std::array<double, 36> reference
          {0.5, -0.5, 0,   0,    0,   0,    -0.5, 1,   -0.5, 0.5, -1,   0.5,
           0,   -0.5, 0.5, -0,   0.5, -0.5, -0.5, 0.5, -0,   0.5, -0.5, 0,
           0,   0,    0,   -0.5, 1,   -0.5, 0,    0,   0,    0,   -0.5, 0.5};

      for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
          KRATOS_EXPECT_NEAR(LHS(i,j), reference[6*i+j], 1e-6);
        }
      }
    }


    /** Checks the IncompressiblePotentialFlowElement element.
 * Checks the EquationIdVector.
 */
    KRATOS_TEST_CASE_IN_SUITE(IncompressiblePerturbationPotentialFlowElementEquationIdVector, CompressiblePotentialApplicationFastSuite)
    {

      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateIncompressiblePerturbationElement(model_part);
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
    KRATOS_TEST_CASE_IN_SUITE(IncompressiblePerturbationPotentialFlowElementEquationIdVectorWake, CompressiblePotentialApplicationFastSuite)
    {

      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateIncompressiblePerturbationElement(model_part);
      Element::Pointer p_element = model_part.pGetElement(1);
      p_element->SetValue(WAKE, true);

      BoundedVector<double,3> distances = AssignDistancesToPerturbationElement();
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

  } // namespace Testing
}  // namespace Kratos.
