//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez
//
//

// System includes
#include <set>

// External includes
#include "containers/model.h"

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"
#include "custom_elements/compressible_potential_flow_element.h"

namespace Kratos {
  namespace Testing {

    typedef ModelPart::IndexType IndexType;
    typedef ModelPart::NodeIterator NodeIteratorType;

    void GenerateElement(ModelPart& rModelPart)
    {
      // Variables addition
      rModelPart.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
      rModelPart.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE);

      // Set the element properties
      Properties::Pointer pElemProp = rModelPart.pGetProperties(0);

      // Geometry creation
      rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
      rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
      rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
      std::vector<ModelPart::IndexType> elemNodes{ 1, 2, 3 };
      rModelPart.CreateNewElement("CompressiblePotentialFlowElement2D3N", 1, elemNodes, pElemProp);
    }

    /** Checks the CompressiblePotentialFlowElement element.
     * Checks the LHS and RHS computation.
     */
    KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElement_CalculateLocalSystem, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);
      //ModelPart model_part("Main");
      GenerateElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);

      // Define the nodal values
      Vector potential(3);
      potential(0) = 1.0;
      potential(1) = 2.0;
      potential(2) = 3.0;

      for (unsigned int i = 0; i < 3; i++)
        pElement->GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE) = potential(i);

      // Compute RHS and LHS
      Vector RHS = ZeroVector(3);
      Matrix LHS = ZeroMatrix(3, 3);

      pElement->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

      // Check the RHS values (the RHS is computed as the LHS x previous_solution, 
      // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
      KRATOS_CHECK_NEAR(RHS(0), 0.5, 1e-7);
      KRATOS_CHECK_NEAR(RHS(1), 0.0, 1e-7);
      KRATOS_CHECK_NEAR(RHS(2), -0.5, 1e-7);
    }

    /** Checks the CompressiblePotentialFlowElement element.
 * Checks the EquationIdVector.
 */
    KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElement_EquationIdVector, CompressiblePotentialApplicationFastSuite)
    {

      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);
      //ModelPart model_part("Main");
      GenerateElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);

      for (unsigned int i = 0; i < 3; i++)
        pElement->GetGeometry()[i].AddDof(POSITIVE_FACE_PRESSURE);

      Element::DofsVectorType ElementalDofList;
      pElement->GetDofList(ElementalDofList, model_part.GetProcessInfo());

      for (int i = 0; i < 3; i++)
        ElementalDofList[i]->SetEquationId(i);

      Element::EquationIdVectorType EquationIdVector;
      pElement->EquationIdVector(EquationIdVector, model_part.GetProcessInfo());

      // Check the EquationIdVector values
      KRATOS_CHECK(EquationIdVector[0] == 0);
      KRATOS_CHECK(EquationIdVector[1] == 1);
      KRATOS_CHECK(EquationIdVector[2] == 2);
    }

    /** Checks the CompressiblePotentialFlowElement element.
 * Checks the EquationIdVector for the Wake.
 */
    KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElement_EquationIdVector_Wake, CompressiblePotentialApplicationFastSuite)
    {

      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);
      //ModelPart model_part("Main");
      GenerateElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);
      pElement->Set(MARKER, true);

      array_1d<double, 3> distances;
      distances[0] = -0.5;
      distances[1] = -0.5;
      distances[2] = 0.5;
      pElement->SetValue(ELEMENTAL_DISTANCES, distances);

      for (unsigned int i = 0; i < 3; i++) {
        pElement->GetGeometry()[i].AddDof(POSITIVE_FACE_PRESSURE);
        pElement->GetGeometry()[i].AddDof(NEGATIVE_FACE_PRESSURE);
      }

      Element::DofsVectorType ElementalDofList;
      pElement->GetDofList(ElementalDofList, model_part.GetProcessInfo());

      for (int i = 0; i < 6; i++)
        ElementalDofList[i]->SetEquationId(i);

      Element::EquationIdVectorType EquationIdVector;
      pElement->EquationIdVector(EquationIdVector, model_part.GetProcessInfo());

      //Check the EquationIdVector values
      KRATOS_CHECK(EquationIdVector[0] == 0);
      KRATOS_CHECK(EquationIdVector[1] == 1);
      KRATOS_CHECK(EquationIdVector[2] == 2);
      KRATOS_CHECK(EquationIdVector[3] == 3);
      KRATOS_CHECK(EquationIdVector[4] == 4);
      KRATOS_CHECK(EquationIdVector[5] == 5);
    }
  } // namespace Testing
}  // namespace Kratos.
