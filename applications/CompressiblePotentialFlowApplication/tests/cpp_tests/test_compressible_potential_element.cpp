//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-POSITIVE_POTENTIALhysics
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
#include "containers/model.h"
#include "includes/model_part.h"
#include "custom_elements/compressible_potential_flow_element.h"

namespace Kratos {
  namespace Testing {

    typedef ModelPart::IndexType IndexType;
    typedef ModelPart::NodeIterator NodeIteratorType;
    void CreateNodes(ModelPart& rModelPart){         
      rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
      rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
      rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    }

    void GenerateElement(ModelPart& rModelPart)
    {
      // Variables addition
      rModelPart.AddNodalSolutionStepVariable(POSITIVE_POTENTIAL);
      rModelPart.AddNodalSolutionStepVariable(NEGATIVE_POTENTIAL);

      // Set the element properties
      Properties::Pointer pProp = rModelPart.pGetProperties(0);

      // Geometry creation
      CreateNodes(rModelPart);
      std::vector<ModelPart::IndexType> elemNodes{ 1, 2, 3 };
      rModelPart.CreateNewElement("IncompressiblePotentialFlowElement2D3N", 1, elemNodes, pProp);
      std::vector<ModelPart::IndexType> condNodes1{ 1, 2};
      std::vector<ModelPart::IndexType> condNodes2{ 2, 3};
      std::vector<ModelPart::IndexType> condNodes3{ 3, 1};
      rModelPart.CreateNewCondition("IncompressiblePotentialWallCondition2D2N", 1, condNodes1, pProp);
      rModelPart.CreateNewCondition("IncompressiblePotentialWallCondition2D2N", 2, condNodes2, pProp);
      rModelPart.CreateNewCondition("IncompressiblePotentialWallCondition2D2N", 3, condNodes3, pProp);
      
    }  
    void GenerateElementStresses(ModelPart& rModelPart)
    {
      // Variables addition
      rModelPart.AddNodalSolutionStepVariable(POSITIVE_POTENTIAL);
      rModelPart.AddNodalSolutionStepVariable(NEGATIVE_POTENTIAL);
	  rModelPart.AddNodalSolutionStepVariable(WAKE_DISTANCE);

      // Set the element properties
      Properties::Pointer pProp = rModelPart.pGetProperties(0);

      // Geometry creation
      CreateNodes(rModelPart);
      std::vector<ModelPart::IndexType> elemNodes{ 1, 2, 3 };
      std::vector<ModelPart::IndexType> condNodes1{ 1, 2};
      std::vector<ModelPart::IndexType> condNodes2{ 2, 3};
      std::vector<ModelPart::IndexType> condNodes3{ 3, 1};
      rModelPart.CreateNewElement("IncompressibleStressesPotentialFlowElement2D3N", 1, elemNodes, pProp);
      rModelPart.CreateNewCondition("IncompressiblePotentialWallCondition2D2N", 1, condNodes1, pProp);
      rModelPart.CreateNewCondition("IncompressiblePotentialWallCondition2D2N", 2, condNodes2, pProp);
      rModelPart.CreateNewCondition("IncompressiblePotentialWallCondition2D2N", 3, condNodes3, pProp);
    }

    /** Checks the CompressiblePotentialFlowElement element.
     * Checks the LHS and RHS computation.
     */
    KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElement_CalculateLocalSystemStresses, CompressiblePotentialApplicationFastSuite)
    {
      std::cout<<std::endl;
		Model this_model;
		ModelPart& model_part = this_model.CreateModelPart("Main", 3);
		ModelPart& model_part_ref = this_model.CreateModelPart("Main_ref", 3);
		//ModelPart model_part("Main");
		GenerateElementStresses(model_part);
		GenerateElement(model_part_ref);
		Element::Pointer pElement = model_part.pGetElement(1);
		Element::Pointer pElement_ref = model_part_ref.pGetElement(1);
		model_part.GetProcessInfo().SetValue(INITIAL_PENALTY,2.0);
		// Define the nodal values
		Vector potential(3);
		Vector distances(3);
    int NumNodes = 3;
    for (unsigned int i = 0; i < 3; i++){
      potential(i) = pElement -> GetGeometry()[i].X()+pElement->GetGeometry()[i].Y();
    }

		distances(0) = 4.0;
		distances(1) = -4.0;
		distances(2) = -1.0;

		for (unsigned int i = 0; i < 3; i++){
      if(distances[i] > 0){
        pElement->GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL)= potential[i];
        pElement_ref->GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL)= potential[i];
      }
      else{
        pElement->GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL)= potential[i];
				pElement_ref->GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL)= potential[i];
		  }
      pElement->GetGeometry()[i].FastGetSolutionStepValue(WAKE_DISTANCE) = distances(i);
    }

        //negative part - sign is opposite to the previous case
    for (unsigned int i = 0; i < 3; i++){
      if(distances[i] < 0){
        pElement->GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL)= potential[i]+5;
        pElement_ref->GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL)= potential[i]+5;
      }
      else{
        pElement->GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL)= potential[i]+5;
				pElement_ref->GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL)= potential[i]+5;
			}
    }

		// for (unsigned int i = 0; i < 3; i++){
		// 	pElement->GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL) = potential(i);
		// 	pElement_ref->GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL) = potential(i);
		// 	pElement->GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL) = potential(i)+5;
		// 	pElement_ref->GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL) = potential(i)+5;
		// 	pElement->GetGeometry()[i].FastGetSolutionStepValue(WAKE_DISTANCE) = distances(i);
		// }
		pElement -> GetValue(ELEMENTAL_DISTANCES) = distances;
		pElement -> Set(MARKER);
		pElement_ref -> GetValue(ELEMENTAL_DISTANCES) = distances;
		pElement_ref -> Set(MARKER);
		// Compute RHS and LHS
		Vector RHS = ZeroVector(6);
    Vector RHS_cond = ZeroVector(6);
		Matrix LHS = ZeroMatrix(6, 6);
		Vector RHS_ref = ZeroVector(6);
		Vector vinf = ZeroVector(3);
		vinf(0)=1.0;
    vinf(1)=+1.0;
		Matrix LHS_ref = ZeroMatrix(6, 6);
		pElement->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());
		pElement_ref->CalculateLocalSystem(LHS_ref, RHS_ref, model_part_ref.GetProcessInfo());
		for (unsigned int i_cond=1; i_cond<4;i_cond++){
			
			Condition::Pointer pCondition = model_part.pGetCondition(i_cond);
			Matrix lhs = ZeroMatrix(2,2);
			Vector rhs = ZeroVector(2);
      if ((i_cond ==1) || (i_cond == 3))
        pCondition -> Set(STRUCTURE);
			pCondition -> GetValue(VELOCITY_INFINITY) = vinf;
			pCondition -> CalculateLocalSystem(lhs,rhs,model_part.GetProcessInfo());
      if (pCondition->Is(STRUCTURE)){
        for (unsigned int i_node=0; i_node<2;i_node++){
          int I = (pCondition -> GetGeometry()[i_node].Id())-1;
          if(distances[I] > 0){
            RHS_cond(I) += rhs(0);
            RHS_cond(I+NumNodes) += rhs(2);
          }else{
            RHS_cond(I) += rhs(1);
            RHS_cond(I+NumNodes) += rhs(3);
          }
        }
      }else{
        for (unsigned int i_node=0;i_node<2;i_node++){
          int I = (pCondition -> GetGeometry()[i_node].Id())-1;
          if(distances[I] > 0)
            RHS_cond(I) += rhs(i_node);          
          else            
            RHS_cond(I+NumNodes) += rhs(i_node);    
        }
      }
      std::cout<<"rhs"<<rhs<<std::endl;
		}
		std::cout<< "Incompressible_stresses" << std::endl;
		KRATOS_WATCH(LHS);
		KRATOS_WATCH(LHS_ref);
		KRATOS_WATCH(RHS);
		KRATOS_WATCH(RHS_ref);
    KRATOS_WATCH(RHS_cond);
		// Check the RHS values (the RHS is computed as the LHS x previous_solution, 
		// hence, it is assumed that if the RHS is correct, the LHS is correct as well)
		KRATOS_CHECK_NEAR(RHS(0), RHS_ref(0), 1e-7);
		KRATOS_CHECK_NEAR(RHS(1), RHS_ref(1), 1e-7);
		KRATOS_CHECK_NEAR(RHS(2), RHS_ref(2), 1e-7);
    }

    KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElement_CalculateLocalSystem, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);
      //ModelPart model_part("Main");
      GenerateElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);

      // Define the nodal values
      Vector potential(3);

      for (unsigned int i = 0; i < 3; i++){
        potential(i) = pElement -> GetGeometry()[i].X()+pElement->GetGeometry()[i].Y();
      }
      Vector vinf = ZeroVector(3);
		  vinf(0)=1.0;
      vinf(1)=+1.0;
      for (unsigned int i = 0; i < 3; i++)
        pElement->GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL) = potential(i);

      // Compute RHS and LHS
      Vector RHS = ZeroVector(3);
      Vector RHS_cond = ZeroVector(3);

      Matrix LHS = ZeroMatrix(3, 3);

      pElement->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

      for (unsigned int i_cond=1; i_cond<4;i_cond++){			
        Condition::Pointer pCondition = model_part.pGetCondition(i_cond);
        Matrix lhs = ZeroMatrix(2,2);
        Vector rhs = ZeroVector(2);
        // if ((i_cond ==1) || (i_cond == 3))
          // pCondition -> Set(STRUCTURE);
        pCondition -> GetValue(VELOCITY_INFINITY) = vinf;
        pCondition -> CalculateLocalSystem(lhs,rhs,model_part.GetProcessInfo());
        for (unsigned int i_node=0; i_node<2;i_node++){
            int I = (pCondition -> GetGeometry()[i_node].Id())-1;
            RHS_cond(I) += rhs(i_node);
        }
        std::cout << rhs << std::endl;  
      }
      std::cout << "RHS: " << RHS << std::endl;
      std::cout << "RHS_cond: "<< RHS_cond<< std::endl;
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
        pElement->GetGeometry()[i].AddDof(POSITIVE_POTENTIAL);

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
        pElement->GetGeometry()[i].AddDof(POSITIVE_POTENTIAL);
        pElement->GetGeometry()[i].AddDof(NEGATIVE_POTENTIAL);
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
