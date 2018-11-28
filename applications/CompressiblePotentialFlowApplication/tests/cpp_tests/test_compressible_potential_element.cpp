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
		std::vector<ModelPart::IndexType> elemNodes1{ 1, 2, 3 };
		rModelPart.CreateNewElement("IncompressiblePotentialFlowElement2D3N", 1, elemNodes1, pProp);
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
		std::vector<ModelPart::IndexType> elemNodes1{ 1, 2, 3 };
		std::vector<ModelPart::IndexType> condNodes1{ 1, 2};
		std::vector<ModelPart::IndexType> condNodes2{ 2, 3};
		std::vector<ModelPart::IndexType> condNodes3{ 3, 1};
		rModelPart.CreateNewElement("IncompressibleStressesPotentialFlowElement2D3N", 1, elemNodes1, pProp);
		rModelPart.CreateNewCondition("IncompressiblePotentialWallCondition2D2N", 1, condNodes1, pProp);
		rModelPart.CreateNewCondition("IncompressiblePotentialWallCondition2D2N", 2, condNodes2, pProp);
		rModelPart.CreateNewCondition("IncompressiblePotentialWallCondition2D2N", 3, condNodes3, pProp);
    }
	void GenerateElementStressesHexagon(ModelPart& rModelPart)
    {
		// Variables addition
		rModelPart.AddNodalSolutionStepVariable(POSITIVE_POTENTIAL);
		rModelPart.AddNodalSolutionStepVariable(NEGATIVE_POTENTIAL);
		rModelPart.AddNodalSolutionStepVariable(WAKE_DISTANCE);

		// Set the element properties
		Properties::Pointer pProp = rModelPart.pGetProperties(0);

		// Geometry creation
		rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
		rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
		rModelPart.CreateNewNode(3, 0.65, 0.65 ,0.0);
		rModelPart.CreateNewNode(4, -0.65, 0.65, 0.0);
		rModelPart.CreateNewNode(5, -1.0, 0.0, 0.0);
		rModelPart.CreateNewNode(6, -0.65, -0.65, 0.0);
		rModelPart.CreateNewNode(7, 0.65, -0.65, 0.0);
		std::vector<ModelPart::IndexType> elemNodes1{ 1, 2, 3 };
		std::vector<ModelPart::IndexType> elemNodes2{ 1, 3, 4 };
		std::vector<ModelPart::IndexType> elemNodes3{ 1, 4, 5 };
		std::vector<ModelPart::IndexType> elemNodes4{ 1, 5, 6 };
		std::vector<ModelPart::IndexType> elemNodes5{ 1, 6, 7 };
		std::vector<ModelPart::IndexType> elemNodes6{ 1, 7, 2 };
		std::vector<ModelPart::IndexType> condNodes1{ 2, 3};
		std::vector<ModelPart::IndexType> condNodes2{ 3, 4};
		std::vector<ModelPart::IndexType> condNodes3{ 4, 5};
		std::vector<ModelPart::IndexType> condNodes4{ 5, 6};
		std::vector<ModelPart::IndexType> condNodes5{ 6, 7};
		std::vector<ModelPart::IndexType> condNodes6{ 7, 2};
		rModelPart.CreateNewElement("IncompressibleStressesPotentialFlowElement2D3N", 1, elemNodes1, pProp);
		rModelPart.CreateNewElement("IncompressibleStressesPotentialFlowElement2D3N", 2, elemNodes2, pProp);
		rModelPart.CreateNewElement("IncompressibleStressesPotentialFlowElement2D3N", 3, elemNodes3, pProp);
		rModelPart.CreateNewElement("IncompressibleStressesPotentialFlowElement2D3N", 4, elemNodes4, pProp);
		rModelPart.CreateNewElement("IncompressibleStressesPotentialFlowElement2D3N", 5, elemNodes5, pProp);
		rModelPart.CreateNewElement("IncompressibleStressesPotentialFlowElement2D3N", 6, elemNodes6, pProp);
		rModelPart.CreateNewCondition("IncompressiblePotentialWallCondition2D2N", 1, condNodes1, pProp);
		rModelPart.CreateNewCondition("IncompressiblePotentialWallCondition2D2N", 2, condNodes2, pProp);
		rModelPart.CreateNewCondition("IncompressiblePotentialWallCondition2D2N", 3, condNodes3, pProp);	
		rModelPart.CreateNewCondition("IncompressiblePotentialWallCondition2D2N", 4, condNodes4, pProp);
		rModelPart.CreateNewCondition("IncompressiblePotentialWallCondition2D2N", 5, condNodes5, pProp);
		rModelPart.CreateNewCondition("IncompressiblePotentialWallCondition2D2N", 6, condNodes6, pProp);
    }

    /** Checks the CompressiblePotentialFlowElement element.
     * Checks the LHS and RHS computation.
     */
    KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElement_CalculateLocalSystemStresses, CompressiblePotentialApplicationFastSuite)
    {
      	std::cout<<std::endl;
		Model this_model;
		ModelPart& model_part = this_model.CreateModelPart("Main", 3);
		GenerateElementStresses(model_part);

		Element::Pointer pElement = model_part.pGetElement(1);
		model_part.GetProcessInfo().SetValue(INITIAL_PENALTY,2.0);
		// Define the nodal values
		Vector potential(3);
		Vector distances(3);
    	int NumNodes = 3;
		for (unsigned int i = 0; i < 3; i++){
			potential(i) = pElement -> GetGeometry()[i].X()+pElement->GetGeometry()[i].Y();
		}

		distances(0) = 1.0;
		distances(1) = -1.0;
		distances(2) = -1.0;

		for (unsigned int i = 0; i < 3; i++){
      	if(distances[i] > 0)
      	 	pElement->GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL)= potential[i];      	
      	else
        	pElement->GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL)= potential[i];
		
      	pElement->GetGeometry()[i].FastGetSolutionStepValue(WAKE_DISTANCE) = distances(i);
    	}
		//negative part - sign is opposite to the previous case
		for (unsigned int i = 0; i < 3; i++){
			if(distances[i] < 0)
				pElement->GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL)= potential[i]+5;
			else
				pElement->GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL)= potential[i]+5;			
		}
		pElement -> GetValue(ELEMENTAL_DISTANCES) = distances;
		pElement -> Set(MARKER);

		// Compute RHS and LHS
		Vector RHS = ZeroVector(6);
    	Vector RHS_cond = ZeroVector(6);
		Matrix LHS = ZeroMatrix(6, 6);

		Vector vinf = ZeroVector(3);
		vinf(0)=1.0;
    	vinf(1)=+1.0;

		pElement->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

		for (unsigned int i_cond=1; i_cond<4;i_cond++){			
			Condition::Pointer pCondition = model_part.pGetCondition(i_cond);
			Matrix lhs = ZeroMatrix(2,2);
			Vector rhs = ZeroVector(2);
     	 	if ((i_cond ==1) || (i_cond == 3)){
				pCondition -> Set(STRUCTURE);
			}
			pCondition -> GetValue(VELOCITY_INFINITY) = vinf;
			pCondition -> CalculateLocalSystem(lhs,rhs,model_part.GetProcessInfo());
      		if (pCondition->Is(STRUCTURE)){        
        		int I_0 = (pCondition -> GetGeometry()[0].Id())-1;
        		int I_1 = (pCondition -> GetGeometry()[1].Id())-1;        
				RHS_cond(I_0) += rhs(0);
				RHS_cond(I_1) += rhs(1);
				RHS_cond(I_0+NumNodes) += rhs(2);
				RHS_cond(I_1+NumNodes) += rhs(3);        
			}else{
				for (unsigned int i_node=0;i_node<2;i_node++){
				int I = (pCondition -> GetGeometry()[i_node].Id())-1;
				if(distances[I] > 0)
					RHS_cond(I) += rhs(i_node);          
				else            
					RHS_cond(I+NumNodes) += rhs(i_node);    
				}
     		}
		}
		std::cout<< "Incompressible_stresses" << std::endl;
		KRATOS_WATCH(LHS);

		KRATOS_WATCH(RHS);
	
    	KRATOS_WATCH(RHS_cond);

		KRATOS_CHECK_NEAR(RHS(0), RHS_cond(0), 1e-7);
		KRATOS_CHECK_NEAR(RHS(1), RHS_cond(1), 1e-7);
		KRATOS_CHECK_NEAR(RHS(2), RHS_cond(2), 1e-7);
    }
    KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElement_CalculateLocalSystemStressesHexagon, CompressiblePotentialApplicationFastSuite)
    {
    	std::cout<<std::endl;
		Model this_model;
		ModelPart& model_part = this_model.CreateModelPart("Main", 3);
		GenerateElementStressesHexagon(model_part);
		Vector vinf = ZeroVector(3);
		vinf(0)=1.0;
		vinf(1)=-1.0;
		int number_of_nodes=model_part.NumberOfNodes();
		KRATOS_WATCH(number_of_nodes);
		Matrix LHS = ZeroMatrix(2*number_of_nodes,2*number_of_nodes);
		Vector RHS = ZeroVector(2*number_of_nodes);
		Vector RHS_cond = ZeroVector(2*number_of_nodes);
		model_part.GetProcessInfo().SetValue(INITIAL_PENALTY,2.0);
		unsigned int NumNodes = 3;
		Vector distances_all(model_part.NumberOfNodes()	);
		for (unsigned int i_node=0;i_node<model_part.NumberOfNodes();i_node++){
			if (i_node==0)
				distances_all(i_node)=-1.0;
			else
				distances_all(i_node)=1.0;
		} 

    	for (unsigned int i_elem=1;i_elem<model_part.NumberOfElements()+1;++i_elem){
			Element::Pointer pElement = model_part.pGetElement(i_elem);

			
			// Define the nodal values
			Vector potential(3);
			Vector distances(3);			
			for (unsigned int i_node=0;i_node<NumNodes;i_node++){
				int id=pElement->GetGeometry()[i_node].Id()-1;
				distances(i_node)=distances_all(id);
				potential(i_node) = pElement -> GetGeometry()[i_node].X()-pElement->GetGeometry()[i_node].Y();
			}

			for (unsigned int i = 0; i < NumNodes; i++){
				if(distances[i] > 0){
					pElement->GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL)= potential[i];
				}
				else{
					pElement->GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL)= potential[i];
				}
				pElement->GetGeometry()[i].FastGetSolutionStepValue(WAKE_DISTANCE) = distances(i);
			}

				//negative part - sign is opposite to the previous case
			for (unsigned int i = 0; i < NumNodes; i++){
				if(distances[i] < 0){
					pElement->GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL)= potential[i]+5;
				}
				else{
					pElement->GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL)= potential[i]+5;
				}
			}

			pElement -> GetValue(ELEMENTAL_DISTANCES) = distances;
			pElement -> Set(MARKER);

			Vector rhs = ZeroVector(6);			
			Matrix lhs = ZeroMatrix(6, 6);	
			pElement->CalculateLocalSystem(lhs, rhs, model_part.GetProcessInfo());
			// Assembling rhs
			for (unsigned int i_node=0;i_node<NumNodes;i_node++){
				int I=pElement->GetGeometry()[i_node].Id()-1;
				if(distances_all[I] > 0){					
					RHS(I) += rhs(i_node);  
					RHS(I+number_of_nodes) += rhs(i_node+NumNodes);     
				}   
				else{  
		 			RHS(I+number_of_nodes) += rhs(i_node); 
					RHS(I) += rhs(i_node+NumNodes);
				}
			}
    	}
		for (unsigned int i_cond=1; i_cond<model_part.NumberOfConditions()+1;i_cond++){
				
				Condition::Pointer pCondition = model_part.pGetCondition(i_cond);
				Matrix lhs_cond = ZeroMatrix(2,2);
				Vector rhs_cond = ZeroVector(2);	
				pCondition -> GetValue(VELOCITY_INFINITY) = vinf;
				pCondition -> CalculateLocalSystem(lhs_cond,rhs_cond,model_part.GetProcessInfo());
				if (pCondition->Is(STRUCTURE)){				
					int I_0 = (pCondition -> GetGeometry()[0].Id())-1;
					int I_1 = (pCondition -> GetGeometry()[1].Id())-1;
					
					RHS_cond(I_0) += rhs_cond(0);
					RHS_cond(I_1) += rhs_cond(1);
					RHS_cond(I_0+NumNodes) += rhs_cond(2);
					RHS_cond(I_1+NumNodes) += rhs_cond(3);				
				}
				else{
					for (unsigned int i_node=0;i_node<2;i_node++){
						int I = (pCondition -> GetGeometry()[i_node].Id())-1;
						if(distances_all[I] > 0)
							RHS_cond(I) += rhs_cond(i_node);          
						else            
							RHS_cond(I+number_of_nodes) += rhs_cond(i_node);    
					}
				}
		}
		// std::cout<< "Hexagon" << std::endl;
		// KRATOS_WATCH(LHS);
		// KRATOS_WATCH(LHS_ref);
		// KRATOS_WATCH(RHS);
		// KRATOS_WATCH(RHS_ref);
    	// KRATOS_WATCH(RHS_cond);
		// Check the RHS values (the RHS is computed as the LHS x previous_solution, 
		// hence, it is assumed that if the RHS is correct, the LHS is correct as well)
		KRATOS_CHECK_NEAR(RHS(0), RHS_cond(0), 1e-7);
		KRATOS_CHECK_NEAR(RHS(1), RHS_cond(1), 1e-7);
		KRATOS_CHECK_NEAR(RHS(2), RHS_cond(2), 1e-7);
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
