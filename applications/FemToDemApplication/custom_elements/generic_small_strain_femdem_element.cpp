//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

// System includes

// External includes

// Project includes
#include "generic_small_strain_femdem_element.hpp"
#include "fem_to_dem_application_variables.h"
#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericSmallStrainFemDemElement<TDim, TyieldSurf>::GenericSmallStrainFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry)
	: SmallDisplacementElement(NewId, pGeometry)
{
	//DO NOT ADD DOFS HERE!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericSmallStrainFemDemElement<TDim, TyieldSurf>::GenericSmallStrainFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
	: SmallDisplacementElement(NewId, pGeometry, pProperties)
{
	// BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
	mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

	// Each component == Each edge
	mNonConvergedThresholds = ZeroVector(NumberOfEdges);   // Equivalent stress
	mThresholds = ZeroVector(NumberOfEdges); // Stress mThreshold on edge
	mDamages = ZeroVector(NumberOfEdges); // Converged mDamage on each edge
	mNonConvergedDamages = ZeroVector(NumberOfEdges); // mDamages on edges of "i" iteration
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericSmallStrainFemDemElement<TDim,TyieldSurf>::GenericSmallStrainFemDemElement(GenericSmallStrainFemDemElement const &rOther)
	: SmallDisplacementElement(rOther)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericSmallStrainFemDemElement<TDim, TyieldSurf> &GenericSmallStrainFemDemElement<TDim,TyieldSurf>::operator=(GenericSmallStrainFemDemElement const &rOther)
{
	SmallDisplacementElement::operator=(rOther);
	return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericSmallStrainFemDemElement<TDim,TyieldSurf>::Create(IndexType NewId, NodesArrayType const &rThisNodes, PropertiesType::Pointer pProperties) const
{
	return Element::Pointer(new GenericSmallStrainFemDemElement(NewId, GetGeometry().Create(rThisNodes), pProperties));
}

//************************************CLONE*******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericSmallStrainFemDemElement<TDim,TyieldSurf>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
{
	GenericSmallStrainFemDemElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
	return Element::Pointer(new GenericSmallStrainFemDemElement(NewElement));
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericSmallStrainFemDemElement<TDim,TyieldSurf>::~GenericSmallStrainFemDemElement()
{
}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::InitializeSolutionStep(
    ProcessInfo& rCurrentProcessInfo
    )
{
	this->ComputeEdgeNeighbours(rCurrentProcessInfo);
	// this->InitializeInternalVariablesAfterMapping();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericSmallStrainFemDemElement<3,0>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}
template<>
void GenericSmallStrainFemDemElement<3,1>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}
template<>
void GenericSmallStrainFemDemElement<3,2>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}
template<>
void GenericSmallStrainFemDemElement<3,3>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}
template<>
void GenericSmallStrainFemDemElement<3,4>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}
template<>
void GenericSmallStrainFemDemElement<3,5>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}
template<>
void GenericSmallStrainFemDemElement<3,6>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericSmallStrainFemDemElement<2,0>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}
template<> 
void GenericSmallStrainFemDemElement<2,1>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}
template<> 
void GenericSmallStrainFemDemElement<2,2>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}
template<> 
void GenericSmallStrainFemDemElement<2,3>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}
template<> 
void GenericSmallStrainFemDemElement<2,4>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}
template<> 
void GenericSmallStrainFemDemElement<2,5>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}
template<> 
void GenericSmallStrainFemDemElement<2,6>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}


/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::AuxComputeEdgeNeighbours(
    ProcessInfo& rCurrentProcessInfo
    )
{
	std::vector<std::vector<Element*>> edge_neighbours_container;
	Geometry<Node<3>>& r_nodes_current_element = this->GetGeometry();

	Node<3> &pNode0 = r_nodes_current_element[0];
	Node<3> &pNode1 = r_nodes_current_element[1];
	Node<3> &pNode2 = r_nodes_current_element[2];
	Node<3> &pNode3 = r_nodes_current_element[3];

	// Neighbour elements of each node of the current element
	GlobalPointersVector<Element>& neigh_node_0 = pNode0.GetValue(NEIGHBOUR_ELEMENTS);
	GlobalPointersVector<Element>& neigh_node_1 = pNode1.GetValue(NEIGHBOUR_ELEMENTS);
	GlobalPointersVector<Element>& neigh_node_2 = pNode2.GetValue(NEIGHBOUR_ELEMENTS);
	GlobalPointersVector<Element>& neigh_node_3 = pNode3.GetValue(NEIGHBOUR_ELEMENTS);

	// Nodal neighbours container
	std::vector<GlobalPointersVector<Element>> nodal_neighbours;
	nodal_neighbours.push_back(neigh_node_0);
	nodal_neighbours.push_back(neigh_node_1);
	nodal_neighbours.push_back(neigh_node_2);
	nodal_neighbours.push_back(neigh_node_3);

	// Aux indexes
	Matrix nodes_indexes = ZeroMatrix(6, 2);
	this->SetNodeIndexes(nodes_indexes);

	// Loop over EDGES to assign the elements that share that edge -> Fill mEdgeNeighboursContainer
	for (unsigned int edge = 0; edge < 6; edge++) {
		const int node_index_1 = nodes_indexes(edge, 0);
		const int node_index_2 = nodes_indexes(edge, 1);

		// Neigh elements of local node 1 and 2  //
		auto& r_neigh_of_node_1 = nodal_neighbours[node_index_1];
		auto& r_neigh_of_node_2 = nodal_neighbours[node_index_2];

		const int node_id_1 = r_nodes_current_element[node_index_1].Id();
		const int node_id_2 = r_nodes_current_element[node_index_2].Id();

		std::vector<Element*> edge_shared_elements_node_1;
		// Loop over neigh elements of the node 1
		for (unsigned int neigh_elem = 0; neigh_elem < r_neigh_of_node_1.size(); neigh_elem++) {
			// Nodes of the neigh element
			Geometry<Node<3>>& r_nodes_neigh_elem = r_neigh_of_node_1[neigh_elem].GetGeometry();

			// Loop over the nodes of the neigh element
			for (unsigned int neigh_elem_node = 0; neigh_elem_node < 4; neigh_elem_node++) {
				const int neigh_element_node_id = r_nodes_neigh_elem[neigh_elem_node].Id();

				if (neigh_element_node_id == node_id_2 && this->Id() != r_neigh_of_node_1[neigh_elem].Id()) {
					edge_shared_elements_node_1.push_back(&r_neigh_of_node_1[neigh_elem]); // ( [] returns an Element object!!)
				}
			}
		}

		std::vector<Element *> edge_shared_elements_node_2;
		// Loop over neigh elements of the node 2
		for (unsigned int neigh_elem = 0; neigh_elem < r_neigh_of_node_2.size(); neigh_elem++) {
			// Nodes of the neigh element
			Geometry<Node<3>> &r_nodes_neigh_elem = r_neigh_of_node_2[neigh_elem].GetGeometry();

			// Loop over the nodes of the neigh element
			for (unsigned int neigh_elem_node = 0; neigh_elem_node < 4; neigh_elem_node++) {
				const int neigh_element_node_id = r_nodes_neigh_elem[neigh_elem_node].Id();

				if (neigh_element_node_id == node_id_1 && this->Id() != r_neigh_of_node_2[neigh_elem].Id()) {
					edge_shared_elements_node_2.push_back(&r_neigh_of_node_2[neigh_elem]);
				}
			}
		}
		// Let's create the vector of neighbour elements for this edge
		std::vector<Element *> edge_shared_elements = edge_shared_elements_node_1;
		// Add the neigh elements from the node 2
		for (unsigned int i = 0; i < edge_shared_elements_node_2.size(); i++) {
			int aux = 0;

			for (unsigned int j = 0; j < edge_shared_elements.size(); j++) {
				if (edge_shared_elements_node_2[i]->Id() == edge_shared_elements[j]->Id())
					aux++;
			}
			if (aux == 0)
				edge_shared_elements.push_back(edge_shared_elements_node_2[i]);
		}
		edge_neighbours_container.push_back(edge_shared_elements);
	} // End loop edges

	// Storages the information inside the element
	this->SaveEdgeNeighboursContainer(edge_neighbours_container);
}
/***********************************************************************************/
/***********************************************************************************/

template class GenericSmallStrainFemDemElement<2,0>;
template class GenericSmallStrainFemDemElement<2,1>;
template class GenericSmallStrainFemDemElement<2,2>;
template class GenericSmallStrainFemDemElement<2,3>;
template class GenericSmallStrainFemDemElement<2,4>;
template class GenericSmallStrainFemDemElement<2,5>;
template class GenericSmallStrainFemDemElement<2,6>;
template class GenericSmallStrainFemDemElement<3,0>;
template class GenericSmallStrainFemDemElement<3,1>;
template class GenericSmallStrainFemDemElement<3,2>;
template class GenericSmallStrainFemDemElement<3,3>;
template class GenericSmallStrainFemDemElement<3,4>;
template class GenericSmallStrainFemDemElement<3,5>;
template class GenericSmallStrainFemDemElement<3,6>;
} // namespace Kratos