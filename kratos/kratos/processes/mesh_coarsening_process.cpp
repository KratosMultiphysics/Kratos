//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   
//

// Project includes
#include "processes/mesh_coarsening_process.h"

namespace Kratos
{

	KRATOS_CREATE_LOCAL_FLAG(MeshCoarseningProcess, COARSE_MESH_NODE, 0);

	MeshCoarseningProcess::MeshCoarseningProcess(ModelPart& rModelPart)
	: Process(), mrModelPart(rModelPart){

	}

	MeshCoarseningProcess:: ~MeshCoarseningProcess() {

	}

	void MeshCoarseningProcess::Execute() {
		SelectCoarseMeshNodes();
		CollapseNodes();

		// TO_ERASE flag is already set in CollapseNodes method
		mrModelPart.RemoveElements();
		mrModelPart.RemoveNodes();
	}

	std::string MeshCoarseningProcess::Info() const {
		return "MeshCoarseningProcess";
	}

	void MeshCoarseningProcess::PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
	}

	void MeshCoarseningProcess::PrintData(std::ostream& rOStream) const {

	}

	void MeshCoarseningProcess::SelectCoarseMeshNodes() {
		// initializing and selecting all boundary nodes as seed
		for (auto i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++) {
			i_node->GetValue(FATHER_NODES).clear();
			i_node->Set(COARSE_MESH_NODE, false);
			if (i_node->Is(BOUNDARY)) // We want to keep all the boundary nodes 
				i_node->Set(COARSE_MESH_NODE);
		}

		// Now selecting the seeds
		// A node without any seed in its neighbour is a seed
		for (ModelPart::NodeIterator i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++) {
			if (i_node->IsNot(COARSE_MESH_NODE)) {
				bool is_seed = true; // as a candidate
				auto& r_neighbours = i_node->GetValue(NEIGHBOUR_NODES);
				for (auto i_neighbour_node = r_neighbours.begin(); i_neighbour_node != r_neighbours.end(); i_neighbour_node++) {
					if ((i_neighbour_node->Is(COARSE_MESH_NODE))) {
						is_seed = false;
						break;
					}
				}
				if (is_seed)
					i_node->Set(COARSE_MESH_NODE);
			}
		}
	}

	void MeshCoarseningProcess::CollapseNodes()	{
		for (auto i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++)
			if (i_node->IsNot(COARSE_MESH_NODE))
				CollapseNode(*i_node);
	}

	void MeshCoarseningProcess::CollapseNode(Node<3>& rThisNode) {
		auto& r_neighbours = rThisNode.GetValue(NEIGHBOUR_NODES);
		auto i_coarse_node = r_neighbours.end();
		double current_quality = std::numeric_limits<double>::epsilon();

		for (auto i_neighbour_node = r_neighbours.begin(); i_neighbour_node != r_neighbours.end(); i_neighbour_node++) {
			if ((i_neighbour_node->Is(COARSE_MESH_NODE))) {
				auto quality = CalculateQualityIfNodeCollapses(rThisNode, *i_neighbour_node);
				if (quality > current_quality) {
					current_quality = quality;
					i_coarse_node = i_neighbour_node;
				}
			}
		}

		if (i_coarse_node == r_neighbours.end()) {
			KRATOS_ERROR << "The node " << rThisNode.Id() << " cannot be collapsed" << std::endl;
		}
		else {
			rThisNode.Set(TO_ERASE);
			auto& neighbour_elements = rThisNode.GetValue(NEIGHBOUR_ELEMENTS);
			for (auto i_element = neighbour_elements.begin(); i_element != neighbour_elements.end(); i_element++)
				if (ElementHas(*i_element, *i_coarse_node)) 
					i_element->Set(TO_ERASE);
				else 
					SwapElementNode(*i_element, rThisNode, i_coarse_node.base()->lock());
		}
	}

	double MeshCoarseningProcess::CalculateQualityIfNodeCollapses(Node<3>& rThisNode, Node<3> const& rCoarseNode) {
		Point<3> original_coordinates = rThisNode;
		rThisNode.Coordinates() = rCoarseNode.Coordinates();
		double min_quality = CalculateMinQualityOfNeighbourElements(rThisNode, rCoarseNode);
		rThisNode.Coordinates() = original_coordinates;
		return min_quality;
	}

	double MeshCoarseningProcess::CalculateMinQualityOfNeighbourElements(Node<3>& rThisNode, Node<3> const& rCoarseNode) {
		auto& neighbour_elements = rThisNode.GetValue(NEIGHBOUR_ELEMENTS);
		double min_quality = std::numeric_limits<double>::max();

		for (auto i_element = neighbour_elements.begin(); i_element != neighbour_elements.end(); i_element++) 
			if(!ElementHas(*i_element, rCoarseNode)) {
			double domain_size = i_element->GetGeometry().DomainSize();
			min_quality = std::min(min_quality, domain_size);
		}
		return min_quality;
	}

	bool MeshCoarseningProcess::ElementHas(Element& rElement, Node<3> const& rCoarseNode) {
		for (auto i_node = rElement.GetGeometry().begin(); i_node != rElement.GetGeometry().end(); i_node++)
			if (i_node->Id() == rCoarseNode.Id())
				return true;

		return false;
	}

	void MeshCoarseningProcess::SwapElementNode(Element& rElement, Node<3> const& rThisNode, Node<3>::Pointer pCoarseNode) {
		std::size_t number_of_nodes = rElement.GetGeometry().size();
		auto& geometry = rElement.GetGeometry();
		for (std::size_t i = 0; i < number_of_nodes; i++)
			if (geometry[i].Id() == rThisNode.Id())
				geometry(i) = pCoarseNode;
	}

}  // namespace Kratos.


