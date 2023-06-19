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
	           
// System includes


// External includes 


// Project includes
#include "processes/mesh_node_collapsing_process.h"
#include "includes/global_pointer_variables.h"


namespace Kratos
{
	KRATOS_CREATE_LOCAL_FLAG(MeshNodeCollapsingProcess, TO_COLLAPSE, 0);

	MeshNodeCollapsingProcess::MeshNodeCollapsingProcess(ModelPart& rModelPart) : mrModelPart(rModelPart) {

	}

	MeshNodeCollapsingProcess::~MeshNodeCollapsingProcess() {

	}

	void MeshNodeCollapsingProcess::Execute() {
		CollapseNodes();

		// TO_ERASE flag is already set in CollapseNodes method
		mrModelPart.RemoveElements();
		mrModelPart.RemoveNodes();

	}

	std::string MeshNodeCollapsingProcess::Info() const {
		return "MeshNodeCollapsingProcess";
	}

	/// Print information about this object.
	void MeshNodeCollapsingProcess::PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
	}

	/// Print object's data.
	void MeshNodeCollapsingProcess::PrintData(std::ostream& rOStream) const {

	}

	void MeshNodeCollapsingProcess::CollapseNodes() {
		for (auto i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++)
			if (i_node->Is(TO_COLLAPSE))
				CollapseNode(*i_node);
	}

	void MeshNodeCollapsingProcess::CollapseNode(Node& rThisNode) {
		auto& r_neighbours = rThisNode.GetValue(NEIGHBOUR_NODES);
		auto i_coarse_node = r_neighbours.end();

		//TetrahedraBall node_ball(rThisNode);
		//double current_quality = node_ball.CalculateMinQuality(Geometry<Node >::QualityCriteria::AVERAGE_LENGTH_VOLUME_RATIO);

		// initializing the min quality of the current mesh as the treshold and also check if there is an
		// element in the ball which is already set to be erased
		double current_quality = std::numeric_limits<double>::max();
		auto& neighbour_elements = rThisNode.GetValue(NEIGHBOUR_ELEMENTS);
		for (auto i_element = neighbour_elements.begin(); i_element != neighbour_elements.end(); i_element++) {
			if (i_element->IsNot(TO_ERASE)) {
				double domain_size = i_element->GetGeometry().DomainSize();
				current_quality = std::min(current_quality, domain_size);
			}
		}
		current_quality *= .1;

		for (auto i_neighbour_node = r_neighbours.begin(); i_neighbour_node != r_neighbours.end(); i_neighbour_node++) {
			if ((i_neighbour_node->IsNot(TO_COLLAPSE))) {
				auto quality = CalculateQualityIfNodeCollapses(rThisNode, *i_neighbour_node);
				if (quality > current_quality) {
					current_quality = quality;
					i_coarse_node = i_neighbour_node;
				}
			}
		}

		if (i_coarse_node != r_neighbours.end()) {
			rThisNode.Set(TO_ERASE);
			auto& neighbour_elements = rThisNode.GetValue(NEIGHBOUR_ELEMENTS);
			for (auto i_element = neighbour_elements.begin(); i_element != neighbour_elements.end(); i_element++)
				if (ElementHas(*i_element, *i_coarse_node))
					i_element->Set(TO_ERASE);
				else
					SwapElementNode(*i_element, rThisNode, i_coarse_node->shared_from_this());
		}
	}

	double MeshNodeCollapsingProcess::CalculateQualityIfNodeCollapses(Node& rThisNode, Node const& rCoarseNode) {
		Point original_coordinates = rThisNode;
		rThisNode.Coordinates() = rCoarseNode.Coordinates();
		double min_quality = CalculateMinQualityOfNeighbourElements(rThisNode, rCoarseNode);
		rThisNode.Coordinates() = original_coordinates;
		return min_quality;
	}

	double MeshNodeCollapsingProcess::CalculateMinQualityOfNeighbourElements(Node& rThisNode, Node const& rCoarseNode) {
		auto& neighbour_elements = rThisNode.GetValue(NEIGHBOUR_ELEMENTS);
		double min_quality = std::numeric_limits<double>::max();

		for (auto i_element = neighbour_elements.begin(); i_element != neighbour_elements.end(); i_element++)
			if (!ElementHas(*i_element, rCoarseNode)) {
				double domain_size = i_element->GetGeometry().DomainSize();
				min_quality = std::min(min_quality, domain_size);
			}
		return min_quality;
	}

	bool MeshNodeCollapsingProcess::ElementHas(Element& rElement, Node const& rCoarseNode) {
		for (auto i_node = rElement.GetGeometry().begin(); i_node != rElement.GetGeometry().end(); i_node++)
			if (i_node->Id() == rCoarseNode.Id())
				return true;

		return false;
	}

	void MeshNodeCollapsingProcess::SwapElementNode(Element& rElement, 
		Node const& rThisNode, Node::Pointer pCoarseNode) {
		std::size_t number_of_nodes = rElement.GetGeometry().size();
		auto& geometry = rElement.GetGeometry();
		for (std::size_t i = 0; i < number_of_nodes; i++)
			if (geometry[i].Id() == rThisNode.Id())
				geometry(i) = pCoarseNode;
	}

}  // namespace Kratos.


