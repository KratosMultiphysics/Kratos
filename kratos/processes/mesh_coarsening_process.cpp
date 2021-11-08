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
#include "modeler/tetrahedra_ball.h"
#include "includes/global_pointer_variables.h"

namespace Kratos
{

	KRATOS_CREATE_LOCAL_FLAG(MeshCoarseningProcess, COARSE_MESH_NODE, 1);

	MeshCoarseningProcess::MeshCoarseningProcess(ModelPart& rModelPart)
	: MeshNodeCollapsingProcess(rModelPart){

	}

	MeshCoarseningProcess:: ~MeshCoarseningProcess() {

	}

	void MeshCoarseningProcess::Execute() {
		SelectCoarseMeshNodes();
		MeshNodeCollapsingProcess::Execute();
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

		for (auto i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++) {
			if (i_node->Is(COARSE_MESH_NODE)) // We want to keep all the boundary nodes 
				i_node->Set(TO_COLLAPSE,false);
			else
				i_node->Set(TO_COLLAPSE, true);
		}
	}

}  // namespace Kratos.


