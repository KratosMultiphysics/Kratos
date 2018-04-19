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
//                   Riccardo Rossi
//                   
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "processes/tetrahedra_mesh_edge_swapping_process.h"



namespace Kratos
{

TetrahedraMeshEdgeSwappingProcess::TetrahedraMeshEdgeSwappingProcess(ModelPart & rModelPart): mrModelPart(rModelPart){

}

TetrahedraMeshEdgeSwappingProcess::~TetrahedraMeshEdgeSwappingProcess(){

}

void TetrahedraMeshEdgeSwappingProcess::Execute(){
	std::cout << std::endl;

	EdgesContainerType edges;
	for(auto i_element = mrModelPart.ElementsBegin() ; i_element != mrModelPart.ElementsEnd() ; i_element++){
		auto& element_geometry = i_element->GetGeometry();
	// 	std::cout << "Processing element #" << i_element->Id() << "[" << element_geometry[0].Id() << "," << element_geometry[1].Id() << ","
	//  << element_geometry[2].Id() << "," << element_geometry[3].Id()	<< "]" << std::endl;
		for(int i = 0 ; i < 4 ; i++)
			for(int j = i+1; j < 4; j++ ){
				auto i_edge = edges.find(Edge(&(element_geometry[i]), &(element_geometry[j])));
				if(i_edge == edges.end())
					i_edge = edges.emplace(std::make_pair(Edge(&(element_geometry[i]), &(element_geometry[j])), TetrahedraEdgeShell(element_geometry[i],element_geometry[j]))).first;

				// std::cout << "Before: edge " << i_edge->first.GetPoint1()->Id() << " -> " << i_edge->first.GetPoint2()->Id() << " has " << i_edge->second.GetNumberOfShellPoints() << " points" << std::endl;
				i_edge->second.AddTetrahedron(&element_geometry);
				// std::cout << "After : edge " << i_edge->first.GetPoint1()->Id() << " -> " << i_edge->first.GetPoint2()->Id() << " has " << i_edge->second.GetNumberOfShellPoints() << " points" << std::endl;
			}
	}
	std::array<int, 100> edge_counter;
	for (auto& i : edge_counter)
		i = 0;

	for(auto& edge : edges){
		auto size = edge.second.GetNumberOfTetrahedra();
		if(size < 100)
			edge_counter[size]++;
		// if(size == 0)
		// std::cout << "edge " << edge.first.GetPoint1()->Id() << " -> " << edge.first.GetPoint2()->Id() << " has " << size << " points" << std::endl;
	}
	for(std::size_t i = 0 ; i < edge_counter.size() ; i++)
		if(edge_counter[i] > 0)
			std::cout << edge_counter[i] << " edges with " << i << " tetrahedra" << std::endl;
	KRATOS_WATCH(edges.size());
}

std::string TetrahedraMeshEdgeSwappingProcess::Info() const{
	 return "TetrahedraMeshEdgeSwappingProcess";
 }

/// Print information about this object.
void TetrahedraMeshEdgeSwappingProcess::PrintInfo(std::ostream& rOStream) const {
	rOStream << Info();
}

/// Print object's data.
void TetrahedraMeshEdgeSwappingProcess::PrintData(std::ostream& rOStream) const {

}

}  // namespace Kratos.
