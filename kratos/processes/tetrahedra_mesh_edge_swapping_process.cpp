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
#include "includes/define.h"
#include "processes/tetrahedra_mesh_edge_swapping_process.h"
#include "geometries/tetrahedra_3d_4.h"
#include "processes/element_erase_process.h"



namespace Kratos
{

TetrahedraMeshEdgeSwappingProcess::TetrahedraMeshEdgeSwappingProcess(ModelPart & rModelPart): mrModelPart(rModelPart), mEdges(){

}

TetrahedraMeshEdgeSwappingProcess::~TetrahedraMeshEdgeSwappingProcess(){

}

void TetrahedraMeshEdgeSwappingProcess::Execute(){
	std::cout << std::endl;

	constexpr int tetrahedra_edges[6][2] = { { 0,1 },{ 1,2 },{ 2,0 },{ 0,3 },{ 1,3 },{ 2,3 } };
	for(auto i_element = mrModelPart.ElementsBegin() ; i_element != mrModelPart.ElementsEnd() ; i_element++){
		auto& element_geometry = i_element->GetGeometry();
	 	//std::cout << "Processing element #" << i_element->Id() << "[" << element_geometry[0].Id() << "," << element_geometry[1].Id() << ","
	  //<< element_geometry[2].Id() << "," << element_geometry[3].Id()	<< "]" << std::endl;
		for(int i = 0 ; i < 6 ; i++){
				auto i_edge = mEdges.find(Edge(&(element_geometry[tetrahedra_edges[i][0]]), &(element_geometry[tetrahedra_edges[i][1]])));
				if(i_edge == mEdges.end())
					i_edge = mEdges.emplace(std::make_pair(Edge(&(element_geometry[tetrahedra_edges[i][0]]), &(element_geometry[tetrahedra_edges[i][1]])), TetrahedraEdgeShell(element_geometry(tetrahedra_edges[i][0]),element_geometry(tetrahedra_edges[i][1])))).first;

				//std::cout << "Before: edge " << i_edge->first.GetPoint1()->Id() << " -> " << i_edge->first.GetPoint2()->Id() << " has ";
				//i_edge->second.PrintData(std::cout);
				//std::cout << std::endl;
				i_edge->second.AddElement((*i_element.base()).get(),i);
				//std::cout << "After : edge " << i_edge->first.GetPoint1()->Id() << " -> " << i_edge->first.GetPoint2()->Id() << " has ";
				//i_edge->second.PrintData(std::cout);
				//std::cout << std::endl;
				auto i_edge1 = mEdges.find(Edge(&(element_geometry[tetrahedra_edges[i][0]]), &(element_geometry[tetrahedra_edges[i][1]])));
				//std::cout << "After :edge1 " << i_edge1->first.GetPoint1()->Id() << " -> " << i_edge1->first.GetPoint2()->Id() << " has ";
				//i_edge1->second.PrintData(std::cout);
				//std::cout << std::endl;
			}
	}
	for (auto& edge : mEdges) {
		//std::cout << "Before: edge " << edge.first.GetPoint1()->Id() << " -> " << edge.first.GetPoint2()->Id() << " has ";
		//edge.second.PrintData(std::cout);
		//std::cout << std::endl;
		edge.second.AddShellPoints();
		//std::cout << "After : edge " << edge.first.GetPoint1()->Id() << " -> " << edge.first.GetPoint2()->Id() << " has ";
		//edge.second.PrintData(std::cout);
		//std::cout << std::endl;
	}

	std::array<int, 100> edge_counter;
	for (auto& i : edge_counter)
		i = 0;

	for(auto& edge : mEdges){
		auto size = edge.second.GetNumberOfShellPoints();
		auto tet_numbers = edge.second.GetNumberOfTetrahedra();
		if(size < 100)
			edge_counter[size]++;
		//if ((tet_numbers == 1 && size !=2) || (tet_numbers == 2 && size != 3) || (tet_numbers == 3 && size != 3) || (tet_numbers == 4 && size != 4) || (tet_numbers == 5 && size != 5) || (tet_numbers == 6 && size != 6))
		//{
		//	std::cout << "edge " << edge.first.GetPoint1()->Id() << " -> " << edge.first.GetPoint2()->Id() << " has " << edge.second.GetNumberOfTetrahedra() << " tetrahedras and " << size << " points ";
		//	edge.second.PrintData(std::cout);
		//	std::cout << std::endl;
		//}
	}
	for(std::size_t i = 0 ; i < edge_counter.size() ; i++)
		if(edge_counter[i] > 0)
			std::cout << edge_counter[i] << " edges with " << i << " points" << std::endl;
	KRATOS_WATCH(mEdges.size());

	for (auto& edge : mEdges) {
		if (edge.second.IsClosed()) {
			if (edge.second.GetNumberOfShellPoints() == 3)
				EdgeSwapping3(edge.second);
			if (edge.second.GetNumberOfShellPoints() == 4)
				EdgeSwapping4(edge.second);
		}

	}
	ElementEraseProcess(mrModelPart).Execute();
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
namespace Internals {
	class EdgeSwappingCase {
	public:
		using TetrahedronConnectivityType = std::array<std::size_t, 4>;
		using ConnectivityType = std::vector<TetrahedronConnectivityType>;
		EdgeSwappingCase() = delete;
		EdgeSwappingCase(ConnectivityType const&& TheConnectivities) : mConnectivites(TheConnectivities) {}
		TetrahedronConnectivityType const& GetTetrahedraConnectivites(std::size_t Index) { return mConnectivites[Index]; }
	private:
		ConnectivityType mConnectivites;
	};
	class EdgeSwappingCases {
		int mNumberOfCases;
	protected:
		std::vector<EdgeSwappingCase> mCases;
	public:
		EdgeSwappingCases() = delete;
		EdgeSwappingCases(int NumberOfCases) : mNumberOfCases(NumberOfCases) {}
		int GetNumberOfCases() { return mNumberOfCases; }
		std::vector<EdgeSwappingCase> const& GetCases() { return mCases; }
	};

	class EdgeSwappingCases3 : public EdgeSwappingCases {
	public:
		EdgeSwappingCases3() : EdgeSwappingCases(1) {
			mCases.push_back(EdgeSwappingCase({ { 0, 1, 2 },{ 0, 2, 1 } }));
		}
	};
}

void TetrahedraMeshEdgeSwappingProcess::EdgeSwapping3(TetrahedraEdgeShell & EdgeShell) {
	Internals::EdgeSwappingCases3 SwappingCases;
	for (auto swapping_case : SwappingCases.GetCases()) {
		Tetrahedra3D4<Node<3>> tetrahedra_1(EdgeShell.Point1(), EdgeShell.ShellPoint(swapping_case.GetTetrahedraConnectivites(0)[0]), EdgeShell.ShellPoint(swapping_case.GetTetrahedraConnectivites(0)[1]), EdgeShell.ShellPoint(swapping_case.GetTetrahedraConnectivites(0)[2]));
		Tetrahedra3D4<Node<3>> tetrahedra_2(EdgeShell.Point2(), EdgeShell.ShellPoint(swapping_case.GetTetrahedraConnectivites(1)[0]), EdgeShell.ShellPoint(swapping_case.GetTetrahedraConnectivites(1)[1]), EdgeShell.ShellPoint(swapping_case.GetTetrahedraConnectivites(1)[2]));
		auto quality_criteria = Geometry<Node<3> >::QualityCriteria::VOLUME_TO_AVERAGE_EDGE_LENGTH;

		double original_min_quality = EdgeShell.CalculateMinQuality(quality_criteria);
		double min_quality = std::min(tetrahedra_1.Quality(quality_criteria), tetrahedra_2.Quality(quality_criteria));
		if (min_quality > original_min_quality) {
			EdgeShell.pGetElement(0)->GetGeometry() = tetrahedra_1;
			EdgeShell.pGetElement(1)->GetGeometry() = tetrahedra_2;
			EdgeShell.pGetElement(2)->Set(TO_ERASE);
		}
	}
	//else
	//	std::cout << min_quality << " is worst respect to " << original_min_quality << std::endl;

}
void TetrahedraMeshEdgeSwappingProcess::EdgeSwapping4(TetrahedraEdgeShell & EdgeShell) {

}

}  // namespace Kratos.
