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
		EdgeSwappingCase() {};
		EdgeSwappingCase(std::vector<std::size_t> const&& TheTrianglesIndices) : mTriangleIndices(TheTrianglesIndices), mMinQuality(2.00){}
		std::size_t GetTringleIndex(std::size_t Index) const { return mTriangleIndices[Index]; }
		double GetMinQuality() const { return mMinQuality; }
		void SetMinQuality(double NewMinQuality) { mMinQuality = NewMinQuality; }
	private:
		std::vector<std::size_t> mTriangleIndices;
		double mMinQuality;
	};
	template<std::size_t TNumberOfCases, std::size_t TNumberOfTriangles, std::size_t TNumberOfTrianglesPerCase>
	class EdgeSwappingCases {
	public:
		using TriangleConnectivityType = std::array<std::size_t, 3>;
		static std::size_t GetNumberOfCases() { return TNumberOfCases; }
		static std::size_t NumberOfTriangles() { return TNumberOfTriangles; }
		static std::size_t NumberOfTrianglesPerCase() { return TNumberOfTrianglesPerCase; }
		std::array<EdgeSwappingCase, TNumberOfCases> const& GetCases() { return mCases; }
		TriangleConnectivityType const& GetTriangleConectivity(std::size_t TheIndex) { return mTriangles[TheIndex]; }
		void SetTetrahedraForCase(EdgeSwappingCase const& TheCase,std::size_t TriangleIndex, TetrahedraEdgeShell& EdgeShell, Tetrahedra3D4<Node<3>>& rTetrahedra1, Tetrahedra3D4<Node<3>>& rTetrahedra2) {
			auto const& triangle = GetTriangleConectivity(TheCase.GetTringleIndex(TriangleIndex));
			rTetrahedra1(0) = EdgeShell.Point1();
			rTetrahedra1(1) = EdgeShell.ShellPoint(triangle[0]);
			rTetrahedra1(2) = EdgeShell.ShellPoint(triangle[1]);
			rTetrahedra1(3) = EdgeShell.ShellPoint(triangle[2]);

			rTetrahedra2(0) = EdgeShell.Point2();
			rTetrahedra2(1) = EdgeShell.ShellPoint(triangle[0]);
			rTetrahedra2(2) = EdgeShell.ShellPoint(triangle[2]);
			rTetrahedra2(3) = EdgeShell.ShellPoint(triangle[1]);
		}
	protected:
		EdgeSwappingCases() {}
		std::array<TriangleConnectivityType, TNumberOfTriangles> mTriangles;
		std::array<EdgeSwappingCase, TNumberOfCases> mCases;
	};

	class EdgeSwappingCases3 : public EdgeSwappingCases<1, 1, 1> {
	public:
		EdgeSwappingCases3() : EdgeSwappingCases() {
			mCases[0] = EdgeSwappingCase({ 0 });
			mTriangles[0] = { 0,1,2 };
		}
	};

	class EdgeSwappingCases4 : public EdgeSwappingCases<2, 4, 2> {
	public:
		EdgeSwappingCases4() : EdgeSwappingCases() {
			mCases[0] = EdgeSwappingCase({ 0,1 });
			mCases[1] = EdgeSwappingCase({ 2,3 });
			mTriangles[0] = { 0,1,2 };
			mTriangles[1] = { 0,2,3 };
			mTriangles[2] = { 1,2,3 };
			mTriangles[3] = { 0,1,3 };
		}
	};
}

void TetrahedraMeshEdgeSwappingProcess::EdgeSwapping3(TetrahedraEdgeShell & EdgeShell) {
	Internals::EdgeSwappingCases3 SwappingCases;
	auto swapping_case = SwappingCases.GetCases()[0];
	auto const& triangle = SwappingCases.GetTriangleConectivity(swapping_case.GetTringleIndex(0));
	Tetrahedra3D4<Node<3>> tetrahedra_1(EdgeShell.Point1(), EdgeShell.ShellPoint(triangle[0]), EdgeShell.ShellPoint(triangle[1]), EdgeShell.ShellPoint(triangle[2]));
	Tetrahedra3D4<Node<3>> tetrahedra_2(EdgeShell.Point2(), EdgeShell.ShellPoint(triangle[0]), EdgeShell.ShellPoint(triangle[2]), EdgeShell.ShellPoint(triangle[1]));
	auto quality_criteria = Geometry<Node<3> >::QualityCriteria::SHORTEST_TO_LONGEST_EDGE;

	double original_min_quality = EdgeShell.CalculateMinQuality(quality_criteria);
	double min_quality = std::min(tetrahedra_1.Quality(quality_criteria), tetrahedra_2.Quality(quality_criteria));
	if (min_quality > original_min_quality) {
		EdgeShell.pGetElement(0)->GetGeometry() = tetrahedra_1;
		EdgeShell.pGetElement(1)->GetGeometry() = tetrahedra_2;
		EdgeShell.pGetElement(2)->Set(TO_ERASE);
	}
	//else
	//	std::cout << min_quality << " is worst respect to " << original_min_quality << std::endl;
}

void TetrahedraMeshEdgeSwappingProcess::EdgeSwapping4(TetrahedraEdgeShell & EdgeShell) {
	Internals::EdgeSwappingCases4 SwappingCases;
	auto quality_criteria = Geometry<Node<3> >::QualityCriteria::SHORTEST_TO_LONGEST_EDGE;
	double original_min_quality = EdgeShell.CalculateMinQuality(quality_criteria);
	KRATOS_WATCH(original_min_quality)
	Tetrahedra3D4<Node<3>> tetrahedra_1 = EdgeShell.pGetElement(0)->GetGeometry(); // This initialization is to avoid creating a dummy
	Tetrahedra3D4<Node<3>> tetrahedra_2 = EdgeShell.pGetElement(0)->GetGeometry(); // It will be reinitialized afterward
	std::size_t best_case = 0; 
	double max_cases_quality = original_min_quality; 

	for (auto i_case = SwappingCases.GetCases().begin(); i_case != SwappingCases.GetCases().end(); i_case++) {
		if (i_case->GetMinQuality() > original_min_quality) {// There are no previously calculated tetrahedra with worse quality
			for (std::size_t i = 0; i < SwappingCases.NumberOfTrianglesPerCase(); i++) {
				SwappingCases.SetTetrahedraForCase(*i_case,i, EdgeShell, tetrahedra_1, tetrahedra_2);
				double min_quality = std::min(tetrahedra_1.Quality(quality_criteria), tetrahedra_2.Quality(quality_criteria));
				if (min_quality > max_cases_quality) {
					best_case = i;
					max_cases_quality = min_quality;
					// Todo: break if apt quality reached.
				}
			}
		}
	}
	KRATOS_WATCH(max_cases_quality)
	if (max_cases_quality > original_min_quality + std::numeric_limits<double>::epsilon()) {
		Tetrahedra3D4<Node<3>> tetrahedra_3 = EdgeShell.pGetElement(0)->GetGeometry(); // This initialization is to avoid creating a dummy
		Tetrahedra3D4<Node<3>> tetrahedra_4 = EdgeShell.pGetElement(0)->GetGeometry(); // It will be reinitialized afterward
		SwappingCases.SetTetrahedraForCase(SwappingCases.GetCases()[best_case], 0, EdgeShell, tetrahedra_1, tetrahedra_2);
		SwappingCases.SetTetrahedraForCase(SwappingCases.GetCases()[best_case], 1, EdgeShell, tetrahedra_3, tetrahedra_4);
		EdgeShell.pGetElement(0)->GetGeometry() = tetrahedra_1;
		EdgeShell.pGetElement(1)->GetGeometry() = tetrahedra_2;
		EdgeShell.pGetElement(2)->GetGeometry() = tetrahedra_3;
		EdgeShell.pGetElement(3)->GetGeometry() = tetrahedra_4;
		//for (std::size_t i = 0; i < SwappingCases.NumberOfTrianglesPerCase(); i++) {
		//	SwappingCases.SetTetrahedraForCase(i, EdgeShell, tetrahedra_1, tetrahedra_2);
		//	EdgeShell.pGetElement(2*i)->GetGeometry() = tetrahedra_1;
		//	EdgeShell.pGetElement((2*i)+1)->GetGeometry() = tetrahedra_2;
		//}
	}
}


}  // namespace Kratos.
