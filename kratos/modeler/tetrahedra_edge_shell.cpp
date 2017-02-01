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

// System includes

// External includes

// Project includes
#include "modeler/tetrahedra_edge_shell.h"
#include "includes/element.h"
 
namespace Kratos
{
	TetrahedraEdgeShell::TetrahedraEdgeShell(PointType& EdgePoint1, PointType& EdgePoint2)
	: mrPoint1(EdgePoint1) , mrPoint2(EdgePoint2){
		//auto ball1 = mrPoint1.GetValue(NEIGHBOUR_ELEMENTS);
		//for (auto i_element = ball1.begin(); i_element != ball1.end(); i_element++) {
		//	for (auto i_node = i_element->GetGeometry().begin(); i_node = i_element->GetGeometry().end(); i_node++)
		//		if(i_node->Id() == mrPoint2.Id()) // Element belongs to the shell

		//}
	}

	TetrahedraEdgeShell::~TetrahedraEdgeShell() {

	}

	void TetrahedraEdgeShell::AddTetrahedron(GeomertyType* pTheTetrahedron){
		mTetrahedra.push_back(pTheTetrahedron);
		//constexpr int number_of_tetrahedron_points = 4;
		//constexpr int tetrahedra_connectivity[number_of_tetrahedron_points][3] = { {3,2,1},{2,3,0},{0,3,1},{0,1,2} };
		//for(int i = 0 ; i < number_of_tetrahedron_points ; i++){
		//	if(&(TheTetrahedron[i]) == &mrPoint1){
		//		auto face = tetrahedra_connectivity[i];
		//		if (&(TheTetrahedron[face[0]]) == &mrPoint2){
		//			AddShellPoints(&(TheTetrahedron[face[1]]), &(TheTetrahedron[face[2]]));
		//		} else if (&(TheTetrahedron[face[1]]) == &mrPoint2){
		//			AddShellPoints(&(TheTetrahedron[face[2]]), &(TheTetrahedron[face[0]]));
		//		} else {
		//			AddShellPoints(&(TheTetrahedron[face[0]]), &(TheTetrahedron[face[1]]));
		//		}
		//		break;
		//	}
		//}
	}

	void TetrahedraEdgeShell::AddShellPoints(PointType* pPoint1, PointType* pPoint2){
		mShellPoints.push_back(std::make_pair(pPoint1, pPoint2));
	}

	std::string TetrahedraEdgeShell::Info() const {
		return "TetrahedraEdgeShell";
	}

	/// Print information about this object.
	void TetrahedraEdgeShell::PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
	}

	/// Print object's data.
	void TetrahedraEdgeShell::PrintData(std::ostream& rOStream) const {

	}


}  // namespace Kratos.
