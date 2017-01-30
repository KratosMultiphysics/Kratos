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


