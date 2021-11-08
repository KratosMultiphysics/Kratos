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
#include "modeler/tetrahedra_ball.h"
#include "includes/element.h"


namespace Kratos
{

	TetrahedraBall::TetrahedraBall(NodeType& rThisNode) {
		auto& neighbour_elements = rThisNode.GetValue(NEIGHBOUR_ELEMENTS);
		for (auto i_element = neighbour_elements.begin(); i_element != neighbour_elements.end(); i_element++) {
			mTetrahedra.push_back(i_element->pGetGeometry().get());
		}

	}

	std::size_t TetrahedraBall::Size() const {
		return mTetrahedra.size();
	}

	double TetrahedraBall::CalculateMinQuality(const Geometry<NodeType>::QualityCriteria QualityCriteria) const {
		if (Size() == 0)
			return 0.00;

		double min_quality = std::numeric_limits<double>::max();
		for (auto i_tetraheron = mTetrahedra.begin(); i_tetraheron != mTetrahedra.end(); i_tetraheron++) {
			min_quality = std::min(min_quality, (*i_tetraheron)->Quality(QualityCriteria));
			}

		return min_quality;
	}

	std::string TetrahedraBall::Info() const {
		return "TetrahedraBall";
	}

	void TetrahedraBall::PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
	}

	void TetrahedraBall::PrintData(std::ostream& rOStream) const {
		rOStream << "Size : " << Size();
	}

}  // namespace Kratos.


