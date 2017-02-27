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
	TetrahedraEdgeShell::TetrahedraEdgeShell(PointType::Pointer pEdgePoint1, PointType::Pointer pEdgePoint2)
	: mpPoint1(pEdgePoint1) , mpPoint2(pEdgePoint2), mShellPoints(), mTetrahedraEdge(), mIsClosed(false) {
	}

	TetrahedraEdgeShell::~TetrahedraEdgeShell() {

	}

	void TetrahedraEdgeShell::AddElement(Element* pTheElement, char EdgeIndex) {
		mTetrahedraEdge.push_back(std::make_pair(pTheElement, EdgeIndex));
	}

	Element* TetrahedraEdgeShell::pGetElement(std::size_t ElementLocalIndex) {
		return mTetrahedraEdge[ElementLocalIndex].first;
	}

	void TetrahedraEdgeShell::AddShellPoints() {
		if (mTetrahedraEdge.empty())
			return;
		//std::cout << "Processing edge " << mpPoint1.Id() << " -> " << mpPoint2.Id() << ": ";
		//PrintData(std::cout);
		//std::cout << std::endl;
		constexpr int tetrahedra_edges[6][2] = { { 0,1 },{ 1,2 },{ 2,0 },{ 0,3 },{ 1,3 },{ 2,3 } };
		constexpr int tetrahedra_edge_other_points[6][2] = { {2,3}, {0,3}, {1,3}, {1,2}, {2,0}, {0,1} };
		auto& tetrahedra_edge = mTetrahedraEdge.begin();
		auto& tetrahedra = (mTetrahedraEdge[0]).first->GetGeometry();
		auto edge_index = (mTetrahedraEdge[0]).second;
		//std::cout << "    First two points are ";
		if (tetrahedra(tetrahedra_edges[edge_index][0]) == mpPoint1) {
			mShellPoints.push_back(tetrahedra(tetrahedra_edge_other_points[edge_index][0]));
			mShellPoints.push_back(tetrahedra(tetrahedra_edge_other_points[edge_index][1]));
			//std::cout << tetrahedra[tetrahedra_edge_other_points[edge_index][0]].Id() << " and " << tetrahedra[tetrahedra_edge_other_points[edge_index][1]].Id() << std::endl;
		}
		else
		{
			mShellPoints.push_back(tetrahedra(tetrahedra_edge_other_points[edge_index][1]));
			mShellPoints.push_back(tetrahedra(tetrahedra_edge_other_points[edge_index][0]));
			//std::cout << tetrahedra[tetrahedra_edge_other_points[edge_index][1]].Id() << " and " << tetrahedra[tetrahedra_edge_other_points[edge_index][0]].Id() << std::endl;

		}
		std::vector<bool> is_inserted(mTetrahedraEdge.size(), false);
		is_inserted[0] = true;
		for (std::size_t i = 1; i < mTetrahedraEdge.size(); i++) {
			//std::cout << "    i : " << i << std::endl;
			for (std::size_t j = 1; j < mTetrahedraEdge.size(); j++) {
				//std::cout << "        j : " << j ;
				if (is_inserted[j] == false) {
					auto& tetrahedra = (mTetrahedraEdge[j].first->GetGeometry());
					//std::cout << " [ ";
					//for (int k = 0; k < 4; k++)
					//	std::cout << tetrahedra[k].Id() << " ";


					edge_index = mTetrahedraEdge[j].second;
					if (tetrahedra(tetrahedra_edges[edge_index][0]) == mpPoint1) {
						if (tetrahedra(tetrahedra_edge_other_points[edge_index][0]) == mShellPoints.back()) {
							mShellPoints.push_back(tetrahedra(tetrahedra_edge_other_points[edge_index][1]));
							is_inserted[j] = true;
						}
						else if (tetrahedra(tetrahedra_edge_other_points[edge_index][1]) == mShellPoints.front()) {
							mShellPoints.insert(mShellPoints.begin(), tetrahedra(tetrahedra_edge_other_points[edge_index][0]));
							is_inserted[j] = true;
						}
					}
					else
					{
						if (tetrahedra(tetrahedra_edge_other_points[edge_index][1]) == mShellPoints.back()){
							mShellPoints.push_back(tetrahedra(tetrahedra_edge_other_points[edge_index][0]));
							is_inserted[j] = true;
						}
						else if (tetrahedra(tetrahedra_edge_other_points[edge_index][0]) == mShellPoints.front()) {
							mShellPoints.insert(mShellPoints.begin(), tetrahedra(tetrahedra_edge_other_points[edge_index][1]));
							is_inserted[j] = true;
						}
					}
					//std::cout << " ] and points : [  ";
					//for (auto i_point = mShellPoints.begin(); i_point != mShellPoints.end(); i_point++)
					//	std::cout << (*i_point)->Id() << " ";
					//std::cout << "]";
				}
			}
			if (mShellPoints.front() == mShellPoints.back()) {
				mIsClosed = true;
				mShellPoints.pop_back();
				break;
			}
		}
		//std::cout << std::endl;
	}

	double TetrahedraEdgeShell::CalculateMinQuality(const Geometry<Node<3>>::QualityCriteria QualityCriteria) const {
		if (mTetrahedraEdge.size() == 0)
			return 0.00;

		double min_quality = std::numeric_limits<double>::max();
		for (auto i_tetraheron = mTetrahedraEdge.begin(); i_tetraheron != mTetrahedraEdge.end(); i_tetraheron++) {
			min_quality = std::min(min_quality, (i_tetraheron->first)->GetGeometry().Quality(QualityCriteria));
			KRATOS_WATCH((i_tetraheron->first)->GetGeometry())
				KRATOS_WATCH((i_tetraheron->first)->GetGeometry().Quality(QualityCriteria))
				KRATOS_WATCH(min_quality)
		}

		return min_quality;
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
		rOStream << "tetrahedra : [  ";
		for (auto i_tetrahedra = mTetrahedraEdge.begin(); i_tetrahedra != mTetrahedraEdge.end(); i_tetrahedra++) {
			rOStream << " [ ";
				for (int i = 0; i < 4; i++)
					rOStream << (i_tetrahedra->first->GetGeometry())[i].Id() << " ";
				rOStream << "]";

		}
		rOStream << " ] and points : [  ";
		for (auto i_point = mShellPoints.begin(); i_point != mShellPoints.end(); i_point++)
			rOStream << (*i_point)->Id() << " ";
		rOStream << "]";
	}


}  // namespace Kratos.
