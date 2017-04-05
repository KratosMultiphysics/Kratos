//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Davand
//
	           
// System includes


// External includes 


// Project includes
#include "processes/find_intersected_geometrical_objects_process.h"


namespace Kratos
{

	FindIntersectedGeometricalObjectsProcess::FindIntersectedGeometricalObjectsProcess(ModelPart& rPart1, ModelPart& rPart2) 
	:mrModelPart1(rPart1), mrModelPart2(rPart2) {

	}

	void FindIntersectedGeometricalObjectsProcess::Initialize() {

	}

	void FindIntersectedGeometricalObjectsProcess::Execute() {
		GenerateOctree();

		std::vector<OctreeType::cell_type*> leaves;
		for (auto p_element_1 : mrModelPart1.ElementsArray()) {
			leaves.clear();
			mOctree.GetIntersectedLeaves(p_element_1, leaves);
			MarkIfIntersected(*p_element_1, leaves);
		}
	}


	/// Turn back information as a string.
	std::string FindIntersectedGeometricalObjectsProcess::Info() const {
		return "FindIntersectedGeometricalObjectsProcess";
	}

	/// Print information about this object.
	void FindIntersectedGeometricalObjectsProcess::PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
	}

	/// Print object's data.
	void FindIntersectedGeometricalObjectsProcess::PrintData(std::ostream& rOStream) const {

	}

	void FindIntersectedGeometricalObjectsProcess::GenerateOctree() {
		SetOctreeBoundingBox();

		// Adding mrModelPart2 to the octree
		for (auto i_node = mrModelPart2.NodesBegin(); i_node != mrModelPart2.NodesEnd(); i_node++) {
			double temp_point[3];
			temp_point[0] = i_node->X();
			temp_point[1] = i_node->Y();
			temp_point[2] = i_node->Z();
			mOctree.Insert(temp_point);
		}
		
		for (auto i_element = mrModelPart2.ElementsBegin(); i_element != mrModelPart2.ElementsEnd(); i_element++) {
			mOctree.Insert(*(i_element).base());
		}
	}

	void  FindIntersectedGeometricalObjectsProcess::SetOctreeBoundingBox() {
		Point<3> low(mrModelPart1.NodesBegin()->Coordinates());
		Point<3> high(mrModelPart1.NodesBegin()->Coordinates());

		// loop over all nodes in first modelpart
		for (auto i_node = mrModelPart1.NodesBegin(); i_node != mrModelPart1.NodesEnd(); i_node++) {
			for (int i = 0; i < 3; i++) {
				low[i] = i_node->Coordinate(i + 1) < low[i] ? i_node->Coordinate(i + 1) : low[i];
				high[i] = i_node->Coordinate(i + 1) > high[i] ? i_node->Coordinate(i + 1) : high[i];
			}
		}

		// loop over all skin nodes
		for (auto i_node = mrModelPart2.NodesBegin(); i_node != mrModelPart2.NodesEnd(); i_node++) {
			for (int i = 0; i < 3; i++) {
				low[i] = i_node->Coordinate(i + 1) < low[i] ? i_node->Coordinate(i + 1) : low[i];
				high[i] = i_node->Coordinate(i + 1) > high[i] ? i_node->Coordinate(i + 1) : high[i];
			}
		}

		// TODO: Octree needs refactoring to work with BoundingBox. Pooyan.
		mOctree.SetBoundingBox(low.data().data(), high.data().data());
	}

	void  FindIntersectedGeometricalObjectsProcess::MarkIfIntersected(Element& rElement1, std::vector<OctreeType::cell_type*>& leaves) {
		for (auto p_leaf : leaves) {
			for (auto p_element_2 : *(p_leaf->pGetObjects())) {
				if (HasIntersection(rElement1.GetGeometry(),p_element_2->GetGeometry())) {
					rElement1.Set(SELECTED);
					return;
				}
			}
		}
	}

	bool FindIntersectedGeometricalObjectsProcess::HasIntersection(Element::GeometryType& rFirstGeometry, Element::GeometryType& rSecondGeometry) {
		auto faces = rFirstGeometry.Faces();
		for (auto& face : faces) {
			if (face.HasIntersection(rSecondGeometry))
				return true;
		}
		// Let check second geometry is inside the first one.
		// Considering that there are no intersection, if one point is inside all of it is inside.
		array_1d<double, 3> local_point;
		if (rFirstGeometry.IsInside(rSecondGeometry.GetPoint(0), local_point))
			return true;
		return false;
	}


}  // namespace Kratos.


