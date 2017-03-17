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
		double low[3];
		double high[3];

		for (int i = 0; i < 3; i++)
		{
			low[i] = high[i] = mrModelPart1.NodesBegin()->Coordinate(i + 1);
		}

		// loop over all nodes in the bounding box
		for (ModelPart::NodeIterator i_node = mrModelPart1.NodesBegin();
		i_node != mrModelPart1.NodesEnd();
			i_node++)
		{
			for (int i = 0; i < 3; i++)
			{
				low[i] = i_node->Coordinate(i + 1) < low[i] ? i_node->Coordinate(i + 1) : low[i];
				high[i] = i_node->Coordinate(i + 1) > high[i] ? i_node->Coordinate(i + 1) : high[i];
			}
		}

		// loop over all skin nodes
		for (ModelPart::NodeIterator i_node = mrModelPart2.NodesBegin();
		i_node != mrModelPart2.NodesEnd();
			i_node++)
		{
			for (int i = 0; i < 3; i++)
			{
				low[i] = i_node->Coordinate(i + 1) < low[i] ? i_node->Coordinate(i + 1) : low[i];
				high[i] = i_node->Coordinate(i + 1) > high[i] ? i_node->Coordinate(i + 1) : high[i];
			}
		}

		mOctree.SetBoundingBox(low, high);

		// Adding mrModelPart2 to the octree
		for (ModelPart::NodeIterator i_node = mrModelPart2.NodesBegin();
		i_node != mrModelPart2.NodesEnd();
			i_node++)
		{
			double temp_point[3];
			temp_point[0] = i_node->X();
			temp_point[1] = i_node->Y();
			temp_point[2] = i_node->Z();
			mOctree.Insert(temp_point);
		}
		
		for (ModelPart::ElementIterator i_element = mrModelPart2.ElementsBegin();
		i_element != mrModelPart2.ElementsEnd();
			i_element++)
		{
			mOctree.Insert(*(i_element).base());
		}
	}


}  // namespace Kratos.


