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

	FindIntersectedGeometricalObjectsProcess::FindIntersectedGeometricalObjectsProcess(ModelPart& rPart1, ModelPart const& rPart2) 
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

	}

  
}  // namespace Kratos.


