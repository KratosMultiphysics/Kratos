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
#include "processes/calculate_distance_to_skin_process.h"


namespace Kratos
{

	CalculateDistanceToSkinProcess::CalculateDistanceToSkinProcess(ModelPart& rVolumePart, ModelPart& rSkinPart)
		: CalculateDiscontinuousDistanceToSkinProcess(rVolumePart, rSkinPart) {

	}

	CalculateDistanceToSkinProcess::~CalculateDistanceToSkinProcess() {}


	void CalculateDistanceToSkinProcess::Execute() {
		for (auto& node : GetModelPart1().Nodes()) {
			node.GetSolutionStepValue(DISTANCE) =  std::numeric_limits<double>::max();
		}

		CalculateDiscontinuousDistanceToSkinProcess::Execute();

	}

	/// Turn back information as a string.
	std::string CalculateDistanceToSkinProcess::Info() const {
		return "CalculateDistanceToSkinProcess";
	}

	/// Print information about this object.
	void CalculateDistanceToSkinProcess::PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
	}

	/// Print object's data.
	void CalculateDistanceToSkinProcess::PrintData(std::ostream& rOStream) const {
	}

  
  
}  // namespace Kratos.


