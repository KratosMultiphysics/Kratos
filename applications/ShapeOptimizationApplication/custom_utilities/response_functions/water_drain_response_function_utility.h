// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Schmölz David, https://github.com/dschmoelz
//
// ==============================================================================

#ifndef WATER_DRAIN_FUNCTION_UTILITY_H
#define WATER_DRAIN_FUNCTION_UTILITY_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "custom_utilities/geometry_utilities.h"


// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

typedef Node NodeType;
typedef NodeType::Pointer NodeTypePointer;
typedef std::vector<NodeType::Pointer> NodeVector;
typedef ModelPart::NodesContainerType NodesArrayType;

struct Volume {
	double mValue = 0.0;
	NodesArrayType mListOfNodes;
	NodesArrayType mNeighbourNodes;
	std::vector<std::pair<double, NodeTypePointer>> mNeighbourNodesSorted;
	Vector mHighestPoint;
	Vector mLowestPoint;
	bool isGrowing = true;
	bool isMerged = false;
};

/// Short class definition.
/** Detail class definition.

 */

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) WaterDrainResponseFunctionUtility
{
public:
	///@name Type Definitions
	///@{

	typedef array_1d<double, 3> array_3d;

	/// Pointer definition of WaterDrainResponseFunctionUtility
	KRATOS_CLASS_POINTER_DEFINITION(WaterDrainResponseFunctionUtility);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	WaterDrainResponseFunctionUtility(ModelPart& rModelPart, Parameters ResponseSettings);

	/// Destructor.
	virtual ~WaterDrainResponseFunctionUtility()
	{}
	///@}
	///@name Operators
	///@{

	///@}
	///@name Operations
	///@{

	void Initialize();

	void InitializeSolutionStep();

	double CalculateValue();

	void CalculateGradient();

	// ==============================================================================

	///@name Input and output
	///@{

	/// Turn back information as a string.
	std::string Info() const
	{
		return "WaterDrainResponseFunctionUtility";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "WaterDrainResponseFunctionUtility";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
	{
	}

	///@}

protected:

	// ==============================================================================

private:

	///@name Operations
	///@{

	void SearchWaterVolumes();

	void GrowVolume(Volume& rVolume);

	void MergeVolumes();

	void MergeTwoVolumes(Volume& rVolume1, Volume& rVolume2);

	void SearchLowPoints();

	double CalculateConditionValue(const Condition& rFace);

	///@}
	///@name Static Member Variables
	///@{

	///@}
	///@name Member Variables
	///@{

	ModelPart &mrModelPart;
	GeometryUtilities mGeometryUtilities;
	double mDelta;
	Vector mGravityDirection;
	int mMaxIterations;
	double mValue;
	std::vector<Volume> mListOfVolumes;
	bool mContinuousSens;
	bool mQuadraticHeightPenalization;
	std::string mEdgeSubModelPartName;

	///@}

}; // Class WaterDrainResponseFunctionUtility

///@}

} // namespace Kratos.

#endif // WATER_DRAIN_FUNCTION_UTILITY_H
