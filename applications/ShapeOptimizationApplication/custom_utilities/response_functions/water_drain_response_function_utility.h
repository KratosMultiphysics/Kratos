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

/// Short class definition.
/**
 Helper structure storing the information of the pond, such as lowest node id and coordinates, maximal water level,
 neighbouring ponds and if the pond has been merged into another pond.
 */
struct Volume {
	IndexType Id;
	IndexType mLowPointId;
	Vector mLowestPoint;
	double mMaxWaterLevel = 0;
	bool isMerged = false;
	std::vector<int> mNeighbourVolumes;
};

/// Short class definition.
/** Water Drain response function to prevent ponding.
    It computes the water volume by integrating the water level of the ponds.
    For this, a node search is conducted which detects the nodes which are forming a pond.
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

	/// @brief Determines the water volumes or ponds on the surface and the water level of these.
	void SearchWaterVolumes();

	/// @brief Finds lowest points of the ponds (elliptic points of the surface).
	void SearchLowPoints();

	/// @brief Determines whereto water flows from each node and adds nodes to their respective ponds.
	void FindSteepestDescentFromEachNode();

	/// @brief Flattens the ponds so that their maximum water level is correct.
	void LevelVolumes();

	/// @brief Merges neighbouring ponds if necessary.
	void MergeVolumes();

	///@}
	///@name Static Member Variables
	///@{

	///@}
	///@name Member Variables
	///@{

	ModelPart &mrModelPart;
	GeometryUtilities mGeometryUtilities;

	Vector mGravityDirection;
	int mMaxIterations;
	bool mContinuousSens;
	bool mQuadraticHeightPenalization;
	bool mDetectFreeEdgeAutomatic = false;
	std::string mFreeEdgeSubModelPartName;

	std::vector<Volume> mListOfVolumes;
	double mValue;
	///@}

}; // Class WaterDrainResponseFunctionUtility

///@}

} // namespace Kratos.

#endif // WATER_DRAIN_FUNCTION_UTILITY_H
