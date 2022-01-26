// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//
// ==============================================================================

#ifndef THICKNESS_RESPONSE_FUNCTION_UTILITY_H
#define THICKNESS_RESPONSE_FUNCTION_UTILITY_H

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
#include "spatial_containers/spatial_containers.h"

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

/// Short class definition.
/** Detail class definition.

 */

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) ThicknessResponseFunctionUtility
{
public:
	///@name Type Definitions
	///@{

	typedef array_1d<double, 3> array_3d;
    // Type definitions for tree-search
    typedef ModelPart::NodeType NodeType;
    typedef ModelPart::NodeType::Pointer NodeTypePointer;
    typedef std::vector<NodeTypePointer> NodeVector;
    typedef std::vector<NodeTypePointer>::iterator NodeVectorIterator;
    typedef std::vector<double> DoubleVector;
    typedef std::vector<double>::iterator DoubleVectorIterator;
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeVectorIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;	

	/// Pointer definition of ThicknessResponseFunctionUtility
	KRATOS_CLASS_POINTER_DEFINITION(ThicknessResponseFunctionUtility);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	ThicknessResponseFunctionUtility(ModelPart& rModelPart, Parameters ResponseSettings);

	/// Destructor.
	virtual ~ThicknessResponseFunctionUtility()
	{}
	///@}
	///@name Operators
	///@{

	///@}
	///@name Operations
	///@{

	void Initialize()
	{}

	double CalculateValue();

	void CalculateGradient();

	// ==============================================================================

	///@name Input and output
	///@{

	/// Turn back information as a string.
	std::string Info() const
	{
		return "ThicknessResponseFunctionUtility";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "ThicknessResponseFunctionUtility";
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

	double CalculateConditionValue(const Condition& rFace,array_3d & rQualPointCoords);
	bool FindConditionQaulPoint(const Condition& rFace,array_3d & rQualPointCoords);

	///@}
	///@name Static Member Variables
	///@{

	///@}
	///@name Member Variables
	///@{

	ModelPart &mrModelPart;
	KDTree::Pointer mpSearchTree;
	NodeVector mListOfNodesInModelPart;
	const unsigned int mBucketSize = 100;
    unsigned int mMaxNeighborNodes = 10000;
	double mDelta;
	double mMinThickness;
	double mValue;

	///@}

}; // Class ThicknessResponseFunctionUtility

///@}

} // namespace Kratos.

#endif // THICKNESS_RESPONSE_FUNCTION_UTILITY_H
