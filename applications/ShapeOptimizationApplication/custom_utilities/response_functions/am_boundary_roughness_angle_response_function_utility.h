// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//
// ==============================================================================

#ifndef AM_BOUNDARY_ROUGHNESS_ANGLE_RESPONSE_FUNCTION_UTILITY_H
#define AM_BOUNDARY_ROUGHNESS_ANGLE_RESPONSE_FUNCTION_UTILITY_H

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

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) AMBoundaryRoughnessAngleResponseFunctionUtility
{
public:
	///@name Type Definitions
	///@{

	typedef array_1d<double, 3> array_3d;

	/// Pointer definition of AMBoundaryRoughnessAngleResponseFunctionUtility
	KRATOS_CLASS_POINTER_DEFINITION(AMBoundaryRoughnessAngleResponseFunctionUtility);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	AMBoundaryRoughnessAngleResponseFunctionUtility(ModelPart& rModelPart, Parameters ResponseSettings);

	/// Destructor.
	virtual ~AMBoundaryRoughnessAngleResponseFunctionUtility()
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
		return "AMBoundaryRoughnessAngleResponseFunctionUtility";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "AMBoundaryRoughnessAngleResponseFunctionUtility";
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

	double CalculateConditionValue(const Condition& rFace);

	///@}
	///@name Static Member Variables
	///@{

	///@}
	///@name Member Variables
	///@{

	ModelPart &mrModelPart;
	double mDelta;
	array_3d mMainDirection;
	double mCosMaxAngle;
	double mValue;

	///@}

}; // Class AMBoundaryRoughnessAngleResponseFunctionUtility

///@}

} // namespace Kratos.

#endif // AM_BOUNDARY_ROUGHNESS_ANGLE_RESPONSE_FUNCTION_UTILITY_H
