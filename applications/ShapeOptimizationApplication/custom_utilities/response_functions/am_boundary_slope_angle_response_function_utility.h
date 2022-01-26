// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//
// ==============================================================================

#ifndef AM_BOUNDARY_SLOPE_ANGLE_RESPONSE_FUNCTION_UTILITY_H
#define AM_BOUNDARY_SLOPE_ANGLE_RESPONSE_FUNCTION_UTILITY_H

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

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) AMBoundarySlopeAngleResponseFunctionUtility
{
public:
	///@name Type Definitions
	///@{

	typedef array_1d<double, 3> array_3d;

	/// Pointer definition of AMBoundarySlopeAngleResponseFunctionUtility
	KRATOS_CLASS_POINTER_DEFINITION(AMBoundarySlopeAngleResponseFunctionUtility);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	AMBoundarySlopeAngleResponseFunctionUtility(ModelPart& rModelPart, Parameters ResponseSettings);

	/// Destructor.
	virtual ~AMBoundarySlopeAngleResponseFunctionUtility()
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
		return "AMBoundarySlopeAngleResponseFunctionUtility";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "AMBoundarySlopeAngleResponseFunctionUtility";
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
	double mCosMinAngle;
	double mValue;

	///@}

}; // Class AMBoundarySlopeAngleResponseFunctionUtility

///@}

} // namespace Kratos.

#endif // AM_BOUNDARY_SLOPE_ANGLE_RESPONSE_FUNCTION_UTILITY_H
