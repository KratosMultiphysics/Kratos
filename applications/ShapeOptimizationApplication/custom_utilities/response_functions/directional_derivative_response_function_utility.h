// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Schmölz David, https://github.com/dschmoelz
//
// ==============================================================================

#ifndef DIRECTIONAL_DERIVATIVE_RESPONSE_FUNCTION_UTILITY_H
#define DIRECTIONAL_DERIVATIVE_RESPONSE_FUNCTION_UTILITY_H

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

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) DirectionalDerivativeResponseFunctionUtility
{
public:
	///@name Type Definitions
	///@{

	typedef array_1d<double, 3> array_3d;

    typedef Node < 3 > NodeType;
    typedef std::vector<std::vector<NodeType::Pointer>> ConditionConnectionVector;

	/// Pointer definition of DirectionalDerivativeResponseFunctionUtility
	KRATOS_CLASS_POINTER_DEFINITION(DirectionalDerivativeResponseFunctionUtility);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	DirectionalDerivativeResponseFunctionUtility(ModelPart& rModelPart, Parameters ResponseSettings);

	/// Destructor.
	virtual ~DirectionalDerivativeResponseFunctionUtility()
	{}
	///@}
	///@name Operators
	///@{

	///@}
	///@name Operations
	///@{

	void Initialize();

	double CalculateValue();

	void CalculateGradient();

	// ==============================================================================

	///@name Input and output
	///@{

	/// Turn back information as a string.
	std::string Info() const
	{
		return "DirectionalDerivativeResponseFunctionUtility";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "DirectionalDerivativeResponseFunctionUtility";
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

	double CalculateConditionValue(Condition& rFace);

	void FindAdjacentConditions();

	///@}
	///@name Static Member Variables
	///@{

	///@}
	///@name Member Variables
	///@{

	ModelPart &mrModelPart;
	array_3d mMainDirection;
	double mValue;
	double mTanAngle;
	double mMaxSearchAngle;

	///@}

}; // Class DirectionalDerivativeResponseFunctionUtility

///@}

} // namespace Kratos.

#endif // DIRECTIONAL_DERIVATIVE_RESPONSE_FUNCTION_UTILITY_H
