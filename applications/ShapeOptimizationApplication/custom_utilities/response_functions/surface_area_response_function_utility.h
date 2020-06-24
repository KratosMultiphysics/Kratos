// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Geiser Armin, https://github.com/armingeiser
//

#ifndef SURFACE_AREA_RESPONSE_FUNCTION_UTILITY_H
#define SURFACE_AREA_RESPONSE_FUNCTION_UTILITY_H

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
#include "utilities/variable_utils.h"

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

class SurfaceAreaResponseFunctionUtility
{
public:
	///@name Type Definitions
	///@{

	typedef array_1d<double, 3> array_3d;

	/// Pointer definition of SurfaceAreaResponseFunctionUtility
	KRATOS_CLASS_POINTER_DEFINITION(SurfaceAreaResponseFunctionUtility);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	SurfaceAreaResponseFunctionUtility(ModelPart& model_part, Parameters responseSettings)
	: mrModelPart(model_part)
	{
		mDelta = responseSettings["perturbation_size"].GetDouble();
	}

	/// Destructor.
	~SurfaceAreaResponseFunctionUtility() = default;

	///@}
	///@name Operators
	///@{

	///@}
	///@name Operations
	///@{

	// --------------------------------------------------------------------------
	double CalculateValue()
	{
		KRATOS_TRY;

		double total_area = 0.0;
		const std::size_t domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

		for (auto& cond_i : mrModelPart.Conditions()){
			total_area += cond_i.GetGeometry().Area();
		}

		return total_area;

		KRATOS_CATCH("");
	}

	// --------------------------------------------------------------------------
	void CalculateGradient()
	{
		KRATOS_TRY;

		// Formula computed in general notation:
		// \frac{da_{total}}{dx}

		VariableUtils().SetHistoricalVariableToZero(SHAPE_SENSITIVITY, mrModelPart.Nodes());

		const std::size_t domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

		for(auto& cond_i : mrModelPart.Conditions())
		{
			const double area_before_fd = cond_i.GetGeometry().Area();

			for (auto& node_i : cond_i.GetGeometry())
			{
				array_3d gradient(3, 0.0);

				double area_after_fd;
				node_i.X() += mDelta;
				node_i.X0() += mDelta;
				area_after_fd = cond_i.GetGeometry().Area();
				gradient[0] = (area_after_fd - area_before_fd) / mDelta;
				node_i.X() -= mDelta;
				node_i.X0() -= mDelta;

				// Apply pertubation in Y-direction and recompute total area of all neighbor elements
				area_after_fd = 0.0;
				node_i.Y() += mDelta;
				node_i.Y0() += mDelta;
				area_after_fd = cond_i.GetGeometry().Area();
				gradient[1] = (area_after_fd - area_before_fd) / mDelta;
				node_i.Y() -= mDelta;
				node_i.Y0() -= mDelta;

				// Apply pertubation in Z-direction and recompute total area of all neighbor elements
				area_after_fd = 0.0;
				node_i.Z() += mDelta;
				node_i.Z0() += mDelta;
				area_after_fd = cond_i.GetGeometry().Area();
				gradient[2] = (area_after_fd - area_before_fd) / mDelta;
				node_i.Z() -= mDelta;
				node_i.Z0() -= mDelta;

				// Compute sensitivity
				noalias(node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY)) += gradient;
			}
		}

		KRATOS_CATCH("");
	}

	// ==============================================================================

	///@}
	///@name Access
	///@{

	///@}
	///@name Inquiry
	///@{

	///@}
	///@name Input and output
	///@{

	/// Turn back information as a string.
	std::string Info() const
	{
		return "SurfaceAreaResponseFunctionUtility";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "SurfaceAreaResponseFunctionUtility";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
	{
	}

	///@}
	///@name Friends
	///@{

	///@}

protected:
	///@name Protected static Member Variables
	///@{

	///@}
	///@name Protected member Variables
	///@{

	///@}
	///@name Protected Operators
	///@{

	///@}
	///@name Protected Operations
	///@{

	// ==============================================================================

	///@}
	///@name Protected  Access
	///@{

	///@}
	///@name Protected Inquiry
	///@{

	///@}
	///@name Protected LifeCycle
	///@{

	///@}

private:
	///@name Static Member Variables
	///@{

	///@}
	///@name Member Variables
	///@{

	ModelPart &mrModelPart;
	double mDelta;

	///@}
///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{

	///@}
	///@name Private  Access
	///@{

	///@}
	///@name Private Inquiry
	///@{

	///@}
	///@name Un accessible methods
	///@{

	/// Assignment operator.
	//      SurfaceAreaResponseFunctionUtility& operator=(SurfaceAreaResponseFunctionUtility const& rOther);

	/// Copy constructor.
	//      SurfaceAreaResponseFunctionUtility(SurfaceAreaResponseFunctionUtility const& rOther);

	///@}

}; // Class SurfaceAreaResponseFunctionUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // SURFACE_AREA_RESPONSE_FUNCTION_UTILITY_H
