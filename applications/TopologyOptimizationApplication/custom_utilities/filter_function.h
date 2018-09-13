// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

#if !defined(KRATOS_FILTER_FUNCTION_H_INCLUDED)
#define  KRATOS_FILTER_FUNCTION_H_INCLUDED

// System includes
#include <iostream>
#include <string>
#include <algorithm>
#include <iomanip>      // for std::setprecision

// External includes
#include <boost/python.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/process_info.h"


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

/// Solution utility to filter results.
/** Detail class definition.

 */

class FilterFunction
{
public:

	///@name Type Definitions
	///@{

	/// Pointer definition of FilterFunction
	KRATOS_CLASS_POINTER_DEFINITION(FilterFunction);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	FilterFunction( std::string FilterFunctionType, const double SearchRadius )
	:mSearchRadius(SearchRadius)
	{
		// Set precision for output
		std::cout.precision(12);

		// Set type of weighting function

		// Type 1: Gaussian function
		std::string linear("linear");
		if(FilterFunctionType.compare(linear)==0)
			mFilterFunctionType = 1;
	}
	/// Destructor.
	virtual ~FilterFunction()
	{
	}


	///@}
	///@name Operators
	///@{


	///@}
	///@name Operations
	///@{

	double ComputeWeight( const double distance )
	{
		KRATOS_TRY;

		// Depending on which weighting function is chosen, compute weight
		double weight = 0.0;
		switch(mFilterFunctionType)
		{
		// Compute weight for linear filter function
		case 1:
		{
			weight = std::max(0.0,mSearchRadius - distance);
			break;
		}

		}

		return weight;

		KRATOS_CATCH("");
	}

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
	virtual std::string Info() const
	{
		return "FilterFunction";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const
	{
		rOStream << "FilterFunction";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream& rOStream) const
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

	const double mSearchRadius;
	unsigned int mFilterFunctionType;

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
	//FilterFunction& operator=(FilterFunction const& rOther);

	/// Copy constructor.
	//FilterFunction(FilterFunction const& rOther);


	///@}

}; // Class FilterFunction

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif	/* KRATOS_FILTER_FUNCTION_H_INCLUDED */
