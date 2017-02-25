#ifndef FACE_TRIM_H
#define FACE_TRIM_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
//#include "ControlPoint.h"
//#include "b_spline_utilities.h"

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

class FaceTrim
{
public:
	///@name Type Definitions
	///@{

	// For matrix / vector operations
	typedef std::vector<double> DoubleVector;
	typedef std::vector<std::vector<double>> ParameterVector;

	///@}

	/// Pointer definition of Trim
	//    KRATOS_CLASS_POINTER_DEFINITION[FaceTrim];

	/// Default constructor.
	FaceTrim(unsigned int face_id, unsigned int trim_index, 
		ParameterVector boundary_parameters, bool relative_direction)
	: m_face_id(face_id),
	  m_trim_index(trim_index),
	  m_boundary_parameters(boundary_parameters),
	  m_relative_direction(relative_direction)
	{
		//m_n_u = m_knot_vector_u.size() - m_p - 1;
	}

	/// Destructor.
	virtual ~FaceTrim()
	{
	}

	// ==============================================================================
	/// Turn back information as a string.
	virtual std::string Info() const
	{
		return "FaceTrim";
	}

	// ==============================================================================
	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "FaceTrim";
	}

	// ==============================================================================
	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
	{
	}


private:
	// ==============================================================================
	// Initialized by class constructor
	// ==============================================================================
	unsigned int m_face_id;
	unsigned int m_trim_index;
	ParameterVector m_boundary_parameters;
	bool m_relative_direction;

	// ==============================================================================
	// General working arrays
	// ==============================================================================
	/// Assignment operator.
	//      FaceTrim& operator=[FaceTrim const& rOther];

	/// Copy constructor.
	//      FaceTrim[FaceTrim const& rOther];

}; // Class FaceTrim

} // namespace Kratos.

#endif // FACE_TRIM_H
