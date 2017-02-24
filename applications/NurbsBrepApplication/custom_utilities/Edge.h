#ifndef EDGE_H
#define EDGE_H

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
#include "ControlPoint.h"
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

class Edge
{
public:
	///@name Type Definitions
	///@{
	typedef std::vector<double> DoubleVector;
	// For matrix / vector operations
	typedef std::vector<ControlPoint> ControlPointVector;

	///@}

	/// Pointer definition of Edge
	//    KRATOS_CLASS_POINTER_DEFINITION[Edge];

	/// Default constructor.
	Edge(DoubleVector knot_vector_u, unsigned int p, ControlPointVector control_points) 
	: m_knot_vector_u(knot_vector_u),
	  m_p(p),
	  m_control_points(control_points)
	{
		//m_n_u = m_knot_vector_u.size() - m_p - 1;
	}

	/// Destructor.
	virtual ~Edge()
	{
	}

	// --------------------------------------------------------------------------
	ControlPointVector& GetControlPoints()
	{
		return m_control_points;
	}

	// ==============================================================================
	/// Turn back information as a string.
	virtual std::string Info() const
	{
		return "Edge";
	}

	// ==============================================================================
	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "Edge";
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
	DoubleVector m_knot_vector_u;
	unsigned int m_p;
	ControlPointVector m_control_points;
	//unsigned int m_n_u; // number of control points in u-direction
	//double m_epsilon = 1e-10; // Tolerance value

	// ==============================================================================
	// General working arrays
	// ==============================================================================
	/// Assignment operator.
	//      Edge& operator=[Edge const& rOther];

	/// Copy constructor.
	//      Edge[Edge const& rOther];

}; // Class Edge

} // namespace Kratos.

#endif // EDGE_H
