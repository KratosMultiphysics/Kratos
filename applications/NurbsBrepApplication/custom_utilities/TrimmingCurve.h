#ifndef TRIMMING_CURVE_H
#define TRIMMING_CURVE_H

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

class TrimmingCurve
{
public:
  ///@name Type Definitions
  ///@{

  // For matrix / vector operations
  typedef std::vector<double> DoubleVector;
  typedef std::vector<int> IntVector;

  ///@}

  /// Pointer definition of TrimmingCurve
  //    KRATOS_CLASS_POINTER_DEFINITION[TrimmingCurve];

  /// Default constructor.
  TrimmingCurve(unsigned int cp_id, DoubleVector knot_vector_u, unsigned int p, 
    IntVector control_points, DoubleVector boundary_vertices)
  : m_cp_id(cp_id),
    m_knot_vector_u(knot_vector_u),
    m_p(p),
    m_control_points(control_points),
    m_boundary_vertices(boundary_vertices)
  {
    //m_n_u = m_knot_vector_u.size() - m_p - 1;
  }

  /// Destructor.
  virtual ~TrimmingCurve()
  {
  }

  // --------------------------------------------------------------------------
  IntVector& GetControlPoints()
  {
    return m_control_points;
  }

  // ==============================================================================
  /// Turn back information as a string.
  virtual std::string Info() const
  {
    return "TrimmingCurve";
  }

  // ==============================================================================
  /// Print information about this object.
  virtual void PrintInfo(std::ostream &rOStream) const
  {
    rOStream << "TrimmingCurve";
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
  unsigned int m_cp_id;
  DoubleVector m_knot_vector_u;
  unsigned int m_p;
  IntVector m_control_points;
  DoubleVector m_boundary_vertices;
  //unsigned int m_n_u; // number of control points in u-direction
  //double m_epsilon = 1e-10; // Tolerance value

  // ==============================================================================
  // General working arrays
  // ==============================================================================
  /// Assignment operator.
  //      TrimmingCurve& operator=[TrimmingCurve const& rOther];

  /// Copy constructor.
  //      TrimmingCurve[TrimmingCurve const& rOther];

}; // Class TrimmingCurve

} // namespace Kratos.

#endif // TRIMMING_CURVE_H
