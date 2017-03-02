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

//// ------------------------------------------------------------------------------
//// External includes
//// ------------------------------------------------------------------------------
//#include <boost/python.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/io.hpp>
// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"

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

//TODO: how is this different from an edge? shouldn't you derive it from edge?
//TODO: make it "public IndexedObject, public Flags"
class TrimmingCurve : public IndexedObject, public Flags
{
public:
  ///@name Type Definitions
  ///@{

  // For matrix / vector operations
  //typedef std::vector<double> DoubleVector;
  typedef std::vector<array_1d<double, 4>> ControlPointVector;

  ///@}

  /// Pointer definition of TrimmingCurve
  //    KRATOS_CLASS_POINTER_DEFINITION[TrimmingCurve];

  /// Default constructor.
//TODO: pass by reference not by value
//TODO: why control points have size 4??? pass them as kratos nodes
  TrimmingCurve(unsigned int trim_index, bool curve_direction, Vector& knot_vector_u,
    unsigned int p, ControlPointVector& control_points, 
    Vector& active_range)
  : m_knot_vector_u(knot_vector_u),
    m_curve_direction(curve_direction),
    m_p(p),
    m_control_points(control_points),
    m_active_range(active_range),
    IndexedObject(trim_index),
    Flags()
  {
    //array_1d<double, 4> m_control_points;
    //m_n_u = m_knot_vector_u.size() - m_p - 1;
  }

  /// Destructor.
  virtual ~TrimmingCurve()
  {
  }
//TODO: you need to give reading access to your internals through the Calculate function

  // --------------------------------------------------------------------------
  //ControlPointVector& GetControlPoints()
  //{
  //  return m_control_points;
  //}

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
  //unsigned int m_cp_id;
  bool m_curve_direction;
  Vector m_knot_vector_u;
  unsigned int m_p;
  ControlPointVector m_control_points;
  Vector m_active_range;
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
