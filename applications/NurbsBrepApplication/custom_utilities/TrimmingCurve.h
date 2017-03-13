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

#include "nurbs_utilities.h"

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
  // Trimming Curve does not have a ID - it only carries an index, thus no use
  // of IndexedObject
class TrimmingCurve
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
    m_trim_index(trim_index)
  {
    //array_1d<double, 4> m_control_points;
    //m_n_u = m_knot_vector_u.size() - m_p - 1;
  }

  /// Destructor.
  virtual ~TrimmingCurve()
  {
  }
//TODO: you need to give reading access to your internals through the Calculate function


  std::vector<array_1d<double, 2>> CreatePolygon(unsigned int& number_polygon_points)
  {
    std::vector<array_1d<double, 2>> polygon;
    polygon.resize(number_polygon_points);

    // Variables needed
    unsigned int counter = 0;

    double u_min = m_knot_vector_u[0];
    double u_max = m_knot_vector_u[m_knot_vector_u.size() - 1];
    double delta_u = (u_max - u_min) / number_polygon_points;

    // Add points of edge to polygon
    double u_i = u_min;
    for (unsigned int i = 0; i<number_polygon_points; i++)
    {
      u_i += delta_u;
      Point<3> curve_point;

      EvaluateCurvePoint(curve_point, u_i);

      polygon[counter][0] = curve_point.X();
      polygon[counter][1] = curve_point.Y();

      counter++;
    }
    return polygon;
  }

  void EvaluateCurvePoint(Point<3>& rCurvePoint, double parameter_u)
  {
    const unsigned int R_Dim = 3;

    Vector tmpHomCtrlPt;
    Vector nBasis;
    Vector resulting_point;
    Vector homPoi;
    matrix<double> homCtrlPts;

    unsigned int span_u = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, parameter_u);
    NurbsUtilities::eval_nonzero_basis_function(nBasis, m_knot_vector_u, parameter_u, span_u, m_p);
    homCtrlPts.resize(R_Dim + 1, m_p + 1);
    for (unsigned int i = 0; i <= m_p; i++)
    {
      int control_point_index = span_u - m_p + i;

      // tmpHomCtrlPt = m_control_points[control_point_index].getWeight();
      // tmpHomCtrlPt = Ctrl_Pt[span-m_p+i]->get_Ctrl_Pt_Coo_w();
      Vector tmpHomCtrlPt(R_Dim + 1);

      double x = m_control_points[control_point_index][0];
      double y = m_control_points[control_point_index][1];
      double z = m_control_points[control_point_index][2];
      double w = m_control_points[control_point_index][3];

      tmpHomCtrlPt[0] = x*w;
      tmpHomCtrlPt[1] = y*w;
      tmpHomCtrlPt[2] = z*w;
      tmpHomCtrlPt[3] = w;

      for (unsigned int j = 0; j <= R_Dim; j++)
      {
        homCtrlPts(j, i) = tmpHomCtrlPt(j);
      }
    }
    homPoi = prod(homCtrlPts, nBasis);
    resulting_point = (1 / homPoi(R_Dim))*homPoi;

    rCurvePoint.X() = resulting_point[0];
    rCurvePoint.Y() = resulting_point[1];
    rCurvePoint.Z() = resulting_point[2];

    //if (std::abs(resulting_point(R_Dim) - 1.00) > m_epsilon)
    //{
    //  KRATOS_THROW_ERROR(std::logic_error, "NURBS 1D: evalutation curve point failed!!!", "");
    //}
  }



  unsigned int& GetIndex()
  {
    return m_trim_index;
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
  unsigned int m_trim_index;
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
