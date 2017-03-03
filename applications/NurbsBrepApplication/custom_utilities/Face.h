#ifndef FACE_H
#define FACE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
//#include "utilities/math_utils.h"

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
//#include "ControlPoint.h"
#include "TrimmingCurve.h"

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

//TODO: make it public IndexedObject, public Flags
class Face
{
public:
  ///@name Type Definitions
  ///@{
  typedef std::vector<double> DoubleVector;
  typedef std::vector<int> IntVector;
  typedef std::vector<std::vector<int>> TrimmingLoopVector;
  typedef std::vector<TrimmingCurve> TrimmingCurveVector;
  // For matrix vector operations
  //typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
  //typedef typename SparseSpaceType::MatrixType SparseMatrixType;
  //typedef typename SparseSpaceType::VectorType VectorType;
  //typedef std::vector<int> IntVector;
  //typedef std::vector<double> DoubleVector;
  //typedef boost::python::extract<double> takeDouble;
  //typedef boost::python::extract<int> takeInt;
  //typedef boost::python::extract<bool> takeBool;
  //typedef std::vector<ControlPoint> ControlPointVector;

  ///@}

  /// Pointer definition of Face
  //    KRATOS_CLASS_POINTER_DEFINITION[Face];


  /// Default constructor.
  Face(unsigned int brep_id, TrimmingCurveVector trimming_curves,
    TrimmingLoopVector trimming_loops,
    DoubleVector knot_vector_u, DoubleVector knot_vector_v,
    int p, int q, IntVector control_point_ids)
  : m_trimming_curves(trimming_curves),
    m_trimming_loops(trimming_loops),
    m_brep_id(brep_id),
    m_knot_vector_u(knot_vector_u),
    m_knot_vector_v(knot_vector_v),
    m_p(p),
    m_q(q),
    m_control_points_ids(control_point_ids)
  {
    unsigned int m_n_u = m_knot_vector_u.size() - m_p - 1;
    unsigned int m_n_v = m_knot_vector_v.size() - m_q - 1;

    if (m_control_points_ids.size() != m_n_u * m_n_v)
    {
      std::cout << "Invalid Face" << std::endl;
    }
  }

  /// Destructor.
  virtual ~Face()
  {
  }

//TODO: you need to give reading access to your internals through the Calculate function

  // ==============================================================================
  /// Turn back information as a string.
  virtual std::string Info() const
  {
    return "Face";
  }

  // ==============================================================================
  /// Print information about this object.
  virtual void PrintInfo(std::ostream &rOStream) const
  {
    rOStream << "Face";
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

  unsigned int m_brep_id;
  TrimmingCurveVector m_trimming_curves;
  TrimmingLoopVector m_trimming_loops;
  DoubleVector m_knot_vector_u;
  DoubleVector m_knot_vector_v;
  unsigned int m_p;
  unsigned int m_q;
  IntVector m_control_points_ids;

  // ==============================================================================
  // General working arrays
  // ==============================================================================
  /// Assignment operator.
  //      Face& operator=[Face const& rOther];

  /// Copy constructor.
  //      Face[Face const& rOther];

}; // Class Face

} // namespace Kratos.

#endif // FACE_H
