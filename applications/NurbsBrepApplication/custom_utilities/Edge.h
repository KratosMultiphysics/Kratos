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
#include "FaceTrim.h"

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

class Edge : public IndexedObject, public Flags
{
public:
  ///@name Type Definitions
  ///@{
  //typedef std::vector<double> DoubleVector;
  // For matrix / vector operations
  typedef std::vector<FaceTrim> FaceTrimVector;
  typedef std::vector<Vector> ParameterVector;

  ///@}

  /// Pointer definition of Edge
  //    KRATOS_CLASS_POINTER_DEFINITION[Edge];

  /// Default constructor.
  //TODO: make it "public IndexedObject, public Flags"
  Edge(unsigned int edge_id, ParameterVector boundary_vertices, FaceTrimVector face_trims_vector)
  : m_boundary_vertices(boundary_vertices),
    m_face_trims_vector(face_trims_vector),
    IndexedObject(edge_id),
    Flags()
  {
    //m_n_u = m_knot_vector_u.size() - m_p - 1;
  }

  /// Destructor.
  virtual ~Edge()
  {
  }
//TODO: you need to give reading access to your internals through the Calculate function
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
  //unsigned int m_edge_id;
  ParameterVector m_boundary_vertices;
  FaceTrimVector m_face_trims_vector;
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
