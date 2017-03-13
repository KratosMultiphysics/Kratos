#ifndef BREP_MODEL_H
#define BREP_MODEL_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <vector>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "Face.h"
#include "Edge.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

  //typedef std::vector<Vertex> VerticesVector;

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

class BrepModel : public IndexedObject, public Flags
{
public:
  ///@name Type Definitions
  ///@{
  typedef std::vector<Face> FacesVector;
  typedef std::vector<Edge> EdgesVector;



  ///@}
  /// Pointer definition of BrepModel
  //    KRATOS_CLASS_POINTER_DEFINITION[BrepModel];

  /// Default constructor.
//TODO: pass by reference not by value
  BrepModel(unsigned int& brep_id, FacesVector& faces, EdgesVector& edges)
  : m_faces(faces),
    m_edges(edges),
    IndexedObject(brep_id),
    Flags()
  {
    std::cout << "m_faces.size(): " << m_faces.size() << std::endl;
  }
  /// Destructor.
  virtual ~BrepModel()
  {
  }
  //TODO: delete explicitly the copy constructor
  /// Copy constructor.
  //BrepModel(const BrepModel&) = delete;
  /// Copy constructor.
  //BrepModel[BrepModel const& rOther];
  // ==============================================================================
  FacesVector& GetFaceVector()
  {
    std::cout << "m_faces.size(): " << m_faces.size() << std::endl;
    return m_faces;
  }

  EdgesVector& GetEdgeVector()
  {
    return m_edges;
  }
  // ==============================================================================
  /*Face& GetFace(const unsigned int& face_id, bool& is_existing)
  {
    is_existing = false;
    for (unsigned int i = 0; i < m_faces.size(); i++)
    {
      if (m_faces[i].GetId() == face_id)
      {
        is_existing = true;
        return m_faces[i];
      }
    }
    return m_faces[0];
  }*/

  // ==============================================================================
  /// Turn back information as a string.
  virtual std::string Info() const
  {
    return "BrepModel";
  }

  // ==============================================================================
  /// Print information about this object.
  virtual void PrintInfo(std::ostream &rOStream) const
  {
    rOStream << "BrepModel";
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
  FacesVector m_faces;
  EdgesVector m_edges;
  //VerticesVector m_vertices;

  // ==============================================================================
  // General working arrays
  // ==============================================================================
  /// Assignment operator.
  //      BrepModel& operator=[BrepModel const& rOther];



}; // Class BrepModel

} // namespace Kratos.

#endif // BREP_MODEL_H
