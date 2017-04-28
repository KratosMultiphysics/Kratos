//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Michael Breitenberger
//

// Project includes
#include "BrepEdge.h"
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"


namespace Kratos
{
  // --------------------------------------------------------------------------
  //std::vector<Node<3>::Pointer> BrepEdge::GetQuadraturePoints(const int& shapefunction_order)
  //{
  //  if (this->isCouplingEdge())
  //  {

  //  }
  //  else
  //  {

  //  }
  //}
  void BrepEdge::GetEdgeInformation(const int& face_trim, int& face_id, int& trim_index)
  {
    face_id = m_brep_face_trims_vector[face_trim].GetFaceId();
    trim_index = m_brep_face_trims_vector[face_trim].GetTrimIndex();
  }

  // --------------------------------------------------------------------------
  bool BrepEdge::isCouplingEdge()
  {
    if (m_brep_face_trims_vector.size() > 1)
      return true;
    else
      return false;
  }

  ///Constructor
  BrepEdge::BrepEdge(unsigned int edge_id, 
    ParameterVector& boundary_vertices, 
    BrepFaceTrimVector& brep_face_trims_vector)
    : m_boundary_vertices(boundary_vertices),
      m_brep_face_trims_vector(brep_face_trims_vector),
      IndexedObject(edge_id),
      Flags()
  {
  }

  ///Destructor
  BrepEdge::~BrepEdge()
  {
  }

}  // namespace Kratos.

