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
#include "BrepModel.h"
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"


namespace Kratos
{
  // --------------------------------------------------------------------------
  std::vector<BrepFace>& BrepModel::GetFaceVector()
  {
    return m_brep_faces;
  }
  std::vector<BrepEdge>& BrepModel::GetEdgeVector()
  {
    return m_brep_edges;
  }
  std::vector<BrepVertex>& BrepModel::GetVertexVector()
  {
	  return m_brep_vertices;
  }
  // --------------------------------------------------------------------------
  ///Constructor
  BrepModel::BrepModel(unsigned int& brep_id, std::vector<BrepFace>& faces,
	  std::vector<BrepEdge>& edges,
	  std::vector<BrepVertex>& vertices)
    : m_brep_faces(faces),
      m_brep_edges(edges),
	  m_brep_vertices(vertices),
      IndexedObject(brep_id),
      Flags()
  {
  }
  ///Destructor
  BrepModel::~BrepModel()
  {
  }
}  // namespace Kratos.

