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
    //std::cout << "m_faces.size(): " << m_faces.size() << std::endl;
    return m_brep_faces;
  }

  std::vector<BrepEdge>& BrepModel::GetEdgeVector()
  {
    return m_brep_edges;
  }
  // --------------------------------------------------------------------------
  ///Constructor
  BrepModel::BrepModel(unsigned int& brep_id, FacesVector& faces, EdgesVector& edges)
    : m_brep_faces(faces),
    m_brep_edges(edges),
    IndexedObject(brep_id),
    Flags()
  {
    //std::cout << "m_faces.size(): " << m_faces.size() << std::endl;
  }

  ///Destructor
  BrepModel::~BrepModel()
  {
  }

}  // namespace Kratos.

