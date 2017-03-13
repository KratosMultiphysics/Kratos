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
#include "Edge.h"
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"


namespace Kratos
{
  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  ///Constructor
  Edge::Edge(unsigned int edge_id, 
    ParameterVector& boundary_vertices, 
    FaceTrimVector& face_trims_vector)
    : m_boundary_vertices(boundary_vertices),
      m_face_trims_vector(face_trims_vector),
      IndexedObject(edge_id),
      Flags()
  {
  }

  ///Destructor
  Edge::~Edge()
  {
  }

}  // namespace Kratos.

