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
#include "brep_face.h"

namespace Kratos
{
    ///Constructor
    BrepFace::BrepFace(
        int& rBrepId,
        bool rIsTrimmed,
        bool rIsRational,
        std::vector<BrepTrimmingCurve>& rTrimmingLoops,
        std::vector<BrepBoundaryLoop>& rEmbeddedLoops,
        std::vector<EmbeddedPoint>& rEmbeddedPoints,
        Vector& rKnotVectorU,
        Vector& rKnotVectorV,
        int& rP,
        int& rQ,
        IntVector& rControlPointIds,
        Kratos::shared_ptr<ModelPart> rModelPart)
        : m_trimming_loops(rTrimmingLoops),
          m_is_trimmed(rIsTrimmed),
          m_is_rational(rIsRational),
          m_embedded_loops(rEmbeddedLoops),
          m_embedded_points(rEmbeddedPoints),
          mp_model_part(rModelPart),
          m_control_points_ids(rControlPointIds),
          IndexedObject(rBrepId),
          Flags()
    {
        //int number_of_nodes_u = rKnotVectorU.size() + rP - 1;
        //int number_of_nodes_v = rKnotVectorV.size() + rQ - 1;


        //NodeSurfaceGeometry3D& m_node_surface_geometry_3d(rP, rQ, number_of_nodes_u, number_of_nodes_v);

        //for (int i = 0; i < number_of_nodes_u; ++i)
        //{
        //    for (int j = 0; j < number_of_nodes_v; ++j)
        //    {
        //        m_node_surface_geometry_3d.SetNode(i, j, rModelPart->pGetNode(rControlPointIds[i*(number_of_nodes_u)+j]));
        //    }
        //}
    }

    ///Destructor
    BrepFace::~BrepFace()
    {}
} // namespace Kratos.