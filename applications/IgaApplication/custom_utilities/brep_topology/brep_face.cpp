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
        ModelPart& rModelPart)
        : m_trimming_loops(rTrimmingLoops),
          m_is_trimmed(rIsTrimmed),
          m_is_rational(rIsRational),
          m_embedded_loops(rEmbeddedLoops),
          m_embedded_points(rEmbeddedPoints),
          m_model_part(rModelPart),
          m_control_points_ids(rControlPointIds),
          IndexedObject(rBrepId),
          Flags()
    {
        int number_of_nodes_u = rKnotVectorU.size() + rP - 1;
        int number_of_nodes_v = rKnotVectorV.size() + rQ - 1;

        m_node_surface_geometry_3d = New<NodeSurfaceGeometry3D>(
            rP, rQ, number_of_nodes_u, number_of_nodes_v);

        for (int i = 0; i < rKnotVectorU.size() - 2; ++i)
        {
            m_node_surface_geometry_3d->SetKnotU(i, rKnotVectorU(i + 1));
        }

        for (int i = 0; i < rKnotVectorV.size() - 2; ++i)
        {
            m_node_surface_geometry_3d->SetKnotV(i, rKnotVectorV(i + 1));
        }

        for (int i = 0; i < number_of_nodes_u; ++i)
        {
            for (int j = 0; j < number_of_nodes_v; ++j)
            {
                Node<3>::Pointer node = rModelPart.pGetNode(rControlPointIds[i*(number_of_nodes_u)+j]);
                m_node_surface_geometry_3d->SetNode(i, j, node);
                if (rIsRational)
                {
                    m_node_surface_geometry_3d->SetWeight(i, j, node->GetValue(NURBS_CONTROL_POINT_WEIGHT));
                }
            }
        }
    }
} // namespace Kratos.