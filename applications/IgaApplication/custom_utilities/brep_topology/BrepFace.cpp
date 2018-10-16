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
#include "BrepFace.h"

namespace Kratos
{
    ///Constructor
    BrepFace::BrepFace(unsigned int rBrepId,
        bool rIsTrimmed,
        bool rIsRational,
        TrimmingLoopVector& rTrimmingLoops,
        TrimmingLoopVector& rEmbeddedLoops,
        std::vector<EmbeddedPoint>& rEmbeddedPoints,
        Vector& rKnotVectorU,
        Vector& rKnotVectorV,
        unsigned int& rP,
        unsigned int& rQ,
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
        int number_of_nodes_u = rKnotVectorU.size() + rP - 1;
        int number_of_nodes_v = rKnotVectorV.size() + rQ - 1;


        m_node_surface_geometry_3d = NodeSurfaceGeometry3D(rP, rQ, number_of_nodes_u, number_of_nodes_v);

        for (int i = 0; i < number_of_nodes_u; ++i)
        {
            for (int j = 0; j < number_of_nodes_v; ++j)
            {
                SetNode(i, j, rModelPart->pGetNode(rControlPointIds[i*(number_of_nodes_u)+j]));
            }
        }
    }

    ///Destructor
    BrepFace::~BrepFace()
    {}
} // namespace Kratos.