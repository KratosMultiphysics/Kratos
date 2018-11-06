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
//                   Thomas Oberbichler
//


// Project includes
#include "brep_edge.h"

namespace Kratos
{
    /* If edge has more then one brep topology object, then the edge is a coupling
    *  edge.
    */
    bool BrepEdge::IsCouplingEdge()
    {
        if (m_brep_edge_topology_vector.size() > 1)
            return true;
        else
            return false;
    }

    std::vector<Node<3>> BrepEdge::GetGeometryIntegrationNodes()
    {
        std::vector<Node<3>> integration_nodes;



        return integration_nodes;
    }

    ///Constructor
    BrepEdge::BrepEdge(
        int rBrepId,
        std::vector<EdgeTopology>& rBrepEdgeTopologyVector,
        std::vector<TrimmingRange>& rTrimmingRangeVector,
        std::vector<EmbeddedPoint>& rEmbeddedPoints,
        int& rDegree,
        Vector& rKnotVector,
        Vector& rActiveRange,
        std::vector<int>& rControlPointIds,
        ModelPart& rModelPart)
        : m_brep_edge_topology_vector(rBrepEdgeTopologyVector),
          m_trimming_range_vector(rTrimmingRangeVector),
          m_embedded_points(rEmbeddedPoints),
          m_active_range(rActiveRange),
          m_control_point_ids(rControlPointIds),
          m_model_part(rModelPart),
          IndexedObject(rBrepId),
          Flags()
    {
        int number_of_nodes = rKnotVector.size() + rDegree - 1;

        m_node_curve_geometry_3d = Kratos::make_unique<NodeCurveGeometry3D>(
            rDegree, number_of_nodes);

        for (int i = 0; i < rControlPointIds.size(); ++i)
        {
            Node<3>::Pointer node = rModelPart.pGetNode(rControlPointIds[i]);
            m_node_curve_geometry_3d->SetNode(i, node);
            m_node_curve_geometry_3d->SetWeight(i, node->GetValue(NURBS_CONTROL_POINT_WEIGHT));
        }

        for (int i = 0; i < rKnotVector.size() - 2; ++i)
        {
            m_node_curve_geometry_3d->SetKnot(i, rKnotVector[i + 1]);
        }
    }
} // namespace Kratos.

