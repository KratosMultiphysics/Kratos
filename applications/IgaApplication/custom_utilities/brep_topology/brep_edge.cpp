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

    ///Constructor
    BrepEdge::BrepEdge(
        int rBrepId,
        std::vector<EdgeTopology>& rBrepEdgeTopologyVector,
        std::vector<TrimmingRange>& rTrimmingRangeVector,
        Vector& rKnotVector,
        Vector& rActiveRange,
        std::vector<int>& rControlPointIds,
        ModelPart& rModelPart)
        : m_brep_edge_topology_vector(rBrepEdgeTopologyVector),
          m_trimming_range_vector(rTrimmingRangeVector),
          m_active_range(rActiveRange),
          m_control_point_ids(rControlPointIds),
          m_model_part(rModelPart),
          IndexedObject(rBrepId),
          Flags()
    {
        //int number_of_nodes = rKnotVector.size() + rP - 1;

        //NodeCurveGeometry3D& m_nurbs_curve_geometry_3d = new NodeCurveGeometry3D(degree, NumberOfNodes);

        //for (int i = 0; i < rControlPointIds.size(); ++i)
        //{
        //    m_nurbs_curve_geometry_3d.SetNode(i, rModelPart.pGetNode(rControlPointIds[i]));
        //}

        //for (int i = 0; i < rKnotVector.size(); ++i)
        //{
        //    m_nurbs_curve_geometry_3d.SetPole(i, rKnotVector[i]);
        //}
    }
} // namespace Kratos.

