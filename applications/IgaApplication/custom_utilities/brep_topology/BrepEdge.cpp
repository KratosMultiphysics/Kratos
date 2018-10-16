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
#include "BrepEdge.h"

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
    BrepEdge::BrepEdge(unsigned int rBrepId,
		std::vector<Topology>& rBrepEdgeTopologyVector,
		std::vector<TrimmingRange>& rTrimmingRangeVector,
		Vector& rKnotVector,
		Vector& rActiveRange,
		std::vector<int>& rControlPointIds,
		Kratos::shared_ptr<ModelPart> rModelPart)
        : m_brep_edge_topology_vector(rBrepEdgeTopologyVector),
          m_trimming_range_vector(rTrimmingRangeVector),
          m_active_range(rActiveRange),
          m_control_point_ids(rControlPointIds),
          mp_model_part(rModelPart),
          IndexedObject(rBrepId),
          Flags()
    {
        m_curve_geometry = NodeCurveGeometry(degree, NumberOfNodes);

        for (int i = 0; i < control_point_ids.size(); ++i)
        {
            m_curve_geometry.SetNode(i, mp_model_part->pGetNode(control_point_ids[i]));
        }

        for (int i = 0; i < knot_vector.size(); ++i)
        {
            m_curve_geometry.SetPole(i, knot_vector[i]);
        }
    }
} // namespace Kratos.

