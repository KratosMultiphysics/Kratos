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

    void BrepEdge::GetGeometryIntegration(ModelPart& rModelPart, 
        const std::string& rType,
        const std::string& rName,
        const int& rPropertiesId,
        const int& rShapeFunctionDerivativesOrder,
        std::vector<std::string> rVariables)
    {
        auto spans = m_node_curve_geometry_3d->Spans();

        Properties::Pointer this_property = rModelPart.pGetProperties(rPropertiesId);

        for (int i = 0; i < spans.size(); ++i)
        {
            ANurbs::Interval<double> domain(spans[i].T0(), spans[i].T1());

            int number_of_points = m_node_curve_geometry_3d->Degree() + 1;
            auto integration_points = ANurbs::IntegrationPoints<double>::Points1(number_of_points, domain);

            ANurbs::CurveShapeEvaluator<double> shape(m_node_curve_geometry_3d->Degree(), rShapeFunctionDerivativesOrder);

            for (int j = 0; j < integration_points.size(); ++j)
            {
                Vector local_coordinates(1);
                local_coordinates[0] = integration_points[j].t;

                shape.Compute(m_node_curve_geometry_3d->Knots(), integration_points[j].t);

                Element::GeometryType::PointsArrayType non_zero_control_points;
                std::vector<Node<3>::Pointer> cps;
                std::vector<int> cp_ids;
                for (int m = shape.FirstNonzeroPole(); m < shape.LastNonzeroPole() + 1; ++m)
                {
                    cp_ids.push_back(m_node_curve_geometry_3d->Node(m)->Id());
                    cps.push_back(m_node_curve_geometry_3d->Node(m));
                    non_zero_control_points.push_back(m_node_curve_geometry_3d->Node(m));
                }

                if (rType == "element")
                {
                    int id = rModelPart.GetRootModelPart().Elements().back().Id() + 1;
                    auto element = rModelPart.CreateNewElement(rName, id, non_zero_control_points, this_property);

                }
                if (rType == "condition")
                {
                    int id = rModelPart.GetRootModelPart().Conditions().back().Id() + 1;
                    auto condition = rModelPart.CreateNewCondition(rName, id, non_zero_control_points, this_property);
                }
            }
        }
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

