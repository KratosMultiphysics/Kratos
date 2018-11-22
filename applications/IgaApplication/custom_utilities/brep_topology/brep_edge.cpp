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

    const BrepEdge::EdgeTopology BrepEdge::GetEdgeTopology(
        const int& rTopologyIndex) const
    {
        KRATOS_ERROR_IF(rTopologyIndex >= m_brep_edge_topology_vector.size()) 
            << "Number of topology references smaller than selected index!" << std::endl;

        return m_brep_edge_topology_vector[rTopologyIndex];
    }

    void BrepEdge::GetIntegrationGeometry(ModelPart& rModelPart,
        const std::string& rType,
        const std::string& rName,
        const int& rPropertiesId,
        const int& rShapeFunctionDerivativesOrder,
        std::vector<std::string> rVariables)
    {
        Properties::Pointer this_property = rModelPart.pGetProperties(rPropertiesId);

        //for (int trims = 0; trims < m_trimming_range_vector.size(); ++trims)
        //{
            //auto this_curve = Kratos::make_unique<Kratos::Curve<3>>(
                //m_node_curve_geometry_3d, m_trimming_range_vector[trims].range);
            auto spans = m_node_curve_geometry_3d->Spans();

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

                    Vector N_0 = ZeroVector(shape.NbNonzeroPoles());
                    Matrix N_1 = ZeroMatrix(shape.NbNonzeroPoles(), 1);

                    Element::GeometryType::PointsArrayType non_zero_control_points;
                    for (int m = shape.FirstNonzeroPole(); m < shape.LastNonzeroPole() + 1; ++m)
                    {
                        non_zero_control_points.push_back(m_node_curve_geometry_3d->Node(m));

                        N_0(m) = shape(0, m);

                        N_1(m, 0) = shape(1, m);
                    }
                    if (rType == "element")
                    {
                        int id = 0;
                        if (rModelPart.GetRootModelPart().Elements().size() > 0)
                            id = rModelPart.GetRootModelPart().Elements().back().Id() + 1;

                        auto element = rModelPart.CreateNewElement(rName, id, non_zero_control_points, this_property);

                        element->SetValue(SHAPE_FUNCTION_VALUES, N_0);
                        element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, N_1);
                        element->SetValue(INTEGRATION_WEIGHT, integration_points[j].weight);
                    }
                    if (rType == "condition")
                    {
                        int id = 0;
                        if (rModelPart.GetRootModelPart().Conditions().size() > 0)
                            int id = rModelPart.GetRootModelPart().Conditions().back().Id() + 1;
                        auto condition = rModelPart.CreateNewCondition(rName, id, non_zero_control_points, this_property);

                        condition->SetValue(SHAPE_FUNCTION_VALUES, N_0);
                        condition->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, N_1);
                        condition->SetValue(INTEGRATION_WEIGHT, integration_points[j].weight);
                    }
                    // for strong application of properties on control point nodes
                    if (rType == "node")
                    {
                        for (int sh_nonzero = 0; sh_nonzero < N_0.size(); ++sh_nonzero)
                        {
                            if (N_0(sh_nonzero) > 0.0)
                            {
                                rModelPart.AddNode(non_zero_control_points(sh_nonzero));
                            }
                        }
                    }
                }
            }
        //}
    }

    void BrepEdge::GetIntegrationBrep(
        ModelPart& rModelPart,
        const int& trim_index,
        const std::string& rType,
        const std::string& rName,
        const int& rPropertiesId,
        const int& rShapeFunctionDerivativesOrder,
        std::vector<std::string> rVariables)
    {
        for (int ep = 0; ep < m_embedded_points.size(); ++ep)
        {
            if (m_embedded_points[ep].trim_index == trim_index)
            {
                ANurbs::CurveShapeEvaluator<double> shape(
                    m_node_curve_geometry_3d->Degree(), 
                    rShapeFunctionDerivativesOrder);

                shape.Compute(m_node_curve_geometry_3d->Knots(), m_embedded_points[ep].local_parameter);
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
        int number_of_nodes = rControlPointIds.size();

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

