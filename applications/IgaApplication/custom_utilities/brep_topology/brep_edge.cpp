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
    void BrepEdge::GetGeometryNodes(
        ModelPart& rModelPart,
        const int& rT) const
    {
        int number_of_cps = mNodeCurveGeometry3D->NbPoles();

        int t_start = 0;
        int t_end = number_of_cps;

        if (rT >= 0)
        {
            t_start = rT * (number_of_cps - 1);
            t_end = rT * (number_of_cps - 1) + 1;
        }

        for (int i = t_start; i < t_end; ++i)
        {
            rModelPart.AddNode(mNodeCurveGeometry3D->GetNode(i));
        }
    }

    void BrepEdge::GetGeometryVariationNodes(
        ModelPart& rModelPart,
        const int& rT) const
    {
        int number_of_cps = mNodeCurveGeometry3D->NbPoles();

        KRATOS_ERROR_IF(number_of_cps < 3) << "BrepEdge::GetGeometryVariationNodes: Not enough control points to get Variation Nodes." << std::endl;

        int t_start = 1;
        int t_end = number_of_cps - 1;

        if (rT == 0)
        {
            rModelPart.AddNode(mNodeCurveGeometry3D->GetNode(1));
        }
        else if (rT == 1)
        {
            rModelPart.AddNode(mNodeCurveGeometry3D->GetNode(number_of_cps - 1));
        }
    }

    bool BrepEdge::IsCouplingEdge()
    {
        return (mBrepEdgeTopologyVector.size() > 1);
    }

    const BrepEdge::EdgeTopology BrepEdge::GetEdgeTopology(
        const int rTopologyIndex) const
    {
        KRATOS_ERROR_IF(rTopologyIndex >= mBrepEdgeTopologyVector.size())
            << "BrepEdge::GetEdgeTopology: Number of topology references smaller than selected index! Selected index: "
            << rTopologyIndex << ", size of edge topologies: " << mBrepEdgeTopologyVector.size() << std::endl;

        return mBrepEdgeTopologyVector[rTopologyIndex];
    }

    const int BrepEdge::GetNumberOfEdgeTopologies() const
    {
        return mBrepEdgeTopologyVector.size();
    }

    void BrepEdge::GetIntegrationGeometry(ModelPart& rModelPart,
        const std::string& rType,
        const std::string& rName,
        const int& rShapeFunctionDerivativesOrder,
        std::vector<std::string> rVariables) const
    {
        //for (int trims = 0; trims < m_trimming_range_vector.size(); ++trims)
        //{
            //auto this_curve = Kratos::make_unique<Kratos::Curve<3>>(
                //mNodeCurveGeometry3D, m_trimming_range_vector[trims].range);
            auto spans = mNodeCurveGeometry3D->Spans();

            for (int i = 0; i < spans.size(); ++i)
            {
                ANurbs::Interval<double> domain(spans[i].T0(), spans[i].T1());

                int number_of_points = mNodeCurveGeometry3D->Degree() + 1;
                auto integration_points = ANurbs::IntegrationPoints<double>::Points1(number_of_points, domain);

                ANurbs::CurveShapeEvaluator<double> shape(mNodeCurveGeometry3D->Degree(), rShapeFunctionDerivativesOrder);

                for (int j = 0; j < integration_points.size(); ++j)
                {
                    Vector local_coordinates(1);
                    local_coordinates[0] = integration_points[j].t;

                    shape.Compute(mNodeCurveGeometry3D->Knots(), integration_points[j].t);

                    Vector N_0 = ZeroVector(shape.NbNonzeroPoles());
                    Matrix N_1 = ZeroMatrix(shape.NbNonzeroPoles(), 1);

                    Element::GeometryType::PointsArrayType non_zero_control_points;
                    for (int m = shape.FirstNonzeroPole(); m < shape.LastNonzeroPole() + 1; ++m)
                    {
                        non_zero_control_points.push_back(mNodeCurveGeometry3D->GetNode(m));

                        N_0(m) = shape(0, m);

                        N_1(m, 0) = shape(1, m);
                    }
                    if (rType == "element")
                    {
                        int id = 0;
                        if (rModelPart.GetRootModelPart().Elements().size() > 0)
                            id = rModelPart.GetRootModelPart().Elements().back().Id() + 1;

                        auto element = rModelPart.CreateNewElement(rName, id, non_zero_control_points, 0);

                        element->SetValue(SHAPE_FUNCTION_VALUES, N_0);
                        element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, N_1);
                        element->SetValue(INTEGRATION_WEIGHT, integration_points[j].weight);
                    }
                    if (rType == "condition")
                    {
                        int id = 0;
                        if (rModelPart.GetRootModelPart().Conditions().size() > 0)
                            int id = rModelPart.GetRootModelPart().Conditions().back().Id() + 1;
                        auto condition = rModelPart.CreateNewCondition(rName, id, non_zero_control_points, 0);

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
        const int& rShapeFunctionDerivativesOrder,
        std::vector<std::string> rVariables) const
    {
        for (int ep = 0; ep < mEmbeddedPoints.size(); ++ep)
        {
            if (mEmbeddedPoints[ep].trim_index == trim_index)
            {
                ANurbs::CurveShapeEvaluator<double> shape(
                    mNodeCurveGeometry3D->Degree(),
                    rShapeFunctionDerivativesOrder);

                shape.Compute(
                    mNodeCurveGeometry3D->Knots(),
                    mEmbeddedPoints[ep].local_parameter);
            }
        }
    }

    const Kratos::shared_ptr<NodeCurveGeometry3D> BrepEdge::GetCurve3d() const
    {
        return mNodeCurveGeometry3D;
    }

    const std::vector<BrepEdge::EdgeTopology>& BrepEdge::GetBrepEdgeTopologyVector() const
    {
        return mBrepEdgeTopologyVector; 
    }

    ///Constructor
    BrepEdge::BrepEdge(
        const int rBrepId,
        std::vector<EdgeTopology>& rBrepEdgeTopologyVector,
        std::vector<TrimmingRange>& rTrimmingRangeVector,
        std::vector<EmbeddedPoint>& rEmbeddedPoints,
        const int rDegree,
        Vector& rKnotVector,
        Vector& rActiveRange,
        std::vector<int>& rControlPointIds,
        ModelPart& rModelPart)
        : mBrepEdgeTopologyVector(rBrepEdgeTopologyVector),
          mTrimmingRangeVector(rTrimmingRangeVector),
          mEmbeddedPoints(rEmbeddedPoints),
          mActiveRange(rActiveRange),
          mModelPart(rModelPart),
          IndexedObject(rBrepId),
          Flags()
    {
        int number_of_nodes = rControlPointIds.size();

        mNodeCurveGeometry3D = Kratos::make_unique<NodeCurveGeometry3D>(
            rDegree, number_of_nodes);

        for (int i = 0; i < rControlPointIds.size(); ++i)
        {
            Node<3>::Pointer node = rModelPart.pGetNode(rControlPointIds[i]);
            mNodeCurveGeometry3D->SetNode(i, node);
            mNodeCurveGeometry3D->SetWeight(i, node->GetValue(NURBS_CONTROL_POINT_WEIGHT));
        }

        for (int i = 0; i < rKnotVector.size() - 2; ++i)
        {
            mNodeCurveGeometry3D->SetKnot(i, rKnotVector[i + 1]);
        }
    }
} // namespace Kratos.

