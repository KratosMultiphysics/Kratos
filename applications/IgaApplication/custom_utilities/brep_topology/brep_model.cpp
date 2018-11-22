//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:    BSD License
//              Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Michael Breitenberger
//                   Thomas Oberbichler
//

// Project includes
#include "brep_model.h"

#include "iga_application_variables.h"


namespace Kratos
{
    bool BrepModel::GetIntegrationDomainGeometry(
        ModelPart& rModelPart, 
        int& brep_id,
        const std::string& rType,
        const std::string& rName,
        const int& rPropertiesId,
        const int& rShapeFunctionDerivativesOrder,
        std::vector<std::string> rVariables)
    {
        bool success = false;

        for (int i = 0; i < m_brep_faces.size(); ++i)
        {
            if (m_brep_faces[i].Id() == brep_id)
            {
                m_brep_faces[i].GetGeometryIntegration(
                    rModelPart, rType, rName, rPropertiesId, 
                    rShapeFunctionDerivativesOrder, rVariables);
                return true;
            }
        }

        for (int i = 0; i < m_brep_edges.size(); ++i)
        {
            if (m_brep_edges[i].Id() == brep_id)
            {
                m_brep_faces[i].GetGeometryIntegration(
                    rModelPart, rType, rName, rPropertiesId, 
                    rShapeFunctionDerivativesOrder, rVariables);
                return true;
            }
        }

        return success;
    }

    bool BrepModel::GetIntegrationDomainBrep(
        ModelPart& rModelPart, 
        int& brep_id,
        const std::string& rType,
        const std::string& rName,
        const int& rPropertiesId,
        const int& rShapeFunctionDerivativesOrder,
        std::vector<std::string> rVariables)
    {
        bool success = false;
        for (int i = 0; i < m_brep_edges.size(); ++i)
        {
            if (m_brep_edges[i].Id() == brep_id)
            {
                BrepEdge::EdgeTopology edge_topology = m_brep_edges[i].GetEdgeTopology(0);

                for (int j = 0; j < m_brep_faces.size(); ++j)
                {
                    if (m_brep_faces[j].Id() == edge_topology.brep_id)
                    {
                        m_brep_faces[j].GetIntegrationBrepEdge(
                            rModelPart, edge_topology.trim_index,
                            rType, rName, rPropertiesId,
                            rShapeFunctionDerivativesOrder,
                            rVariables);
                        return true;
                    }
                }
            }
        }

        for (int i = 0; i < m_brep_vertices.size(); ++i)
        {
            if (m_brep_vertices[i].Id() == brep_id)
            {
                BrepVertex::VertexTopology vertex_topology = m_brep_vertices[i].GetVertexTopology(0);

                for (int j = 0; j < m_brep_edges.size(); ++j)
                {
                    if (m_brep_edges[j].Id() == vertex_topology.brep_id)
                    {
                        m_brep_edges[j].GetIntegrationBrep(
                            rModelPart, 
                            vertex_topology.trim_index, 
                            rType, rName, rPropertiesId, 
                            rShapeFunctionDerivativesOrder, 
                            rVariables);
                        return true;
                    }
                }
            }
        }

        return success;
    }

    bool BrepModel::GetIntegrationDomainBrepCoupling(
        ModelPart& rModelPart,
        const std::string& rType,
        const std::string& rName,
        const int& rPropertiesId,
        const int& rShapeFunctionDerivativesOrder,
        std::vector<std::string> rVariables)
    {
        bool success = false;

        for (int i = 0; i < m_brep_edges.size(); ++i)
        {
            if (m_brep_edges[i].IsCouplingEdge())
            {
                auto master = m_brep_edges[i].GetEdgeTopology(0);
                auto slave = m_brep_edges[i].GetEdgeTopology(1);

                GetIntegrationBrepCouplingEdge(
                    master, slave,
                    rModelPart,
                    rType,
                    rName,
                    rPropertiesId,
                    rShapeFunctionDerivativesOrder,
                    rVariables);

                success = true;
            }
        }
        return success;
    }

    bool BrepModel::GetIntegrationDomainBrepCoupling(
        ModelPart& rModelPart,
        int& brep_id,
        const std::string& rType,
        const std::string& rName,
        const int& rPropertiesId,
        const int& rShapeFunctionDerivativesOrder,
        std::vector<std::string> rVariables)
    {
        bool success = false;
        for (int i = 0; i < m_brep_edges.size(); ++i)
        {
            if (m_brep_edges[i].Id() == brep_id)
            {
                BrepEdge::EdgeTopology edge_topology = m_brep_edges[i].GetEdgeTopology(0);
                if (m_brep_faces[i].Id() == edge_topology.brep_id)
                {
                    m_brep_faces[i].GetIntegrationBrepEdge(
                        rModelPart, 
                        edge_topology.trim_index,
                        rType, 
                        rName, 
                        rPropertiesId,
                        rShapeFunctionDerivativesOrder,
                        rVariables);
                    return true;
                }
            }
        }
        return success;
    }


    void BrepModel::GetIntegrationBrepCouplingEdge(
        const BrepEdge::EdgeTopology& master,
        const BrepEdge::EdgeTopology& slave,
        ModelPart& rModelPart,
        const std::string& rType,
        const std::string& rName,
        const int& rPropertiesId,
        const int& rShapeFunctionDerivativesOrder,
        std::vector<std::string> rVariables) const
    {
        Properties::Pointer this_property = rModelPart.pGetProperties(rPropertiesId);

        const BrepFace& face_1 = GetFace(master.brep_id);
        auto surface_1 = face_1.GetSurface();
        auto curve_2d_1 = face_1.GetTrimCurve(master.trim_index);

        int degree_1 = surface_1->DegreeU() + surface_1->DegreeV();

        const BrepFace& face_2 = GetFace(slave.brep_id);
        auto surface_2 = face_2.GetSurface();
        auto curve_2d_2 = face_2.GetTrimCurve(slave.trim_index);

        int degree_2 = surface_1->DegreeU() + surface_1->DegreeV();

        int number_of_points_per_knot_span = std::max(degree_1, degree_2) + 1;

        auto curve_on_surface_3d_1 = Kratos::make_shared<CurveOnSurface<3>>(
            curve_2d_1->CurveGeometry(), surface_1, curve_2d_1->Domain());
        auto curve_on_surface_3d_2 = Kratos::make_shared<CurveOnSurface<3>>(
            curve_2d_2->CurveGeometry(), surface_2, curve_2d_2->Domain());

        auto projection_1 = ANurbs::PointOnCurveProjection<Kratos::array_1d<double, 3>>(
            curve_on_surface_3d_1, 0.1);
        auto projection_2 = ANurbs::PointOnCurveProjection<Kratos::array_1d<double, 3>>(
            curve_on_surface_3d_2, 0.1);

        auto curve_knot_intersections_1 = curve_on_surface_3d_1->Spans();

        std::vector<double> curve_knot_intersections_vector;
        for (int i = 0; i < curve_knot_intersections_1.size(); ++i)
        {
            curve_knot_intersections_vector.push_back(curve_knot_intersections_1[i].T0());
            auto point_3d = curve_on_surface_3d_2->PointAt(curve_knot_intersections_1[i].T0());
            projection_2.Compute(point_3d);
            curve_knot_intersections_vector.push_back(projection_2.Parameter());
        }
        curve_knot_intersections_vector.push_back(curve_on_surface_3d_1->Domain().T1());
        std::sort(curve_knot_intersections_vector.begin(), curve_knot_intersections_vector.end());

        std::cout << "here4: " << curve_knot_intersections_vector.size() << std::endl;
        for (int i = 0; i < curve_knot_intersections_vector.size() - 1; ++i)
        {
            std::cout << "curve_knot_intersections_vector[i]: " << curve_knot_intersections_vector[i] 
                << ", curve_knot_intersections_vector[i + 1]: " << curve_knot_intersections_vector[i + 1] << std::endl;
            if (std::abs(curve_knot_intersections_vector[i] - curve_knot_intersections_vector[i + 1]) > 1e-6)
            {
                ANurbs::Interval<double> interval(curve_knot_intersections_vector[i], curve_knot_intersections_vector[i + 1]);

                auto integration_points = ANurbs::IntegrationPoints<double>::Points1(number_of_points_per_knot_span, interval);

                std::cout << "here5" << std::endl;
                for (int ip = 0; ip < integration_points.size(); ++ip)
                {
                    auto point_3d = curve_on_surface_3d_2->PointAt(integration_points[ip].t);
                    auto point_3d2 = curve_on_surface_3d_1->PointAt(integration_points[ip].t);
                    auto derivatives = curve_2d_2->DerivativesAt(integration_points[ip].t, 1);
                    Element::GeometryType::PointsArrayType control_points;
                    Vector shape_function;
                    Matrix shape_function_derivative;
                    Matrix shape_function_second_derivative;

                    KRATOS_WATCH(point_3d2)
                    KRATOS_WATCH(point_3d)
                    std::cout << "here6: x:" << derivatives[0][0] << ", y:" << derivatives[0][1] << std::endl;
                    face_2.EvaluatePoint(derivatives[0][0], derivatives[0][1], control_points,
                        shape_function, shape_function_derivative, shape_function_second_derivative);

                    //SLAVE
                    projection_1.Compute(point_3d);

                    Vector shape_function_slave;
                    Matrix shape_function_derivative_slave;
                    Matrix shape_function_second_derivative_slave;

                    std::cout << "here7" << std::endl;
                    auto point_2d_slave = curve_2d_1->DerivativesAt(projection_1.Parameter(), 1);
                    face_1.EvaluatePoint(point_2d_slave[0][0], point_2d_slave[0][1], control_points,
                        shape_function_slave, shape_function_derivative_slave, shape_function_second_derivative_slave);

                    int id = 0;
                    if (rModelPart.GetRootModelPart().Elements().size() > 0)
                        id = rModelPart.GetRootModelPart().Elements().back().Id() + 1;

                    auto element = rModelPart.CreateNewElement(rName, id, control_points, this_property);

                    std::cout << "here8" << std::endl;
                    element->SetValue(SHAPE_FUNCTION_VALUES, shape_function);
                    element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, shape_function_derivative);
                    element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, shape_function_second_derivative);

                    Vector tangents(2);
                    tangents[0] = derivatives[1][0];
                    tangents[1] = derivatives[1][1];
                    element->SetValue(TANGENTS, tangents);

                    element->SetValue(INTEGRATION_WEIGHT, integration_points[ip].weight);

                    element->SetValue(SHAPE_FUNCTION_VALUES_SLAVE, shape_function_slave);
                    element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE, shape_function_derivative_slave);
                    element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE, shape_function_second_derivative_slave);

                    Vector tangents_slave(2);
                    tangents_slave[0] = point_2d_slave[1][0];
                    tangents_slave[1] = point_2d_slave[1][1];
                    element->SetValue(TANGENTS_SLAVE, tangents_slave);
                }
            }
        }
    }

    const BrepFace& BrepModel::GetFace(const int& brep_id) const
    {
        for (int i = 0; i < m_brep_faces.size(); ++i)
        {
            if (m_brep_faces[i].GetId() == brep_id)
                return m_brep_faces[i];
        }
        KRATOS_ERROR << "brep_id not in list of brep faces." << std::endl;
    }


    // --------------------------------------------------------------------------
    const std::vector<BrepFace>& BrepModel::GetFaceVector() const
    {
        return m_brep_faces;
    }
    const std::vector<BrepEdge>& BrepModel::GetEdgeVector() const
    {
        return m_brep_edges;
    }
    const std::vector<BrepVertex>& BrepModel::GetVertexVector() const
    {
        return m_brep_vertices;
    }

    BrepModel::BrepModel(
        int& brep_id,
        double& model_tolerance,
        std::vector<BrepFace>& faces,
        std::vector<BrepEdge>& edges,
        std::vector<BrepVertex>& vertices)
        : m_model_tolerance(model_tolerance),
          m_brep_faces(faces),
          m_brep_edges(edges),
          m_brep_vertices(vertices),
          IndexedObject(brep_id),
          Flags()
    {

    }
}  // namespace Kratos.

