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
    void BrepFace::GetGeometryIntegration(ModelPart& rModelPart,
        const std::string& rType,
        const std::string& rName,
        const int& rPropertiesId,
        const int& rShapeFunctionDerivativesOrder,
        std::vector<std::string> rVariables)
    {
        Properties::Pointer this_property = rModelPart.pGetProperties(rPropertiesId);

        auto spans_u = m_node_surface_geometry_3d->SpansU();
        auto spans_v = m_node_surface_geometry_3d->SpansV();

        ANurbs::SurfaceShapeEvaluator<double> shape(
            m_node_surface_geometry_3d->DegreeU(),
            m_node_surface_geometry_3d->DegreeV(),
            rShapeFunctionDerivativesOrder);

        for (int i = 0; i < spans_u.size(); ++i)
        {
            for (int j = 0; j < spans_v.size(); ++j)
            {
                ANurbs::Interval<double> domain_u(spans_u[i].T0(), spans_u[i].T1());
                ANurbs::Interval<double> domain_v(spans_v[j].T0(), spans_v[j].T1());


                std::cout << "domain_u: " << spans_u[i].T0() << ", " << spans_u[i].T1() << std::endl;
                std::cout << "domain_v: " << spans_v[j].T0() << ", " << spans_v[j].T1() << std::endl;

                auto integration_points = ANurbs::IntegrationPoints<double>::Points2(
                    m_node_surface_geometry_3d->DegreeU() + 1,
                    m_node_surface_geometry_3d->DegreeV() + 1,
                    domain_u,
                    domain_v);
                std::cout << "i: " << i << ", j" << j << std::endl;
                for (int k = 0; k < integration_points.size(); ++k)
                {
                    shape.Compute(
                        m_node_surface_geometry_3d->KnotsU(),
                        m_node_surface_geometry_3d->KnotsV(),
                        m_node_surface_geometry_3d->Weights(),
                        integration_points[k].u, 
                        integration_points[k].v);

                    std::cout << "u: " << integration_points[k].u << std::endl;
                    std::cout << "v: " << integration_points[k].v << std::endl;

                    std::cout << "k: " << k << std::endl;
                    Element::GeometryType::PointsArrayType non_zero_control_points;

                    std::cout << "shape.NbNonzeroPoles(): " << shape.NbNonzeroPoles() << std::endl;
                    Vector N_0 = ZeroVector(shape.NbNonzeroPoles());
                    Matrix N_1 = ZeroMatrix(shape.NbNonzeroPoles(), 2);
                    Matrix N_2 = ZeroMatrix(shape.NbNonzeroPoles(), 3);

                    for (int n = shape.FirstNonzeroPoleU(); n < shape.LastNonzeroPoleU(); ++n)
                    {
                        for (int m = shape.FirstNonzeroPoleV(); m < shape.LastNonzeroPoleV(); ++m)
                        {


                            int index = (n - shape.FirstNonzeroPoleU()) * shape.NbNonzeroPolesU() + (m-shape.FirstNonzeroPoleV());

                            int indexU = n;// -shape.FirstNonzeroPoleU();
                            int indexV = m;// -shape.FirstNonzeroPoleV();

                            N_0[index] = shape(0, indexU, indexV);
                            N_1(index, 0) = shape(1, indexU, indexV);
                            N_1(index, 1) = shape(2, indexU, indexV);
                            N_2(index, 0) = shape(3, indexU, indexV);
                            N_2(index, 1) = shape(5, indexU, indexV);
                            N_2(index, 2) = shape(4, indexU, indexV);
                        }
                    }

                    KRATOS_WATCH(N_0)
                    KRATOS_WATCH(N_1)
                    KRATOS_WATCH(N_2)

                    if (rType == "element")
                    {
                        int id = 0;
                        if (rModelPart.GetRootModelPart().Elements().size() > 0)
                            id = rModelPart.GetRootModelPart().Elements().back().Id() + 1;

                        auto element = rModelPart.CreateNewElement(rName, id, non_zero_control_points, this_property);

                        element->SetValue(SHAPE_FUNCTION_VALUES, N_0);
                        element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, N_1);
                        element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, N_2);
                        element->SetValue(INTEGRATION_WEIGHT, integration_points[j].weight);
                    }
                }
            }
        }
    }

    ///Constructor
    BrepFace::BrepFace(
        int rBrepId,
        bool rIsTrimmed,
        bool rIsRational,
        std::vector<BrepBoundaryLoop>& rTrimmingLoops,
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
          m_control_points_ids(rControlPointIds),
          m_model_part(rModelPart),
          IndexedObject(rBrepId),
          Flags()
    {
        int number_of_nodes_u = rKnotVectorU.size() - rP - 1;
        int number_of_nodes_v = rKnotVectorV.size() - rQ - 1;

        m_node_surface_geometry_3d = Kratos::make_unique<NodeSurfaceGeometry3D>(
            rP, rQ, number_of_nodes_u, number_of_nodes_v);

        for (int i = 0; i < rKnotVectorU.size()-2; ++i)
        {
            m_node_surface_geometry_3d->SetKnotU(i, rKnotVectorU(i+1));
        }

        for (int i = 0; i < rKnotVectorV.size()-2; ++i)
        {
            m_node_surface_geometry_3d->SetKnotV(i, rKnotVectorV(i+1));
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

        KRATOS_WATCH(m_node_surface_geometry_3d->KnotsU())
        KRATOS_WATCH(m_node_surface_geometry_3d->KnotsV())
        KRATOS_WATCH(m_node_surface_geometry_3d->NbPoles())
        //KRATOS_WATCH(m_node_surface_geometry_3d->Weights())
    }
} // namespace Kratos.