//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// Project includes
#include "includes/define.h"
#include "mpm_modeler.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities.h"

namespace Kratos
{

    ///@name Stages
    ///@{

    void MpmModeler::SetupModelPart() {
        KRATOS_ERROR_IF_NOT(mParameters.Has("model_part_name_background"))
            << "MpmModeler: Missing \"model_part_name_background\" section" << std::endl;
        KRATOS_ERROR_IF_NOT(mParameters.Has("background_geometry_name"))
            << "MpmModeler: Missing \"background_geometry_name\" section" << std::endl;

        ModelPart& model_part_background = mpModel->HasModelPart(mParameters["model_part_name_background"].GetString())
            ? mpModel->GetModelPart(mParameters["model_part_name_background"].GetString())
            : mpModel->CreateModelPart(mParameters["model_part_name_background"].GetString());

        auto p_background_geometry = model_part_background.pGetGeometry(mParameters["background_geometry_name"].GetString());
        NurbsSurfaceGeometryPointerType p_nurbs_surface = dynamic_pointer_cast<NurbsSurfaceGeometryType>(p_background_geometry);

        ModelPart& model_part_analysis = mpModel->HasModelPart(mParameters["model_part_name_analysis"].GetString())
            ? mpModel->GetModelPart(mParameters["model_part_name_analysis"].GetString())
            : mpModel->CreateModelPart(mParameters["model_part_name_analysis"].GetString());

        for(auto& element : model_part_analysis.Elements())
        {
            auto p_geometry = element.pGetGeometry();
            array_1d<double, 3> r_global_location = p_geometry->Center();
            double weight = p_geometry->IntegrationPoints()[0].Weight() * p_geometry->DeterminantOfJacobian(0);
            array_1d<double, 3> local_location{ 0,0,0 };
            p_nurbs_surface->ProjectionPointGlobalToLocalSpace(r_global_location, local_location);
            IntegrationPoint<3> integration_point(local_location, weight);
            PointerVector<GeometryType> geometries(1);
            p_nurbs_surface->CreateQuadraturePointGeometries(
                geometries, 2, { integration_point }, p_geometry->GetDefaultIntegrationInfo());
            if (p_geometry->GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Curve_On_Surface_Geometry)
            {
                auto normal = p_geometry->Normal(0);
                KRATOS_WATCH(normal)
                element.SetGeometry(Kratos::make_shared<
                    QuadraturePointCurveOnSurfaceGeometry<NodeType, 2>>(
                        std::move(geometries(0)->Points()),
                        geometries(0)->GetGeometryData().GetGeometryShapeFunctionContainer(),
                        normal[0], normal[1],
                        p_nurbs_surface.get()));
            }
            else
            {
                element.SetGeometry(Kratos::make_shared<
                    QuadraturePointGeometry<NodeType, 2>>(
                        std::move(geometries(0)->Points()),
                        geometries(0)->GetGeometryData().GetGeometryShapeFunctionContainer(),
                        p_nurbs_surface.get()));
            }
            model_part_analysis.AddNodes(element.GetGeometry().begin(), element.GetGeometry().end());
        }
        for (auto& element : model_part_analysis.Conditions())
        {
            auto p_geometry = element.pGetGeometry();
            array_1d<double, 3> r_global_location = p_geometry->Center();
            double weight = p_geometry->IntegrationPoints()[0].Weight() * p_geometry->DeterminantOfJacobian(0);
            array_1d<double, 3> local_location{ 0,0,0 };
            p_nurbs_surface->ProjectionPointGlobalToLocalSpace(r_global_location, local_location);
            IntegrationPoint<3> integration_point(local_location, weight);
            PointerVector<GeometryType> geometries(1);
            p_nurbs_surface->CreateQuadraturePointGeometries(
                geometries, 2, { integration_point }, p_geometry->GetDefaultIntegrationInfo());
            KRATOS_WATCH("here")
            if (p_geometry->GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Curve_On_Surface_Geometry)
            {
                auto normal = p_geometry->Normal(0);
                auto detJ = p_geometry->DeterminantOfJacobian(0);
                KRATOS_WATCH(normal);
                KRATOS_WATCH(detJ);
                array_1d<double, 3> local_tangent;
                p_geometry->Calculate(LOCAL_TANGENT, local_tangent);
                KRATOS_WATCH(local_tangent)
                Matrix J;
                p_geometry->Jacobian(J, 0);
                auto t = prod(J, project(local_tangent, range(0, 2)));
                KRATOS_WATCH(J)
                KRATOS_WATCH(t)
                element.SetGeometry(Kratos::make_shared<
                    QuadraturePointCurveOnSurfaceGeometry<NodeType, 2>>(
                        std::move(geometries(0)->Points()),
                        geometries(0)->GetGeometryData().GetGeometryShapeFunctionContainer(),
                        normal[0], normal[1],
                        p_nurbs_surface.get()));
                Matrix J_inv;
                element.GetGeometry().InverseOfJacobian(J_inv, 0);
                Vector local_t = prod(J_inv, local_tangent);
                //std::swap(local_t[0], local_t[1]);
                element.GetGeometry().Assign(LOCAL_TANGENT, local_t);
                Matrix J_1;
                element.GetGeometry().Jacobian(J_1, 0);
                
                KRATOS_WATCH(element.GetGeometry().DeterminantOfJacobian(0))
                KRATOS_WATCH(J_inv)
                KRATOS_WATCH(J_1)
                KRATOS_WATCH(local_t)
                KRATOS_WATCH("new normal:");
                KRATOS_WATCH(element.GetGeometry().Normal(0));
            }
            else
            {
                element.SetGeometry(Kratos::make_shared<
                    QuadraturePointGeometry<NodeType, 2>>(
                        std::move(geometries(0)->Points()),
                        geometries(0)->GetGeometryData().GetGeometryShapeFunctionContainer(),
                        p_nurbs_surface.get()));
            }
            model_part_analysis.AddNodes(element.GetGeometry().begin(), element.GetGeometry().end());
        }
        if (mParameters.Has("supports"))
        {
            for (SizeType i = 0; i < mParameters["supports"].size(); ++i)
            {
                std::string model_part_name = mParameters["supports"][i]["model_part_name"].GetString();
                ModelPart& model_part_support = mpModel->HasModelPart(model_part_name)
                    ? mpModel->GetModelPart(model_part_name)
                    : mpModel->CreateModelPart(model_part_name);

                array_1d<int, 2> local_coordinates{ 0,0 };
                int support_index = mParameters["supports"][i]["support_index"].GetInt();
                if (support_index == 0) {
                    local_coordinates[0] = -1;
                }
                else if (support_index == 1) {
                    local_coordinates[0] = 1;
                    local_coordinates[1] = -1;
                }
                else if (support_index == 2) {
                    local_coordinates[0] = -1;
                    local_coordinates[1] = 1;
                }
                else if (support_index == 3) {
                    local_coordinates[1] = -1;
                }

                SizeType number_of_cps_u = p_nurbs_surface->PointsNumberInDirection(0);
                SizeType number_of_cps_v = p_nurbs_surface->PointsNumberInDirection(1);

                IndexType u_start = 0;
                IndexType u_end = number_of_cps_u;
                IndexType v_start = 0;
                IndexType v_end = number_of_cps_v;

                if (local_coordinates[0] >= 0) {
                    u_start = local_coordinates[0] * (number_of_cps_u - 1);
                    u_end = local_coordinates[0] * (number_of_cps_u - 1) + 1;
                }
                if (local_coordinates[1] >= 0) {
                    v_start = local_coordinates[1] * (number_of_cps_v - 1);
                    v_end = local_coordinates[1] * (number_of_cps_v - 1) + 1;
                }
                for (IndexType i = u_start; i < u_end; ++i) {
                    for (IndexType j = v_start; j < v_end; ++j) {
                        model_part_support.AddNode(p_nurbs_surface->pGetPoint(i + j * number_of_cps_u));
                    }
                }
            }
        }
    }

    ///@}
} // end namespace kratos
