//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes

// External includes

// Project includes
#include "cad_json_output.h"

namespace Kratos
{
    void CadJsonOutput::GetParameters(
        ModelPart& rModelPart, Parameters& rCadGeometry, IndexType EchoLevel) {

        Parameters tolerances_parameters;
        tolerances_parameters.AddDouble("model_tolerance", 0.01);

        rCadGeometry.AddValue("tolerances", tolerances_parameters);
        rCadGeometry.AddInt("version_number", 1);

        Parameters breps_parameters;
        breps_parameters.AddEmptyArray("faces");
        breps_parameters.AddEmptyArray("edges");
        breps_parameters.AddEmptyArray("vertices");
        for (auto geometry_itr = rModelPart.GeometriesBegin(); geometry_itr != rModelPart.GeometriesEnd(); ++geometry_itr) {
            if (geometry_itr->GetGeometryType() == GeometryData::Kratos_Brep_Surface) {
                GetBrepSurfaceParameters(geometry_itr, breps_parameters, EchoLevel);
            }
        }

        rCadGeometry.AddEmptyArray("breps");
        rCadGeometry["breps"].Append(breps_parameters);
    }

    void CadJsonOutput::GetBrepSurfaceParameters(
        const typename ModelPart::GeometryIterator rGeometryIterator, Parameters& rBrepsParameters, IndexType EchoLevel) {
        const auto& r_aux_geometry = *rGeometryIterator;
        const BrepSurfaceType r_brep_surface_geom = dynamic_cast<const BrepSurfaceType&>(r_aux_geometry);
        Parameters face_parameters;
        face_parameters.AddInt("brep_id", r_brep_surface_geom.Id());
        face_parameters.AddBool("swapped_surface_normal", false);

        auto r_nurbs_surface_geometry = r_brep_surface_geom.pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
        auto r_nurbs_surface_geom = dynamic_pointer_cast<NurbsSurfaceType>(r_nurbs_surface_geometry);

        Parameters surface_parameters;
        surface_parameters.AddBool("is_trimmed", r_brep_surface_geom.IsTrimmed());
        surface_parameters.AddBool("is_rational", r_nurbs_surface_geom->IsRational());
        surface_parameters.AddEmptyArray("degrees");
        surface_parameters["degrees"].Append((int)r_nurbs_surface_geom->PolynomialDegree(0));
        surface_parameters["degrees"].Append((int)r_nurbs_surface_geom->PolynomialDegree(1));
        surface_parameters.AddEmptyArray("knot_vectors");

        Vector knots_u = r_nurbs_surface_geom->KnotsU();
        Vector knots_u_full = ZeroVector(knots_u.size() + 2);
        knots_u_full[0] = knots_u[0];
        knots_u_full[knots_u_full.size() - 1] = knots_u[knots_u.size() - 1];
        for (IndexType i = 0; i < knots_u.size(); ++i)
        {
            knots_u_full[i + 1] = knots_u[i];
        }

        Vector knots_v = r_nurbs_surface_geom->KnotsV();
        Vector knots_v_full = ZeroVector(knots_v.size() + 2);
        knots_v_full[0] = knots_v[0];
        knots_v_full[knots_v_full.size() - 1] = knots_v[knots_v.size() - 1];
        for (IndexType i = 0; i < knots_v.size(); ++i)
        {
            knots_v_full[i + 1] = knots_v[i];
        }

        surface_parameters["knot_vectors"].Append(knots_u_full);
        surface_parameters["knot_vectors"].Append(knots_v_full);
        surface_parameters.AddEmptyArray("control_points");

        Vector r_weights = r_nurbs_surface_geom->Weights();
        if (r_weights.size() == 0) {
            r_weights.resize(r_nurbs_surface_geom->size());
            std::fill(r_weights.begin(), r_weights.end(), 1.0);
        }
        Vector control_point = ZeroVector(4);
        for (IndexType i = 0; i < r_nurbs_surface_geom->size(); ++i)
        {
            auto& r_node = r_nurbs_surface_geom->GetPoint(i);
            const int control_points_id = r_node.GetId();
            control_point[0] = r_node.X();
            control_point[1] = r_node.Y();
            control_point[2] = r_node.Z();
            control_point[3] = r_weights[i];

            Parameters control_points_para;
            control_points_para.AddEmptyArray("cp");
            Parameters control_points_para2;
            control_points_para["cp"].Append(control_points_id);
            control_points_para["cp"].Append(control_point);

            surface_parameters["control_points"].Append(control_points_para["cp"]);
        }
        face_parameters.AddValue("surface", surface_parameters);

        const auto& outer_loops = r_brep_surface_geom.GetOuterLoops();

        face_parameters.AddEmptyArray("boundary_loops");
        for (IndexType loop_i = 0; loop_i < outer_loops.size(); ++loop_i)
        {
            Parameters boundary_loop_parameters;
            boundary_loop_parameters.AddString("loop_type", "outer");
            GetBoundaryLoopParameters(
                outer_loops[loop_i], boundary_loop_parameters, EchoLevel);
            face_parameters["boundary_loops"].Append(boundary_loop_parameters);
        }

        const auto& inner_loops = r_brep_surface_geom.GetInnerLoops();

        face_parameters.AddEmptyArray("boundary_loops");
        for (IndexType loop_i = 0; loop_i < inner_loops.size(); ++loop_i)
        {
            Parameters boundary_loop_parameters;
            boundary_loop_parameters.AddString("loop_type", "inner");
            GetBoundaryLoopParameters(
                inner_loops[loop_i], boundary_loop_parameters, EchoLevel);
            face_parameters["boundary_loops"].Append(boundary_loop_parameters);
        }

        rBrepsParameters["faces"].Append(face_parameters);
    }

    void CadJsonOutput::GetBoundaryLoopParameters(
        const BrepCurveOnSurfaceArrayType& rCurveOnSurfaceArray, Parameters& rBoundaryLoop, IndexType EchoLevel) {
        rBoundaryLoop.AddEmptyArray("trimming_curves");
        for (IndexType i = 0; i < rCurveOnSurfaceArray.size(); ++i) {
            Parameters trimming_curves_parameters;

            trimming_curves_parameters.AddInt("trim_index", rCurveOnSurfaceArray[i]->Id());
            trimming_curves_parameters.AddBool("curve_direction", rCurveOnSurfaceArray[i]->HasSameCurveDirection());

            auto p_curve_on_surface = rCurveOnSurfaceArray[i]->pGetCurveOnSurface();
            auto p_curve = p_curve_on_surface->pGetCurve();

            Parameters parameter_curve_parameters;
            parameter_curve_parameters.AddBool("is_rational", p_curve->IsRational());
            parameter_curve_parameters.AddInt("degree", p_curve->PolynomialDegree(0));
            Vector knots = p_curve->Knots();
            Vector knots_full = ZeroVector(knots.size() + 2);
            knots_full[0] = knots[0];
            knots_full[knots_full.size() - 1] = knots[knots.size() - 1];
            for (IndexType j = 0; j < knots.size(); ++j)
            {
                knots_full[j + 1] = knots[j];
            }

            parameter_curve_parameters.AddVector("knot_vector", knots_full);
            auto nurbs_interval = rCurveOnSurfaceArray[i]->DomainInterval();
            Vector nurbs_interval_vector = Vector(2);
            nurbs_interval_vector[0] = nurbs_interval.GetT0();
            nurbs_interval_vector[1] = nurbs_interval.GetT1();
            parameter_curve_parameters.AddVector("active_range", nurbs_interval_vector);
            parameter_curve_parameters.AddEmptyArray("control_points");

            Vector r_weights = p_curve->Weights();
            if (r_weights.size() == 0) {
                r_weights.resize(p_curve->size());
                std::fill(r_weights.begin(), r_weights.end(), 1.0);
            }
            Vector control_points = ZeroVector(4);
            for (IndexType j = 0; j < p_curve->size(); ++j) {
                const auto& r_node = p_curve->GetPoint(j);
                control_points[0] = r_node[0];
                control_points[1] = r_node[1];
                control_points[2] = r_node[2];
                control_points[3] = r_weights[j];

                Parameters control_points_parameters;
                control_points_parameters.AddEmptyArray("cp");
                control_points_parameters["cp"].Append(0);
                control_points_parameters["cp"].Append(control_points);
                parameter_curve_parameters["control_points"].Append(control_points_parameters["cp"]);
            }
            trimming_curves_parameters.AddValue("parameter_curve", parameter_curve_parameters);

            rBoundaryLoop["trimming_curves"].Append(trimming_curves_parameters);
        }
    }
}  // namespace Kratos.
