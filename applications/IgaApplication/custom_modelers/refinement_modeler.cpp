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
#include "refinement_modeler.h"

namespace Kratos
{
    ///@name Stages
    ///@{

    void RefinementModeler::PrepareGeometryModel()
    {
        const std::string DataFileName = mParameters.Has("refinements_file_name")
            ? mParameters["refinements_file_name"].GetString()
            : "refinements.iga.json";

        KRATOS_INFO_IF("::[RefinementModeler]::", mEchoLevel > 0) << "Refining model by: " << DataFileName << std::endl;

        const Parameters refinements_parameters = ReadParamatersFile(DataFileName);

        ApplyRefinements(refinements_parameters);
    }

    void RefinementModeler::ApplyRefinements(
        const Parameters rParameters) const
    {
        KRATOS_ERROR_IF_NOT(rParameters.Has("refinements"))
            << "Parameters do not have refinements section.\n"
            << rParameters << std::endl;

        KRATOS_ERROR_IF_NOT(rParameters["refinements"].IsArray())
            << "refinements section need to be of type array.\n"
            << rParameters << std::endl;

        for (IndexType i = 0; i < rParameters["refinements"].size(); ++i) {
            ApplyRefinement(rParameters["refinements"][i]);
        }
    }

    void RefinementModeler::ApplyRefinement(
        const Parameters rParameters) const
    {
        KRATOS_ERROR_IF_NOT(rParameters.Has("model_part_name"))
            << "Missing \"model_part_name\" in refinements block.\n"
            << rParameters << std::endl;
        const std::string model_part_name = rParameters["model_part_name"].GetString();
        ModelPart& r_model_part = mpModel->HasModelPart(model_part_name)
            ? mpModel->GetModelPart(model_part_name)
            : mpModel->CreateModelPart(model_part_name);

        // Generate the list of geometries, which are needed, here.
        GeometriesArrayType geometry_list;
        GetGeometryList(geometry_list, r_model_part, rParameters);

        KRATOS_ERROR_IF_NOT(rParameters.Has("parameters"))
            << "Missing \"parameters\" in refinements block.\n"
            << rParameters << std::endl;
        KRATOS_ERROR_IF_NOT(rParameters.Has("geometry_type"))
            << "Missing \"geometry_type\".\n"
            << rParameters << std::endl;
        if (rParameters["geometry_type"].GetString() == "NurbsSurface") {
            for (IndexType n = 0; n < geometry_list.size(); ++n) {
                auto p_nurbs_surface = (geometry_list[n].size() > 0)
                    ? dynamic_pointer_cast<NurbsSurfaceGeometry<3, PointerVector<Node>>>(geometry_list(n))
                    : dynamic_pointer_cast<NurbsSurfaceGeometry<3, PointerVector<Node>>>(geometry_list(n)->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX));

                KRATOS_INFO_IF("[Refinment Modeler]", mEchoLevel > 3 ) << "box_refinement_u_v_min_max" << std::endl;
                // MODIFIED!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                Vector box_refinement_u_v_min_max;
                double delta_u_refined;
                double delta_v_refined; 
                if (rParameters["local_refinement_parameters"]["do_you_want_to_refine_inner"].GetBool()) {
                    box_refinement_u_v_min_max = r_model_part.GetProcessInfo().GetValue(MARKER_MESHES);

                    int nb_per_span_u = rParameters["parameters"]["insert_nb_per_span_u"].GetInt();
                    int nb_per_span_v = rParameters["parameters"]["insert_nb_per_span_v"].GetInt();
                    const int nb_per_span_u_refined = (nb_per_span_u+1) * rParameters["local_refinement_parameters"]["how_many_times_refine_inner"].GetInt()-1;
                    const int nb_per_span_v_refined = (nb_per_span_v+1) * rParameters["local_refinement_parameters"]["how_many_times_refine_inner"].GetInt()-1;

                    std::vector<double> spans_local_space_u;
                    p_nurbs_surface->SpansLocalSpace(spans_local_space_u, 0);
                    std::vector<double> spans_local_space_v;
                    p_nurbs_surface->SpansLocalSpace(spans_local_space_v, 1);
                    
                    delta_u_refined = (spans_local_space_u[1] - spans_local_space_u[0]) / (nb_per_span_u_refined + 1);
                    delta_v_refined = (spans_local_space_v[1] - spans_local_space_v[0]) / (nb_per_span_v_refined + 1);
                }
                KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Refining nurbs surface by elevating degree in u: " << std::endl;
       
                if (rParameters["parameters"].Has("increase_degree_u")) {
                    SizeType increase_degree_u = rParameters["parameters"]["increase_degree_u"].GetInt();

                    if (increase_degree_u > 0) {
                        std::vector<double> spans_local_space;
                        p_nurbs_surface->SpansLocalSpace(spans_local_space, 0);

                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Refining nurbs surface #" << p_nurbs_surface->Id() << " by elevating degree in u: "
                            << increase_degree_u << std::endl;

                        PointerVector<NodeType> PointsRefined;
                        Vector KnotsURefined;
                        Vector WeightsRefined;

                        NurbsSurfaceRefinementUtilities::DegreeElevationU(*(p_nurbs_surface.get()), increase_degree_u,
                            PointsRefined, KnotsURefined, WeightsRefined);

                        // Recreate nodes in model part to ensure correct assignment of dofs
                        IndexType node_id = (r_model_part.GetParentModelPart().NodesEnd() - 1)->Id() + 1;
                        for (IndexType i = 0; i < PointsRefined.size(); ++i) {
                            if (PointsRefined(i)->Id() == 0) {
                                PointsRefined(i) = r_model_part.CreateNewNode(node_id, PointsRefined[i][0], PointsRefined[i][1], PointsRefined[i][2]);
                                node_id++;
                            }
                        }

                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "New knot vector: " << KnotsURefined
                            << ", new weights vector: " << WeightsRefined << std::endl;

                        p_nurbs_surface->SetInternals(PointsRefined,
                            p_nurbs_surface->PolynomialDegreeU() + increase_degree_u, p_nurbs_surface->PolynomialDegreeV(),
                            KnotsURefined, p_nurbs_surface->KnotsV(),
                            WeightsRefined);
                    } else {
                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Trying to refine nurbs surface #" << p_nurbs_surface->Id()
                            << " by elevating degree in u, however increase_degree_u is set to 0." << std::endl;
                    }
                }
                KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Refining nurbs surface by elevating degree in v: " << std::endl;
                            
                if (rParameters["parameters"].Has("increase_degree_v")) {
                    SizeType increase_degree_v = rParameters["parameters"]["increase_degree_v"].GetInt();

                    if (increase_degree_v > 0) {
                        std::vector<double> spans_local_space;
                        p_nurbs_surface->SpansLocalSpace(spans_local_space, 0);

                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Refining nurbs surface #" << p_nurbs_surface->Id() << " by elevating degree in v: "
                            << increase_degree_v << std::endl;

                        PointerVector<NodeType> PointsRefined;
                        Vector KnotsVRefined;
                        Vector WeightsRefined;

                        NurbsSurfaceRefinementUtilities::DegreeElevationV(*(p_nurbs_surface.get()), increase_degree_v,
                            PointsRefined, KnotsVRefined, WeightsRefined);

                        // Recreate nodes in model part to ensure correct assignment of dofs
                        IndexType node_id = (r_model_part.GetParentModelPart().NodesEnd() - 1)->Id() + 1;
                        for (IndexType i = 0; i < PointsRefined.size(); ++i) {
                            if (PointsRefined(i)->Id() == 0) {
                                PointsRefined(i) = r_model_part.CreateNewNode(node_id, PointsRefined[i][0], PointsRefined[i][1], PointsRefined[i][2]);
                                node_id++;
                            }
                        }
                    
                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "New knot vector: " << KnotsVRefined
                            << ", new weights vector: " << WeightsRefined << std::endl;

                        p_nurbs_surface->SetInternals(PointsRefined,
                            p_nurbs_surface->PolynomialDegreeU(), p_nurbs_surface->PolynomialDegreeV() + increase_degree_v,
                            p_nurbs_surface->KnotsU(), KnotsVRefined,
                            WeightsRefined);
                    } else {
                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Trying to refine nurbs surface #" << p_nurbs_surface->Id()
                            << " by elevating degree in v, however increase_degree_v is set to 0." << std::endl;
                    }
                }

                KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "starting u knot insertion: " << std::endl;
                if (rParameters["parameters"].Has("insert_nb_per_span_u")) {
                    const IndexType nb_per_span_u = rParameters["parameters"]["insert_nb_per_span_u"].GetInt();

                    if (nb_per_span_u > 0) {
                        std::vector<double> spans_local_space;
                        p_nurbs_surface->SpansLocalSpace(spans_local_space, 0);

                        std::vector<double> knots_to_insert_u;

                        if (rParameters["local_refinement_parameters"]["do_you_want_to_refine_inner"].GetBool()) {
                            for (IndexType i = 0; i < spans_local_space.size() - 1; ++i) {
                                const double delta_u = (spans_local_space[i + 1] - spans_local_space[i]) / (nb_per_span_u + 1);
                                for (IndexType j = 1; j < nb_per_span_u + 1; ++j) {
                                    double u_local = spans_local_space[i] + delta_u*j;
                                    if ( (u_local > box_refinement_u_v_min_max[0] 
                                            || std::abs(u_local - box_refinement_u_v_min_max[0]) < 1e-13 )
                                            && u_local < box_refinement_u_v_min_max[1]) {
                                        u_local = box_refinement_u_v_min_max[0];
                                        while (u_local < box_refinement_u_v_min_max[1]) {
                                            knots_to_insert_u.push_back(u_local);
                                            u_local += delta_u_refined;
                                        }
                                        if (std::abs(knots_to_insert_u[knots_to_insert_u.size()-1] - box_refinement_u_v_min_max[1]) > 1e-13){
                                            knots_to_insert_u.push_back(box_refinement_u_v_min_max[1]);
                                        }
                                        j = ((u_local - spans_local_space[i]) / delta_u) +1; 
                                        // if (std::abs(spans_local_space[i] + delta_u*j - u_local) < 1e-13) { j++;}
                                    }
                                    knots_to_insert_u.push_back(spans_local_space[i] + delta_u * j);
                                }
                            }
                        } else {
                            for (IndexType i = 0; i < spans_local_space.size() - 1; ++i) {
                                const double delta_u = (spans_local_space[i + 1] - spans_local_space[i]) / (nb_per_span_u + 1);
                                for (IndexType j = 1; j < nb_per_span_u + 1; ++j) {
                                    knots_to_insert_u.push_back(spans_local_space[i] + delta_u * j);
                                }
                            }
                        }
                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Refining nurbs surface #" << p_nurbs_surface->Id() << " by inserting knots in u: "
                            << knots_to_insert_u << std::endl;

                        // KRATOS_WATCH(knots_to_insert_u)
                        PointerVector<NodeType> PointsRefined;
                        Vector KnotsURefined;
                        Vector WeightsRefined;

                        NurbsSurfaceRefinementUtilities::KnotRefinementU(*(p_nurbs_surface.get()), knots_to_insert_u,
                            PointsRefined, KnotsURefined, WeightsRefined);

                        // Recreate nodes in model part to ensure correct assignment of dofs
                        IndexType node_id = (r_model_part.GetParentModelPart().NodesEnd() - 1)->Id() + 1;
                        for (IndexType i = 0; i < PointsRefined.size(); ++i) {
                            if (PointsRefined(i)->Id() == 0) {
                                PointsRefined(i) = r_model_part.CreateNewNode(node_id, PointsRefined[i][0], PointsRefined[i][1], PointsRefined[i][2]);
                                node_id++;
                            }
                        }

                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "New knot vector: " << KnotsURefined
                            << ", new weights vector: " << WeightsRefined << std::endl;

                        p_nurbs_surface->SetInternals(PointsRefined,
                            p_nurbs_surface->PolynomialDegreeU(), p_nurbs_surface->PolynomialDegreeV(),
                            KnotsURefined, p_nurbs_surface->KnotsV(),
                            WeightsRefined);
                    } else {
                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Trying to refine nurbs surface #" << p_nurbs_surface->Id()
                            << " by inserting knots in u, however insert_nb_per_span_u is set to 0." << std::endl;
                    }
                }
                KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "starting v knot insertion: " << std::endl;
                if (rParameters["parameters"].Has("insert_nb_per_span_v")) {
                    const IndexType nb_per_span_v = rParameters["parameters"]["insert_nb_per_span_v"].GetInt();

                    if (nb_per_span_v > 0) {
                    std::vector<double> spans_local_space;
                    p_nurbs_surface->SpansLocalSpace(spans_local_space, 1);

                    std::vector<double> knots_to_insert_v;

                    if (rParameters["local_refinement_parameters"]["do_you_want_to_refine_inner"].GetBool()) {
                        for (IndexType i = 0; i < spans_local_space.size() - 1; ++i) {
                            const double delta_v = (spans_local_space[i + 1] - spans_local_space[i]) / (nb_per_span_v + 1);
                            for (IndexType j = 1; j < nb_per_span_v + 1; ++j) {
                                double v_local = spans_local_space[i] + delta_v*j;
                                if ( (v_local > box_refinement_u_v_min_max[2] 
                                        || std::abs(v_local - box_refinement_u_v_min_max[2]) < 1e-13 )
                                        && v_local < box_refinement_u_v_min_max[3]) {
                                    v_local = box_refinement_u_v_min_max[2];
                                    while (v_local < box_refinement_u_v_min_max[3]) {
                                        knots_to_insert_v.push_back(v_local);
                                        v_local += delta_v_refined;
                                    }
                                    if (std::abs(knots_to_insert_v[knots_to_insert_v.size()-1] - box_refinement_u_v_min_max[3]) > 1e-13){
                                        knots_to_insert_v.push_back(box_refinement_u_v_min_max[3]);
                                    }
                                    j = ((v_local - spans_local_space[i]) / delta_v) +1; 
                                    // if (std::abs(spans_local_space[i] + delta_v*j - v_local) < 1e-13) { j++;}
                                }
                                knots_to_insert_v.push_back(spans_local_space[i] + delta_v * j);
                            }
                        }
                    } else {
                        for (IndexType i = 0; i < spans_local_space.size() - 1; ++i) {
                            const double delta_v = (spans_local_space[i + 1] - spans_local_space[i]) / (nb_per_span_v + 1);
                            for (IndexType j = 1; j < nb_per_span_v + 1; ++j) {
                                knots_to_insert_v.push_back(spans_local_space[i] + delta_v * j);
                            }
                        }
                    }
                    // KRATOS_WATCH(knots_to_insert_v)

                    KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                        << "Refining nurbs surface #" << p_nurbs_surface->Id() << " by inserting knots in v: "
                        << knots_to_insert_v << std::endl;

                    PointerVector<NodeType> PointsRefined;
                    Vector KnotsVRefined;
                    Vector WeightsRefined;

                    NurbsSurfaceRefinementUtilities::KnotRefinementV(*(p_nurbs_surface.get()), knots_to_insert_v,
                        PointsRefined, KnotsVRefined, WeightsRefined);

                    // Recreate nodes in model part to ensure correct assignment of dofs
                    IndexType node_id = (r_model_part.GetParentModelPart().NodesEnd() - 1)->Id() + 1;
                    for (IndexType i = 0; i < PointsRefined.size(); ++i) {
                        if (PointsRefined(i)->Id() == 0) {
                            PointsRefined(i) = r_model_part.CreateNewNode(node_id, PointsRefined[i][0], PointsRefined[i][1], PointsRefined[i][2]);
                            node_id++;
                        }
                    }

                    KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                        << "New knot vector: " << KnotsVRefined
                        << ", new weights vector: " << WeightsRefined << std::endl;

                    p_nurbs_surface->SetInternals(PointsRefined,
                        p_nurbs_surface->PolynomialDegreeU(), p_nurbs_surface->PolynomialDegreeV(),
                        p_nurbs_surface->KnotsU(), KnotsVRefined,
                        WeightsRefined);
                    }
                    else {
                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Trying to refine nurbs surface #" << p_nurbs_surface->Id()
                            << " by inserting knots in v, however insert_nb_per_span_v is set to 0." << std::endl;
                    }
                }
                array_1d<double, 3> rOutput{0.0, 0.0, 0.0};  // Initialize the output array
                p_nurbs_surface->Calculate(CHARACTERISTIC_GEOMETRY_LENGTH, rOutput);
                KRATOS_WATCH(rOutput)
            }
        } 
        // ------------------------------------------------------------------------------------------
        // ##########################################################################################
        // 3D REFINEMENT -> NurbsVolume
        // ##########################################################################################
        // ------------------------------------------------------------------------------------------
        else if (rParameters["geometry_type"].GetString() == "NurbsVolume") {
            KRATOS_INFO_IF("[Refinment Modeler]", mEchoLevel > 3 ) << "Refinement of Nurbs Volume" << std::endl;
            for (IndexType n = 0; n < geometry_list.size(); ++n) {
                auto p_nurbs_volume = (geometry_list[n].size() > 0)
                    ? dynamic_pointer_cast<NurbsVolumeGeometry<PointerVector<Node>>>(geometry_list(n))
                    : dynamic_pointer_cast<NurbsVolumeGeometry<PointerVector<Node>>>(geometry_list(n)->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX));

                KRATOS_INFO_IF("[Refinment Modeler]", mEchoLevel > 3 ) << "box_refinement_u_v_min_max" << std::endl;
                // // MODIFIED!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                Vector box_refinement_u_v_min_max;
                double delta_u_refined;
                double delta_v_refined; 
                // if (rParameters["local_refinement_parameters"]["do_you_want_to_refine_inner"].GetBool()) {
                //     box_refinement_u_v_min_max = r_model_part.GetProcessInfo().GetValue(MARKER_MESHES);

                //     int nb_per_span_u = rParameters["parameters"]["insert_nb_per_span_u"].GetInt();
                //     int nb_per_span_v = rParameters["parameters"]["insert_nb_per_span_v"].GetInt();
                //     const int nb_per_span_u_refined = (nb_per_span_u+1) * rParameters["local_refinement_parameters"]["how_many_times_refine_inner"].GetInt()-1;
                //     const int nb_per_span_v_refined = (nb_per_span_v+1) * rParameters["local_refinement_parameters"]["how_many_times_refine_inner"].GetInt()-1;

                //     std::vector<double> spans_local_space_u;
                //     p_nurbs_surface->SpansLocalSpace(spans_local_space_u, 0);
                //     std::vector<double> spans_local_space_v;
                //     p_nurbs_surface->SpansLocalSpace(spans_local_space_v, 1);
                    
                //     delta_u_refined = (spans_local_space_u[1] - spans_local_space_u[0]) / (nb_per_span_u_refined + 1);
                //     delta_v_refined = (spans_local_space_v[1] - spans_local_space_v[0]) / (nb_per_span_v_refined + 1);
                // }
                KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Refining nurbs surface by elevating degree in u: " << std::endl;
                            
                if (rParameters["parameters"].Has("increase_degree_u")) {
                    SizeType increase_degree_u = rParameters["parameters"]["increase_degree_u"].GetInt();

                    if (increase_degree_u > 0) {
                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Refining nurbs volume #" << p_nurbs_volume->Id() << " by elevating degree in u: "
                            << increase_degree_u << std::endl;

                        PointerVector<NodeType> PointsRefined;
                        Vector KnotsURefined;
                        Vector WeightsRefined;

                        NurbsVolumeRefinementUtilities::DegreeElevationU(*(p_nurbs_volume.get()), increase_degree_u,
                            PointsRefined, KnotsURefined, WeightsRefined);

                        // Recreate nodes in model part to ensure correct assignment of dofs
                        IndexType node_id = (r_model_part.GetParentModelPart().NodesEnd() - 1)->Id() + 1;
                        for (IndexType i = 0; i < PointsRefined.size(); ++i) {
                            if (PointsRefined(i)->Id() == 0) {
                                PointsRefined(i) = r_model_part.CreateNewNode(node_id, PointsRefined[i][0], PointsRefined[i][1], PointsRefined[i][2]);
                                node_id++;
                            }
                        }

                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "New knot vector: " << KnotsURefined
                            << ", new weights vector: " << WeightsRefined << std::endl;

                        p_nurbs_volume->SetInternals(PointsRefined,
                            p_nurbs_volume->PolynomialDegreeU() + increase_degree_u, p_nurbs_volume->PolynomialDegreeV(), p_nurbs_volume->PolynomialDegreeW(),
                            KnotsURefined, p_nurbs_volume->KnotsV(), p_nurbs_volume->KnotsW());
                    } else {
                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Trying to refine nurbs surface #" << p_nurbs_volume->Id()
                            << " by elevating degree in u, however increase_degree_u is set to 0." << std::endl;
                    }
                }
                KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Refining nurbs surface by elevating degree in v: " << std::endl;
                            
                if (rParameters["parameters"].Has("increase_degree_v")) {
                    SizeType increase_degree_v = rParameters["parameters"]["increase_degree_v"].GetInt();

                    if (increase_degree_v > 0) {

                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Refining nurbs surface #" << p_nurbs_volume->Id() << " by elevating degree in v: "
                            << increase_degree_v << std::endl;

                        PointerVector<NodeType> PointsRefined;
                        Vector KnotsVRefined;
                        Vector WeightsRefined;

                        NurbsVolumeRefinementUtilities::DegreeElevationV(*(p_nurbs_volume.get()), increase_degree_v,
                            PointsRefined, KnotsVRefined, WeightsRefined);

                        // Recreate nodes in model part to ensure correct assignment of dofs
                        IndexType node_id = (r_model_part.GetParentModelPart().NodesEnd() - 1)->Id() + 1;
                        for (IndexType i = 0; i < PointsRefined.size(); ++i) {
                            if (PointsRefined(i)->Id() == 0) {
                                PointsRefined(i) = r_model_part.CreateNewNode(node_id, PointsRefined[i][0], PointsRefined[i][1], PointsRefined[i][2]);
                                node_id++;
                            }
                        }
                    
                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "New knot vector: " << KnotsVRefined
                            << ", new weights vector: " << WeightsRefined << std::endl;

                        p_nurbs_volume->SetInternals(PointsRefined,
                            p_nurbs_volume->PolynomialDegreeU(), p_nurbs_volume->PolynomialDegreeV() + increase_degree_v, p_nurbs_volume->PolynomialDegreeW(),
                            p_nurbs_volume->KnotsU(), KnotsVRefined, p_nurbs_volume->KnotsW());
                    } else {
                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Trying to refine nurbs surface #" << p_nurbs_volume->Id()
                            << " by elevating degree in v, however increase_degree_v is set to 0." << std::endl;
                    }
                }

                if (rParameters["parameters"].Has("increase_degree_w")) {
                    SizeType increase_degree_w = rParameters["parameters"]["increase_degree_w"].GetInt();

                    if (increase_degree_w > 0) {

                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Refining nurbs surface #" << p_nurbs_volume->Id() << " by elevating degree in v: "
                            << increase_degree_w << std::endl;

                        PointerVector<NodeType> PointsRefined;
                        Vector KnotsWRefined;
                        Vector WeightsRefined;

                        NurbsVolumeRefinementUtilities::DegreeElevationW(*(p_nurbs_volume.get()), increase_degree_w,
                            PointsRefined, KnotsWRefined, WeightsRefined);

                        // Recreate nodes in model part to ensure correct assignment of dofs
                        IndexType node_id = (r_model_part.GetParentModelPart().NodesEnd() - 1)->Id() + 1;
                        for (IndexType i = 0; i < PointsRefined.size(); ++i) {
                            if (PointsRefined(i)->Id() == 0) {
                                PointsRefined(i) = r_model_part.CreateNewNode(node_id, PointsRefined[i][0], PointsRefined[i][1], PointsRefined[i][2]);
                                node_id++;
                            }
                        }
                    
                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "New knot vector: " << KnotsWRefined
                            << ", new weights vector: " << WeightsRefined << std::endl;

                        p_nurbs_volume->SetInternals(PointsRefined,
                            p_nurbs_volume->PolynomialDegreeU(), p_nurbs_volume->PolynomialDegreeV(), p_nurbs_volume->PolynomialDegreeW() + increase_degree_w,
                            p_nurbs_volume->KnotsU(), p_nurbs_volume->KnotsV(), KnotsWRefined);
                    } else {
                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Trying to refine nurbs surface #" << p_nurbs_volume->Id()
                            << " by elevating degree in v, however increase_degree_v is set to 0." << std::endl;
                    }
                }
                // -----------------------------------------------------------------------------------------------------------------

                KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "starting u knot insertion: " << std::endl;
                if (rParameters["parameters"].Has("insert_nb_per_span_u")) {
                    const IndexType nb_per_span_u = rParameters["parameters"]["insert_nb_per_span_u"].GetInt();

                    if (nb_per_span_u > 0) {
                        std::vector<double> spans_local_space;
                        p_nurbs_volume->SpansLocalSpace(spans_local_space, 0);

                        std::vector<double> knots_to_insert_u;

                        if (rParameters["local_refinement_parameters"]["do_you_want_to_refine_inner"].GetBool()) {
                            for (IndexType i = 0; i < spans_local_space.size() - 1; ++i) {
                                const double delta_u = (spans_local_space[i + 1] - spans_local_space[i]) / (nb_per_span_u + 1);
                                for (IndexType j = 1; j < nb_per_span_u + 1; ++j) {
                                    double u_local = spans_local_space[i] + delta_u*j;
                                    if ( (u_local > box_refinement_u_v_min_max[0] 
                                            || std::abs(u_local - box_refinement_u_v_min_max[0]) < 1e-13 )
                                            && u_local < box_refinement_u_v_min_max[1]) {
                                        u_local = box_refinement_u_v_min_max[0];
                                        while (u_local < box_refinement_u_v_min_max[1]) {
                                            knots_to_insert_u.push_back(u_local);
                                            u_local += delta_u_refined;
                                        }
                                        if (std::abs(knots_to_insert_u[knots_to_insert_u.size()-1] - box_refinement_u_v_min_max[1]) > 1e-13){
                                            knots_to_insert_u.push_back(box_refinement_u_v_min_max[1]);
                                        }
                                        j = ((u_local - spans_local_space[i]) / delta_u) +1; 
                                        // if (std::abs(spans_local_space[i] + delta_u*j - u_local) < 1e-13) { j++;}
                                    }
                                    knots_to_insert_u.push_back(spans_local_space[i] + delta_u * j);
                                }
                            }
                        } else {
                            for (IndexType i = 0; i < spans_local_space.size() - 1; ++i) {
                                const double delta_u = (spans_local_space[i + 1] - spans_local_space[i]) / (nb_per_span_u + 1);
                                for (IndexType j = 1; j < nb_per_span_u + 1; ++j) {
                                    knots_to_insert_u.push_back(spans_local_space[i] + delta_u * j);
                                }
                            }
                        }
                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Refining nurbs surface #" << p_nurbs_volume->Id() << " by inserting knots in u: "
                            << knots_to_insert_u << std::endl;

                        // KRATOS_WATCH(knots_to_insert_u)
                        PointerVector<NodeType> PointsRefined;
                        Vector KnotsURefined;
                        Vector WeightsRefined;

                        NurbsVolumeRefinementUtilities::KnotRefinementU(*(p_nurbs_volume.get()), knots_to_insert_u,
                            PointsRefined, KnotsURefined);

                        // Recreate nodes in model part to ensure correct assignment of dofs
                        IndexType node_id = (r_model_part.GetParentModelPart().NodesEnd() - 1)->Id() + 1;
                        for (IndexType i = 0; i < PointsRefined.size(); ++i) {
                            if (PointsRefined(i)->Id() == 0) {
                                PointsRefined(i) = r_model_part.CreateNewNode(node_id, PointsRefined[i][0], PointsRefined[i][1], PointsRefined[i][2]);
                                node_id++;
                            }
                        }

                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "New knot vector: " << KnotsURefined
                            << ", new weights vector: " << WeightsRefined << std::endl;

                        p_nurbs_volume->SetInternals(PointsRefined,
                            p_nurbs_volume->PolynomialDegreeU(), p_nurbs_volume->PolynomialDegreeV(), p_nurbs_volume->PolynomialDegreeW(),
                            KnotsURefined, p_nurbs_volume->KnotsV(), p_nurbs_volume->KnotsW());
                    } else {
                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Trying to refine nurbs surface #" << p_nurbs_volume->Id()
                            << " by inserting knots in u, however insert_nb_per_span_u is set to 0." << std::endl;
                    }
                }
                KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "starting v knot insertion: " << std::endl;
                if (rParameters["parameters"].Has("insert_nb_per_span_v")) {
                    const IndexType nb_per_span_v = rParameters["parameters"]["insert_nb_per_span_v"].GetInt();

                    if (nb_per_span_v > 0) {
                    std::vector<double> spans_local_space;
                    p_nurbs_volume->SpansLocalSpace(spans_local_space, 1);

                    std::vector<double> knots_to_insert_v;

                    if (rParameters["local_refinement_parameters"]["do_you_want_to_refine_inner"].GetBool()) {
                        for (IndexType i = 0; i < spans_local_space.size() - 1; ++i) {
                            const double delta_v = (spans_local_space[i + 1] - spans_local_space[i]) / (nb_per_span_v + 1);
                            for (IndexType j = 1; j < nb_per_span_v + 1; ++j) {
                                double v_local = spans_local_space[i] + delta_v*j;
                                if ( (v_local > box_refinement_u_v_min_max[2] 
                                        || std::abs(v_local - box_refinement_u_v_min_max[2]) < 1e-13 )
                                        && v_local < box_refinement_u_v_min_max[3]) {
                                    v_local = box_refinement_u_v_min_max[2];
                                    while (v_local < box_refinement_u_v_min_max[3]) {
                                        knots_to_insert_v.push_back(v_local);
                                        v_local += delta_v_refined;
                                    }
                                    if (std::abs(knots_to_insert_v[knots_to_insert_v.size()-1] - box_refinement_u_v_min_max[3]) > 1e-13){
                                        knots_to_insert_v.push_back(box_refinement_u_v_min_max[3]);
                                    }
                                    j = ((v_local - spans_local_space[i]) / delta_v) +1; 
                                    // if (std::abs(spans_local_space[i] + delta_v*j - v_local) < 1e-13) { j++;}
                                }
                                knots_to_insert_v.push_back(spans_local_space[i] + delta_v * j);
                            }
                        }
                    } else {
                        for (IndexType i = 0; i < spans_local_space.size() - 1; ++i) {
                            const double delta_v = (spans_local_space[i + 1] - spans_local_space[i]) / (nb_per_span_v + 1);
                            for (IndexType j = 1; j < nb_per_span_v + 1; ++j) {
                                knots_to_insert_v.push_back(spans_local_space[i] + delta_v * j);
                            }
                        }
                    }
                    // KRATOS_WATCH(knots_to_insert_v)

                    KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                        << "Refining nurbs surface #" << p_nurbs_volume->Id() << " by inserting knots in v: "
                        << knots_to_insert_v << std::endl;

                    PointerVector<NodeType> PointsRefined;
                    Vector KnotsVRefined;
                    Vector WeightsRefined;

                    NurbsVolumeRefinementUtilities::KnotRefinementV(*(p_nurbs_volume.get()), knots_to_insert_v,
                        PointsRefined, KnotsVRefined);

                    // Recreate nodes in model part to ensure correct assignment of dofs
                    IndexType node_id = (r_model_part.GetParentModelPart().NodesEnd() - 1)->Id() + 1;
                    for (IndexType i = 0; i < PointsRefined.size(); ++i) {
                        if (PointsRefined(i)->Id() == 0) {
                            PointsRefined(i) = r_model_part.CreateNewNode(node_id, PointsRefined[i][0], PointsRefined[i][1], PointsRefined[i][2]);
                            node_id++;
                        }
                    }

                    KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                        << "New knot vector: " << KnotsVRefined
                        << ", new weights vector: " << WeightsRefined << std::endl;

                    p_nurbs_volume->SetInternals(PointsRefined,
                        p_nurbs_volume->PolynomialDegreeU(), p_nurbs_volume->PolynomialDegreeV(), p_nurbs_volume->PolynomialDegreeW(),
                        p_nurbs_volume->KnotsU(), KnotsVRefined, p_nurbs_volume->KnotsW());
                    }
                }
                KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                        << "starting w knot insertion: " << std::endl;
                if (rParameters["parameters"].Has("insert_nb_per_span_w")) {
                    const IndexType nb_per_span_w = rParameters["parameters"]["insert_nb_per_span_w"].GetInt();

                    if (nb_per_span_w > 0) {
                    std::vector<double> spans_local_space;
                    p_nurbs_volume->SpansLocalSpace(spans_local_space, 2);

                    std::vector<double> knots_to_insert_w;

                    // if (rParameters["local_refinement_parameters"]["do_you_want_to_refine_inner"].GetBool()) {
                    //     for (IndexType i = 0; i < spans_local_space.size() - 1; ++i) {
                    //         const double delta_v = (spans_local_space[i + 1] - spans_local_space[i]) / (nb_per_span_v + 1);
                    //         for (IndexType j = 1; j < nb_per_span_v + 1; ++j) {
                    //             double v_local = spans_local_space[i] + delta_v*j;
                    //             if ( (v_local > box_refinement_u_v_min_max[2] 
                    //                     || std::abs(v_local - box_refinement_u_v_min_max[2]) < 1e-13 )
                    //                     && v_local < box_refinement_u_v_min_max[3]) {
                    //                 v_local = box_refinement_u_v_min_max[2];
                    //                 while (v_local < box_refinement_u_v_min_max[3]) {
                    //                     knots_to_insert_v.push_back(v_local);
                    //                     v_local += delta_v_refined;
                    //                 }
                    //                 if (std::abs(knots_to_insert_v[knots_to_insert_v.size()-1] - box_refinement_u_v_min_max[3]) > 1e-13){
                    //                     knots_to_insert_v.push_back(box_refinement_u_v_min_max[3]);
                    //                 }
                    //                 j = ((v_local - spans_local_space[i]) / delta_v) +1; 
                    //                 // if (std::abs(spans_local_space[i] + delta_v*j - v_local) < 1e-13) { j++;}
                    //             }
                    //             knots_to_insert_v.push_back(spans_local_space[i] + delta_v * j);
                    //         }
                    //     }
                    // } else {
                    for (IndexType i = 0; i < spans_local_space.size() - 1; ++i) {
                        const double delta_w = (spans_local_space[i + 1] - spans_local_space[i]) / (nb_per_span_w + 1);
                        for (IndexType j = 1; j < nb_per_span_w + 1; ++j) {
                            knots_to_insert_w.push_back(spans_local_space[i] + delta_w * j);
                        }
                    }
                    // }

                    KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                        << "Refining nurbs volume #" << p_nurbs_volume->Id() << " by inserting knots in w: "
                        << knots_to_insert_w << std::endl;

                    PointerVector<NodeType> PointsRefined;
                    Vector KnotsWRefined;
                    Vector WeightsRefined;

                    NurbsVolumeRefinementUtilities::KnotRefinementW(*(p_nurbs_volume.get()), knots_to_insert_w,
                        PointsRefined, KnotsWRefined);
                    

                    // Recreate nodes in model part to ensure correct assignment of dofs
                    IndexType node_id = (r_model_part.GetParentModelPart().NodesEnd() - 1)->Id() + 1;
                    for (IndexType i = 0; i < PointsRefined.size(); ++i) {
                        if (PointsRefined(i)->Id() == 0) {
                            PointsRefined(i) = r_model_part.CreateNewNode(node_id, PointsRefined[i][0], PointsRefined[i][1], PointsRefined[i][2]);
                            node_id++;
                        }
                    }

                    KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                        << "New knot vector: " << KnotsWRefined
                        << ", new weights vector: " << WeightsRefined << std::endl;

                    p_nurbs_volume->SetInternals(PointsRefined,
                        p_nurbs_volume->PolynomialDegreeU(), p_nurbs_volume->PolynomialDegreeV(), p_nurbs_volume->PolynomialDegreeW(),
                        p_nurbs_volume->KnotsU(), p_nurbs_volume->KnotsV(), KnotsWRefined);
                    }
                    else {
                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Trying to refine nurbs surface #" << p_nurbs_volume->Id()
                            << " by inserting knots in v, however insert_nb_per_span_v is set to 0." << std::endl;
                    }
                }
            }
        } 
        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                        << "end refinement modeler" <<  std::endl;
    }

    void RefinementModeler::GetGeometryList(
        GeometriesArrayType& rGeometryList,
        ModelPart& rModelPart,
        const Parameters rParameters) const
    {
        if (rParameters.Has("brep_id")) {
            rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_id"].GetInt()));
        }
        if (rParameters.Has("brep_ids")) {
            for (SizeType i = 0; i < rParameters["brep_ids"].size(); ++i) {
                rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_ids"][i].GetInt()));
            }
        }
        if (rParameters.Has("brep_name")) {
            rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_name"].GetString()));
        }
        if (rParameters.Has("brep_names")) {
            for (SizeType i = 0; i < rParameters["brep_names"].size(); ++i) {
                rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_names"][i].GetString()));
            }
        }

        KRATOS_ERROR_IF(rGeometryList.size() == 0)
            << "Empty geometry list. Either \"brep_id\", \"brep_ids\", \"brep_name\" or \"brep_names\" are the possible options." << std::endl;
    }

    /// Reads in a json formatted file and returns its KratosParameters instance.
    Parameters RefinementModeler::ReadParamatersFile(
        const std::string& rDataFileName) const
    {
        // Check if rDataFileName ends with ".cad.json" and add it if needed.
        const std::string data_file_name = (rDataFileName.compare(rDataFileName.size() - 9, 9, ".iga.json") != 0)
            ? rDataFileName + ".iga.json"
            : rDataFileName;

        std::ifstream infile(data_file_name);
        KRATOS_ERROR_IF_NOT(infile.good()) << "Physics fil: "
            << data_file_name << " cannot be found." << std::endl;
        KRATOS_INFO_IF("ReadParamatersFile", mEchoLevel > 3)
            << "Reading file: \"" << data_file_name << "\"" << std::endl;

        std::stringstream buffer;
        buffer << infile.rdbuf();

        return Parameters(buffer.str());
    }

    ///@}
}
