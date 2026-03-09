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
                        for (auto& r_point : p_nurbs_surface->Points()) {
                            r_point.Set(TO_ERASE, true);
                        }
                        r_model_part.RemoveNodesFromAllLevels(TO_ERASE);
                        IndexType node_id = r_model_part.Nodes().empty() ? 1 : (r_model_part.NodesEnd() - 1)->Id() + 1;

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
                        for (auto& r_point : p_nurbs_surface->Points()) {
                            r_point.Set(TO_ERASE, true);
                        }
                        r_model_part.RemoveNodesFromAllLevels(TO_ERASE);
                        IndexType node_id = r_model_part.Nodes().empty() ? 1 : (r_model_part.NodesEnd() - 1)->Id() + 1;

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
                if (rParameters["parameters"].Has("insert_nb_per_span_u")) {
                    const IndexType nb_per_span_u = rParameters["parameters"]["insert_nb_per_span_u"].GetInt();

                    if (nb_per_span_u > 0) {
                        std::vector<double> spans_local_space;
                        p_nurbs_surface->SpansLocalSpace(spans_local_space, 0);

                        std::vector<double> knots_to_insert_u;
                        for (IndexType i = 0; i < spans_local_space.size() - 1; ++i) {
                            const double delta_u = (spans_local_space[i + 1] - spans_local_space[i]) / (nb_per_span_u + 1);
                            for (IndexType j = 1; j < nb_per_span_u + 1; ++j) {
                                knots_to_insert_u.push_back(spans_local_space[i] + delta_u * j);
                            }
                        }

                        KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                            << "Refining nurbs surface #" << p_nurbs_surface->Id() << " by inserting knots in u: "
                            << knots_to_insert_u << std::endl;

                        PointerVector<NodeType> PointsRefined;
                        Vector KnotsURefined;
                        Vector WeightsRefined;

                        NurbsSurfaceRefinementUtilities::KnotRefinementU(*(p_nurbs_surface.get()), knots_to_insert_u,
                            PointsRefined, KnotsURefined, WeightsRefined);

                        // Recreate nodes in model part to ensure correct assignment of dofs
                        for (auto& r_point : p_nurbs_surface->Points()) {
                            r_point.Set(TO_ERASE, true);
                        }
                        r_model_part.RemoveNodesFromAllLevels(TO_ERASE);
                        IndexType node_id = r_model_part.Nodes().empty() ? 1 : (r_model_part.NodesEnd() - 1)->Id() + 1;
                        
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
                if (rParameters["parameters"].Has("insert_nb_per_span_v")) {
                    const IndexType nb_per_span_v = rParameters["parameters"]["insert_nb_per_span_v"].GetInt();

                    if (nb_per_span_v > 0) {
                    std::vector<double> spans_local_space;
                    p_nurbs_surface->SpansLocalSpace(spans_local_space, 1);

                    std::vector<double> knots_to_insert_v;
                    for (IndexType i = 0; i < spans_local_space.size() - 1; ++i) {
                        const double delta_v = (spans_local_space[i + 1] - spans_local_space[i]) / (nb_per_span_v + 1);
                        for (IndexType j = 1; j < nb_per_span_v + 1; ++j) {
                            knots_to_insert_v.push_back(spans_local_space[i] + delta_v * j);
                        }
                    }

                    KRATOS_INFO_IF("::[RefinementModeler]::ApplyRefinement", mEchoLevel > 1)
                        << "Refining nurbs surface #" << p_nurbs_surface->Id() << " by inserting knots in v: "
                        << knots_to_insert_v << std::endl;

                    PointerVector<NodeType> PointsRefined;
                    Vector KnotsVRefined;
                    Vector WeightsRefined;

                    NurbsSurfaceRefinementUtilities::KnotRefinementV(*(p_nurbs_surface.get()), knots_to_insert_v,
                        PointsRefined, KnotsVRefined, WeightsRefined);

                    // Recreate nodes in model part to ensure correct assignment of dofs
                    for (auto& r_point : p_nurbs_surface->Points()) {
                        r_point.Set(TO_ERASE, true);
                    }
                    r_model_part.RemoveNodesFromAllLevels(TO_ERASE);
                    IndexType node_id = r_model_part.Nodes().empty() ? 1 : (r_model_part.NodesEnd() - 1)->Id() + 1;

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
            }
        }

        // Reassign node ids from 1 to n after refinement
        std::size_t new_id = 1;
        for (auto& r_node : r_model_part.Nodes()) {
            r_node.SetId(new_id);
            new_id++;
        }
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
