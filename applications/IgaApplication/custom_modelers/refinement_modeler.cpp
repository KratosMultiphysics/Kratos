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
                //auto nurbs_surface = (geometry_list[n].size() > 0)
                //    ? static_cast<NurbsSurfaceGeometry<3, PointerVector<Node<3>>>>(geometry_list[n])
                //    : static_cast<NurbsSurfaceGeometry<3, PointerVector<Node<3>>>>(geometry_list[n].GetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX));
                auto p_nurbs_surface = dynamic_pointer_cast<NurbsSurfaceGeometry<3, PointerVector<Node<3>>>>(geometry_list[n].pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX));
                if (rParameters["parameters"].Has("insert_nb_per_span_u")) {
                    const IndexType nb_per_span_u = rParameters["parameters"]["insert_nb_per_span_u"].GetInt();

                    std::vector<double> spans_local_space;
                    p_nurbs_surface->Spans(spans_local_space, 0);

                    std::vector<double> knots_to_insert_u;
                    for (IndexType i = 0; i < spans_local_space.size() - 1; ++i) {
                        const double delta_u = (spans_local_space[i] - spans_local_space[i + 1]) / (nb_per_span_u + 1);
                        for (IndexType j = 0; j < nb_per_span_u; ++j) {
                            knots_to_insert_u.push_back(spans_local_space[i] + delta_u * j);
                        }
                    }

                    PointerVector<NodeType> PointsRefined;
                    Vector KnotsURefined;
                    Vector WeightsRefined;

                    NurbsSurfaceRefinementUtilities::KnotRefinementU(*(p_nurbs_surface.get()), knots_to_insert_u,
                        PointsRefined, KnotsURefined, WeightsRefined);

                    // Recreate nodes in model part to ensure correct assignment of dofs
                    IndexType node_id = (r_model_part.NodesEnd() - 1)->Id() + 1;
                    for (IndexType i = 0; i < PointsRefined.size(); ++i) {
                        if (PointsRefined(i)->Id() == 0) {
                            PointsRefined(i) = r_model_part.CreateNewNode(node_id, PointsRefined[i][0], PointsRefined[i][1], PointsRefined[i][2]);
                            node_id++;
                        }
                    }

                    p_nurbs_surface->SetInternals(PointsRefined,
                        p_nurbs_surface->PolynomialDegreeU(), p_nurbs_surface->PolynomialDegreeV(),
                        KnotsURefined, p_nurbs_surface->KnotsV(),
                        WeightsRefined);
                }
                if (rParameters["parameters"].Has("insert_nb_per_span_v")) {
                    const IndexType nb_per_span_v = rParameters["parameters"]["insert_nb_per_span_v"].GetInt();

                    std::vector<double> spans_local_space;
                    p_nurbs_surface->Spans(spans_local_space, 1);

                    std::vector<double> knots_to_insert_v;
                    for (IndexType i = 0; i < spans_local_space.size() - 1; ++i) {
                        const double delta_v = (spans_local_space[i] - spans_local_space[i + 1]) / (nb_per_span_v + 1);
                        for (IndexType j = 0; j < spans_local_space.size() - 1; ++j) {
                            knots_to_insert_v.push_back(spans_local_space[j] + delta_v);
                        }
                    }

                    PointerVector<NodeType> PointsRefined;
                    Vector KnotsURefined;
                    Vector WeightsRefined;

                    NurbsSurfaceRefinementUtilities::KnotRefinementU(*(p_nurbs_surface.get()), knots_to_insert_v,
                        PointsRefined, KnotsURefined, WeightsRefined);

                    // Recreate nodes in model part to ensure correct assignment of dofs
                    IndexType node_id = (r_model_part.NodesEnd() - 1)->Id() + 1;
                    for (IndexType i = 0; i < PointsRefined.size(); ++i) {
                        if (PointsRefined(i)->Id() == 0) {
                            PointsRefined(i) = r_model_part.CreateNewNode(node_id, PointsRefined[i][0], PointsRefined[i][1], PointsRefined[i][2]);
                            node_id++;
                        }
                    }

                    p_nurbs_surface->SetInternals(PointsRefined,
                        p_nurbs_surface->PolynomialDegreeU(), p_nurbs_surface->PolynomialDegreeV(),
                        KnotsURefined, p_nurbs_surface->KnotsV(),
                        WeightsRefined);
                }
            }
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
