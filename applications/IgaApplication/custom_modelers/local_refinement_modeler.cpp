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
#include "local_refinement_modeler.h"

namespace Kratos
{
    ///@name Stages
    ///@{

    void LocalRefinementModeler::PrepareGeometryModel()
    {
        const std::string DataFileName = mParameters.Has("local_refinements_file_name")
            ? mParameters["local_refinements_file_name"].GetString()
            : "local_refinements.iga.json";

        KRATOS_INFO_IF("::[LocalRefinementModeler]::", mEchoLevel > 0) << "Refining model by: " << DataFileName << std::endl;

        const Parameters local_refinements_parameters = ReadParamatersFile(DataFileName);

        ApplyLocalRefinements(local_refinements_parameters);
    }

    void LocalRefinementModeler::ApplyLocalRefinements(
        const Parameters rParameters) const
    {
        KRATOS_ERROR_IF_NOT(rParameters.Has("local_refinements"))
            << "Parameters do not have local_refinements section.\n"
            << rParameters << std::endl;

        KRATOS_ERROR_IF_NOT(rParameters["local_refinements"].IsArray())
            << "local refinements section need to be of type array.\n"
            << rParameters << std::endl;

        for (IndexType i = 0; i < rParameters["local_refinements"].size(); ++i) {
            ApplyLocalRefinement(rParameters["local_refinements"][i]);
        }
    }

    void LocalRefinementModeler::ApplyLocalRefinement(
        const Parameters rParameters) const
    {
        KRATOS_ERROR_IF_NOT(rParameters.Has("model_part_name"))
            << "Missing \"model_part_name\" in local refinements block.\n"
            << rParameters << std::endl;
        const std::string model_part_name = rParameters["model_part_name"].GetString();
        ModelPart& r_model_part = mpModel->HasModelPart(model_part_name)
            ? mpModel->GetModelPart(model_part_name)
            : mpModel->CreateModelPart(model_part_name);

        // Generate the list of geometries, which are needed, here.
        GeometriesArrayType geometry_list;
        GetGeometryList(geometry_list, r_model_part, rParameters);

        KRATOS_ERROR_IF_NOT(rParameters.Has("parameters"))
            << "Missing \"parameters\" in local refinements block.\n"
            << rParameters << std::endl;
        KRATOS_ERROR_IF_NOT(rParameters.Has("geometry_type"))
            << "Missing \"geometry_type\".\n"
            << rParameters << std::endl;
        if (rParameters["geometry_type"].GetString() == "THBSurface") {
            for (IndexType n = 0; n < geometry_list.size(); ++n) {
                auto p_thb_surface = (geometry_list[n].size() > 0)
                    ? dynamic_pointer_cast<THBSurfaceGeometry<3, PointerVector<Node>>>(geometry_list(n))
                    : dynamic_pointer_cast<THBSurfaceGeometry<3, PointerVector<Node>>>(geometry_list(n)->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX));

                if (rParameters["parameters"].Has("ref_boxes")) {
                    Vector temp_ref_boxes = rParameters["parameters"]["ref_boxes"].GetVector();
                    std::vector<int> ref_boxes(temp_ref_boxes.size()); 

                    for (size_t i = 0; i < temp_ref_boxes.size(); ++i) {
                        ref_boxes[i] = static_cast<int>(temp_ref_boxes[i]); 
                    }

                    if (ref_boxes.size() > 0) {
                        KRATOS_INFO_IF("::[LocalRefinementModeler]::ApplyLocalRefinement", mEchoLevel > 1)
                            << "Refining thb surface #" << p_thb_surface->Id() << " by ref boxes"
                            << std::endl;

                        p_thb_surface->SetInternals(ref_boxes);

                    } else {
                        KRATOS_INFO_IF("::[LocalRefinementModeler]::ApplyLocalRefinement", mEchoLevel > 1)
                            << "Trying to refine THB surface #" << p_thb_surface->Id()
                            << " by elevating degree in u, however increase_degree_u is set to 0." << std::endl;
                    }
                }
            }
        }
    }

    void LocalRefinementModeler::GetGeometryList(
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
    Parameters LocalRefinementModeler::ReadParamatersFile(
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
