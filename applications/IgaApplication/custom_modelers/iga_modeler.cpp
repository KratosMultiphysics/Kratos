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
#include "iga_modeler.h"


namespace Kratos
{
    ///@name Stages
    ///@{

    void IgaModeler::SetupModelPart()
    {
        KRATOS_ERROR_IF_NOT(mParameters.Has("cad_model_part_name"))
            << "Missing \"cad_model_part\" section" << std::endl;
        ModelPart& cad_model_part = mpModel->GetModelPart(mParameters["cad_model_part_name"].GetString());
        KRATOS_ERROR_IF_NOT(mParameters.Has("analysis_model_part_name"))
            << "Missing \"analysis_model_part_name\" section" << std::endl;
        ModelPart& analysis_model_part = mpModel->GetModelPart(mParameters["analysis_model_part_name"].GetString());

        const std::string& rDataFileName = mParameters.Has("physics_file_name")
            ? mParameters["physics_file_name"].GetString()
            : "physics.iga.json";

        const Parameters iga_physics_parameters = ReadParamatersFile(rDataFileName);

        CreateIntegrationDomain(
            cad_model_part,
            analysis_model_part,
            iga_physics_parameters);
    }

    ///@}

    void IgaModeler::CreateIntegrationDomain(
        ModelPart& rCadModelPart,
        ModelPart& rModelPart,
        const Parameters rParameters) const
    {
        KRATOS_ERROR_IF_NOT(rParameters.Has("element_condition_list"))
            << "Missing \"element_condition_list\" section" << std::endl;

        KRATOS_ERROR_IF_NOT(rParameters["element_condition_list"].IsArray())
            << "\"element_condition_list\" needs to be an array." << std::endl;

        for (SizeType i = 0; i < rParameters["element_condition_list"].size(); ++i)
        {
            CreateIntegrationDomainPerUnit(
                rCadModelPart,
                rModelPart,
                rParameters["element_condition_list"][i]);
        }
    }

    void IgaModeler::CreateIntegrationDomainPerUnit(
        ModelPart& rCadModelPart,
        ModelPart& rModelPart,
        const Parameters rParameters) const
    {
        KRATOS_ERROR_IF_NOT(rParameters.Has("iga_model_part"))
            << "\"iga_model_part\" need to be specified." << std::endl;

        KRATOS_ERROR_IF_NOT(rParameters.Has("parameters"))
            << "\"parameters\" need to be specified." << std::endl;

        std::string sub_model_part_name = rParameters["iga_model_part"].GetString();

        ModelPart& sub_model_part = rModelPart.HasSubModelPart(sub_model_part_name)
            ? rModelPart.GetSubModelPart(sub_model_part_name)
            : rModelPart.CreateSubModelPart(sub_model_part_name);

        // Generate the list of geometries, which are needed, here.
        GeometriesArrayType geometry_list;
        GetCadGeometryList(geometry_list, rCadModelPart, rParameters);

        if (!rParameters.Has("geometry_type")) {
            CreateQuadraturePointGeometries(
                geometry_list, sub_model_part, rParameters["parameters"]);
        }
        else {
            std::string geometry_type = rParameters["geometry_type"].GetString();
            if (geometry_type == "GeometrySurfaceNodes"
                || geometry_type == "GeometrySurfaceVariationNodes"
                || geometry_type == "GeometryCurveNodes"
                || geometry_type == "GeometryCurveVariationNodes") {
                GetPointsAt(geometry_list, geometry_type, rParameters["parameters"], sub_model_part);
            }
            else {
                CreateQuadraturePointGeometries(
                    geometry_list, sub_model_part, rParameters["parameters"]);
            }
        }
        KRATOS_INFO_IF("CreateIntegrationDomainElementCondition", mEchoLevel > 3)
            << "Creation of elements/ conditions finished in: " << sub_model_part << std::endl;
    }

    void IgaModeler::CreateQuadraturePointGeometries(
        GeometriesArrayType& rGeometryList,
        ModelPart& rModelPart,
        const Parameters rParameters) const
    {
        KRATOS_ERROR_IF_NOT(rParameters.Has("type"))
            << "\"type\" need to be specified." << std::endl;
        std::string type = rParameters["type"].GetString();
        KRATOS_ERROR_IF_NOT(rParameters.Has("name"))
            << "\"name\" need to be specified." << std::endl;
        std::string name = rParameters["name"].GetString();

        SizeType shape_function_derivatives_order = 1;
        if (rParameters.Has("shape_function_derivatives_order")) {
            shape_function_derivatives_order = rParameters["shape_function_derivatives_order"].GetInt();
        }
        else {
            KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 4)
                << "shape_function_derivatives_order is not provided and thus being considered as 1. " << std::endl;
        }

        KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 0)
            << "Creating " << name << "s of type: " << type
            << " for " << rGeometryList.size() << " geometries"
            << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

        for (SizeType i = 0; i < rGeometryList.size(); ++i)
        {
            GeometriesArrayType geometries;
            IntegrationInfo integration_info = rGeometryList[i].GetDefaultIntegrationInfo();
            rGeometryList[i].CreateQuadraturePointGeometries(
                geometries, shape_function_derivatives_order, integration_info);

            KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 1)
                << geometries.size() << " quadrature point geometries have been created." << std::endl;

            if (type == "element" || type == "Element") {
                SizeType id = 1;
                if (rModelPart.GetRootModelPart().Elements().size() > 0)
                    id = rModelPart.GetRootModelPart().Elements().back().Id() + 1;

                this->CreateElements(
                    geometries.ptr_begin(), geometries.ptr_end(),
                    rModelPart, name, id, PropertiesPointerType());
            }
            else if (type == "condition" || type == "Condition") {
                SizeType id = 1;
                if (rModelPart.GetRootModelPart().Conditions().size() > 0)
                    id = rModelPart.GetRootModelPart().Conditions().back().Id() + 1;

                this->CreateConditions(
                    geometries.ptr_begin(), geometries.ptr_end(),
                    rModelPart, name, id, PropertiesPointerType());
            }
            else {
                KRATOS_ERROR << "\"type\" does not exist: " << type
                    << ". Possible types are \"element\" and \"condition\"." << std::endl;
            }
        }
    }

    ///@}
    ///@name CAD functionalities
    ///@{

    void IgaModeler::GetCadGeometryList(
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

    ///@}
    ///@name Generate Elements and Conditions
    ///@{

    void IgaModeler::CreateElements(
        typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
        typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
        ModelPart& rModelPart,
        std::string& rElementName,
        SizeType& rIdCounter,
        PropertiesPointerType pProperties) const
    {
        KRATOS_ERROR_IF(!KratosComponents<Element>::Has(rElementName))
            << rElementName << " not registered." << std::endl;

        const Element& rReferenceElement = KratosComponents<Element>::Get(rElementName);

        ElementsContainerType new_element_list;

        KRATOS_INFO_IF("CreateElements", mEchoLevel > 2)
            << "Creating elements of type " << rElementName
            << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

        for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it)
        {
            new_element_list.push_back(
                rReferenceElement.Create(rIdCounter, (*it), pProperties));
            for (SizeType i = 0; i < (*it)->size(); ++i) {
                rModelPart.AddNode((*it)->pGetPoint(i));
            }
            rIdCounter++;
        }

        rModelPart.AddElements(new_element_list.begin(), new_element_list.end());
    }

    void IgaModeler::CreateConditions(
        typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
        typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
        ModelPart& rModelPart,
        std::string& rConditionName,
        SizeType& rIdCounter,
        PropertiesPointerType pProperties) const
    {
        const Condition& rReferenceCondition = KratosComponents<Condition>::Get(rConditionName);

        ModelPart::ConditionsContainerType new_condition_list;

        KRATOS_INFO_IF("CreateConditions", mEchoLevel > 2)
            << "Creating conditions of type " << rConditionName
            << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

        for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it)
        {
            new_condition_list.push_back(
                rReferenceCondition.Create(rIdCounter, (*it), pProperties));
            for (SizeType i = 0; i < (*it)->size(); ++i) {
                rModelPart.AddNode((*it)->pGetPoint(i));
            }
            rIdCounter++;
        }

        rModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());
    }

    /// Reads in a json formatted file and returns its KratosParameters instance.
    Parameters IgaModeler::ReadParamatersFile(
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
    ///@name Get Points at Boundaries
    ///@{

    /// Searches points at boundaries of nurbs geometries.
    void IgaModeler::GetPointsAt(
        GeometriesArrayType& rGeometryList,
        const std::string& rGeometryType,
        const Parameters rParameters,
        ModelPart& rModelPart) const
    {
        for (auto& geom : rGeometryList) {

            // local_coordinates 0->Beginn
            //                   1->End
            //                  -1->All nodes in this dimension
            Vector local_coordinates = rParameters["local_parameters"].GetVector();

            auto p_background_geometry = geom.pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);

            if (rGeometryType == "GeometryCurveNodes") {
                KRATOS_DEBUG_ERROR_IF(geom.Dimension() != 1) << "Geometry #" << geom.Id()
                    << " needs to have a dimension of 1 for type GeometryCurveNodes. Dimension: " << geom.Dimension()
                    << ". Geometry" << geom << std::endl;

                SizeType number_of_cps = p_background_geometry->size();

                IndexType t_start = 0;
                IndexType t_end = number_of_cps;

                if (local_coordinates[0] >= 0) {
                    t_start = local_coordinates[0] * (number_of_cps - 1);
                    t_end = local_coordinates[0] * (number_of_cps - 1) + 1;
                }

                for (IndexType i = t_start; i < t_end; ++i) {
                    rModelPart.AddNode(p_background_geometry->pGetPoint(i));
                }
            }
            else if (rGeometryType == "GeometryCurveVariationNodes") {
                KRATOS_DEBUG_ERROR_IF(geom.Dimension() != 1) << "Geometry #" << geom.Id()
                    << " needs to have a dimension of 1 for type GeometryCurveVariationNodes. Dimension: " << geom.Dimension()
                    << ". Geometry" << geom << std::endl;

                SizeType number_of_cps = p_background_geometry->size();

                KRATOS_ERROR_IF(number_of_cps < 3)
                    << "GetPointsAt: Not enough control points to get second row of nodes."
                    << std::endl;

                if (local_coordinates[0] == 0) {
                    rModelPart.AddNode(p_background_geometry->pGetPoint(1));
                }
                else if (local_coordinates[0] == 1) {
                    rModelPart.AddNode(p_background_geometry->pGetPoint(number_of_cps - 1));
                }
                else {
                    KRATOS_ERROR << "GetPointsAt: GeometrySurfaceVariationNodes and local coordinates: " << local_coordinates[0]
                        << " is no available option. Only 0 and 1 are possible with this combination." << std::endl;
                }
            }
            else if (rGeometryType == "GeometrySurfaceNodes"
                || rGeometryType == "GeometrySurfaceVariationNodes") {
                KRATOS_DEBUG_ERROR_IF(geom.Dimension() != 2) << "Geometry #" << geom.Id()
                    << " needs to have a dimension of 2 for type " << rGeometryType << ". Dimension: "
                    << geom.Dimension() << ". Geometry" << geom << std::endl;

                SizeType number_of_cps_u = geom.PointsNumberInDirection(0);
                SizeType number_of_cps_v = geom.PointsNumberInDirection(1);

                IndexType u_start = 0;
                IndexType u_end = number_of_cps_u;
                IndexType v_start = 0;
                IndexType v_end = number_of_cps_v;

                if (rGeometryType == "GeometrySurfaceNodes") {
                    if (local_coordinates[0] >= 0) {
                        u_start = local_coordinates[0] * (number_of_cps_u - 1);
                        u_end = local_coordinates[0] * (number_of_cps_u - 1) + 1;
                    }
                    if (local_coordinates[1] >= 0) {
                        v_start = local_coordinates[1] * (number_of_cps_v - 1);
                        v_end = local_coordinates[1] * (number_of_cps_v - 1) + 1;
                    }
                }
                else if (rGeometryType == "GeometrySurfaceVariationNodes") {
                    if (local_coordinates[0] == 0) {
                        u_start = 1;
                        u_end = 2;
                    }
                    if (local_coordinates[0] == 1) {
                        u_start = number_of_cps_u - 2;
                        u_end = number_of_cps_u - 1;
                    }
                    if (local_coordinates[1] == 0) {
                        v_start = 1;
                        v_end = 2;
                    }
                    if (local_coordinates[1] == 1) {
                        v_start = number_of_cps_v - 2;
                        v_end = number_of_cps_v - 1;
                    }
                }

                for (IndexType i = u_start; i < u_end; ++i) {
                    for (IndexType j = v_start; j < v_end; ++j) {
                        rModelPart.AddNode(p_background_geometry->pGetPoint(i + j * number_of_cps_u));
                    }
                }
            }
        }
    }

    ///@}
}
