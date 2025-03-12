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
#include "integration/integration_point_utilities.h"
#include "iga_application_variables.h"
#include "includes/global_pointer_variables.h"


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

        if (mParameters.Has("integrate_on_true_boundary")) {
            if (mParameters["integrate_on_true_boundary"].GetBool()) {
                // method for computing the integral of the solution along the true boundary
                prepareIntegrationOnTrueBoundary(analysis_model_part);
            } 
        }
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
                geometry_list, sub_model_part, rParameters["parameters"], std::string{});
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
                //Check if is isSBM 
                if (rParameters["parameters"].Has("SBM_parameters")) {
                    CreateQuadraturePointGeometriesSbm(
                        geometry_list, rModelPart, sub_model_part, rParameters["parameters"], geometry_type);
                }
                else {
                    CreateQuadraturePointGeometries(
                        geometry_list, sub_model_part, rParameters["parameters"], geometry_type);
                }
            }
        }
        KRATOS_INFO_IF("CreateIntegrationDomainElementCondition", mEchoLevel > 3)
            << "Creation of elements/ conditions finished in: " << sub_model_part << std::endl;
    }

    void IgaModeler::CreateQuadraturePointGeometries(
        GeometriesArrayType& rGeometryList,
        ModelPart& rModelPart,
        const Parameters rParameters,
        std::string GeometryType) const
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

        std::string quadrature_method = rParameters.Has("quadrature_method")
            ? rParameters["integration_rule"].GetString()
            : "GAUSS";

        KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 0)
            << "Creating " << name << "s of type: " << type
            << " for " << rGeometryList.size() << " geometries"
            << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

        
        for (SizeType i = 0; i < rGeometryList.size(); ++i)
        {
            GeometriesArrayType geometries;
            IntegrationInfo integration_info = rGeometryList[i].GetDefaultIntegrationInfo();
            for (IndexType i = 0; i < integration_info.LocalSpaceDimension(); ++i) {
                if (quadrature_method == "GAUSS") {
                    integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GAUSS);
                }
                else if (quadrature_method == "EXTENDED_GAUSS") {
                    integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::EXTENDED_GAUSS);
                }
                else if (quadrature_method == "GRID") {
                    integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GRID);
                }
                else {
                    KRATOS_INFO("CreateQuadraturePointGeometries") << "Quadrature method: " << quadrature_method
                        << " is not available. Available options are \"GAUSS\" and \"GRID\". Default quadrature method is being considered." << std::endl;
                }
            }

            if (rParameters.Has("number_of_integration_points_per_span")) {
                for (IndexType i = 0; i < integration_info.LocalSpaceDimension(); ++i) {
                    integration_info.SetNumberOfIntegrationPointsPerSpan(i, rParameters["number_of_integration_points_per_span"].GetInt());
                }
            }
            if (GeometryType == "SurfaceEdge"
                && rGeometryList[i].GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Coupling_Geometry)
            {
                rGeometryList[i].GetGeometryPart(0).CreateQuadraturePointGeometries(
                    geometries, shape_function_derivatives_order, integration_info);
            }
            else
            {
                rGeometryList[i].CreateQuadraturePointGeometries(
                    geometries, shape_function_derivatives_order, integration_info);
            }

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
                std::vector<int> listIdClosestCondition(geometries.size());
                this->CreateConditions(
                    geometries.ptr_begin(), geometries.ptr_end(),
                    rModelPart, name, id, PropertiesPointerType());
            }
            else {
                KRATOS_ERROR << "\"type\" does not exist: " << type
                    << ". Possible types are \"element\" and \"condition\"." << std::endl;
            }
        }
        KRATOS_WATCH("finee")
    }




    void IgaModeler::CreateQuadraturePointGeometriesSbm(
        GeometriesArrayType& rGeometryList,
        ModelPart& rIgaModelPart,
        ModelPart& rModelPart,
        const Parameters rParameters,
        std::string GeometryType) const
    {
        KRATOS_ERROR_IF_NOT(rParameters.Has("type"))
            << "\"type\" need to be specified." << std::endl;
        std::string type = rParameters["type"].GetString();
        KRATOS_ERROR_IF_NOT(rParameters.Has("name"))
            << "\"name\" need to be specified." << std::endl;
        std::string name = rParameters["name"].GetString();

        std::string body_fitted_condition_name = rParameters.Has("body_fitted_condition_name")
                                                ? rParameters["body_fitted_condition_name"].GetString() : "";

        SizeType shape_function_derivatives_order = 1;
        if (rParameters.Has("shape_function_derivatives_order")) {
            shape_function_derivatives_order = rParameters["shape_function_derivatives_order"].GetInt();
        }
        else {
            KRATOS_INFO_IF("CreateQuadraturePointGeometriesSbm", mEchoLevel > 4)
                << "shape_function_derivatives_order is not provided and thus being considered as 1. " << std::endl;
        }

        std::string quadrature_method = rParameters.Has("quadrature_method")
            ? rParameters["integration_rule"].GetString()
            : "GAUSS";

        KRATOS_INFO_IF("CreateQuadraturePointGeometriesSbm", mEchoLevel > 0)
            << "Creating " << name << "s of type: " << type
            << " for " << rGeometryList.size() << " geometries"
            << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

        PointVector points;
        std::string skin_model_part_name;
        if (!mParameters.Has("skin_model_part_name")) skin_model_part_name = "skin_model_part";
        else {
            skin_model_part_name = mParameters["skin_model_part_name"].GetString();
        }
        ModelPart& skin_model_part = mpModel->HasModelPart(skin_model_part_name)
                ? mpModel->GetModelPart(skin_model_part_name)
                : KRATOS_ERROR << "::[CreateQuadraturePointGeometriesSbm]::: Sbm case -> skin_model_part has not been defined before."
                               << "Maybe you are not calling the nurbs_modeler_sbm" << std::endl;
        ModelPart& skin_sub_model_part_in = skin_model_part.GetSubModelPart("inner");
        ModelPart& skin_sub_model_part_out = skin_model_part.GetSubModelPart("outer");

        std::string surrogate_sub_model_part_name;
        if (!(rParameters["SBM_parameters"]).Has("surrogate_sub_model_part_name")) surrogate_sub_model_part_name = "da_cambiare";
        else {
            surrogate_sub_model_part_name = rParameters["SBM_parameters"]["surrogate_sub_model_part_name"].GetString();
        }

        ModelPart& surrogate_sub_model_part = rIgaModelPart.HasSubModelPart(surrogate_sub_model_part_name)
                ? rIgaModelPart.GetSubModelPart(surrogate_sub_model_part_name)
                : KRATOS_ERROR << "::[CreateQuadraturePointGeometriesSbm]::: Sbm case -> surrogate_sub_model_part has not been defined before."
                               << "Maybe you are not calling the nurbs_modeler_sbm" << std::endl;

        bool is_inner;     
        Vector meshSizes_uv;
        Vector parameterExternalCoordinates;
        double radius;
        
        is_inner = rParameters["SBM_parameters"]["is_inner"].GetBool();
        double meshSize;
        if (is_inner) { // INNER
            for (auto &i_cond : skin_sub_model_part_in.Conditions()) {
                points.push_back(PointTypePointer(new PointType(i_cond.Id(), i_cond.GetGeometry().Center().X(), i_cond.GetGeometry().Center().Y(), i_cond.GetGeometry().Center().Z())));
            }
        } 
        else { // OUTER
            for (auto &i_cond : skin_sub_model_part_out.Conditions()) {
                points.push_back(PointTypePointer(new PointType(i_cond.Id(), i_cond.GetGeometry().Center().X(), i_cond.GetGeometry().Center().Y(), i_cond.GetGeometry().Center().Z())));
            }
        }
        meshSizes_uv = surrogate_sub_model_part.GetValue(MARKER_MESHES);
        rIgaModelPart.SetValue(MARKER_MESHES, meshSizes_uv);
        parameterExternalCoordinates = surrogate_sub_model_part.GetValue(LOAD_MESHES);
        rIgaModelPart.SetValue(LOAD_MESHES,parameterExternalCoordinates);

        meshSize = meshSizes_uv[0];
        if (meshSizes_uv[1] > meshSize) {meshSize = meshSizes_uv[1];}
        if (meshSizes_uv.size() > 2) {if (meshSizes_uv[2] > meshSize) {meshSize = meshSizes_uv[2];}}
        

        radius = sqrt(2)*(meshSize);  //TODO: sqrt(3) in 3D
        // radius = 30*(meshSize);
        DynamicBins testBins(points.begin(), points.end());
        const int numberOfResults = 1e6; 
        ModelPart::NodesContainerType::ContainerType Results(numberOfResults);
        std::vector<double> list_of_distances(numberOfResults);
        for (SizeType i = 0; i < rGeometryList.size(); ++i)
        {
            GeometriesArrayType geometries;
            IntegrationInfo integration_info = rGeometryList[i].GetDefaultIntegrationInfo();
            for (IndexType i = 0; i < integration_info.LocalSpaceDimension(); ++i) {
                if (quadrature_method == "GAUSS") {
                    integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GAUSS);
                }
                else if (quadrature_method == "EXTENDED_GAUSS") {
                    integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::EXTENDED_GAUSS);
                }
                else if (quadrature_method == "GRID") {
                    integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GRID);
                }
                else {
                    KRATOS_INFO("CreateQuadraturePointGeometries") << "Quadrature method: " << quadrature_method
                        << " is not available. Available options are \"GAUSS\" and \"GRID\". Default quadrature method is being considered." << std::endl;
                }
            }
            if (rParameters.Has("number_of_integration_points_per_span")) {
                for (IndexType i = 0; i < integration_info.LocalSpaceDimension(); ++i) {
                    integration_info.SetNumberOfIntegrationPointsPerSpan(i, rParameters["number_of_integration_points_per_span"].GetInt());
                }
            }
            if (GeometryType == "SurfaceEdge"
                && rGeometryList[i].GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Coupling_Geometry)
            {
                rGeometryList[i].GetGeometryPart(0).CreateQuadraturePointGeometries(
                    geometries, shape_function_derivatives_order, integration_info);
            }
            else
            {
                rGeometryList[i].CreateQuadraturePointGeometries(
                    geometries, shape_function_derivatives_order, integration_info);
            }

            KRATOS_INFO_IF("CreateQuadraturePointGeometriesSbm", mEchoLevel > 1)
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
                std::vector<int> listIdClosestCondition(geometries.size());

                Point gaussPoint = geometries[0].Center(); 

                if (geometries.size() > 1) gaussPoint = geometries[1].Center(); 

                bool isCoincidentToExternalParameterSpace = false;   
                  
                // isCoincidentToExternalParameterSpace = CheckIsOnExternalParameterSpace(gaussPoint, parameterExternalCoordinates); 
                
                if (name == "SBMCondition" && !isCoincidentToExternalParameterSpace) { 
                    for (auto j= 0; j < geometries.size() ; j++) {  

                        Point gaussPoint = geometries[j].Center(); 
                        
                        PointerType pointToSearch = PointerType(new PointType(10000, gaussPoint.X(), gaussPoint.Y(), gaussPoint.Z()));

                        // OLD SEARCH (not working well- maybe in the future)
                        // PointerType nearestPoint = testBins.SearchNearestPoint(*pointToSearch);
                        // listIdClosestCondition[j] = nearestPoint->Id();

                        int obtainedResults = testBins.SearchInRadius(*pointToSearch, radius, Results.begin(), list_of_distances.begin(), numberOfResults);

                        double minimum_distance=1e10;
                        int nearestNodeId;

                        double best_n_dot_distance = 1e10;
                        Node best_node;
                        for (int i_distance = 0; i_distance < obtainedResults; i_distance++) {
                            double new_distance = list_of_distances[i_distance];   
                            if (new_distance < minimum_distance) { 
                                minimum_distance = new_distance;
                                nearestNodeId = i_distance;
                                }
                        }
                        // for (int i_distance = 0; i_distance < obtainedResults; i_distance++) {
                        //     double new_distance = list_of_distances[i_distance];   
                        //     //FIXME:
                        //     auto& closest_condition = skin_model_part.GetCondition(Results[i_distance]->Id());
                        //     auto& closest_point = closest_condition.GetGeometry()[0];

                        //     Vector condition_tangent_vector = geometries[0].Center() - geometries[1].Center();
                        //     condition_tangent_vector /= norm_2(condition_tangent_vector);
                        //     Vector temp_normal = closest_point.GetValue(NORMAL);
                        //     double n_dot_distance = std::abs(inner_prod(temp_normal, condition_tangent_vector));

                        //     // Vector distance_vector = closest_point - gaussPoint;
                        //     // distance_vector /= norm_2(distance_vector);
                        //     // double n_dot_distance = std::abs(inner_prod(temp_normal, distance_vector));
                            
                        //     // KRATOS_WATCH(gaussPoint)
                        //     // KRATOS_WATCH(closest_point)
                        //     // KRATOS_WATCH(temp_normal)
                        //     // KRATOS_WATCH(condition_tangent_vector)
                        //     // KRATOS_WATCH(n_dot_distance)
                        //     // KRATOS_WATCH("---------------")

                        //     if (n_dot_distance < best_n_dot_distance && new_distance < minimum_distance *1.5)
                        //     {
                        //         best_n_dot_distance = n_dot_distance;
                        //         nearestNodeId = i_distance;
                        //         best_node = closest_point;
                        //     }
                        // }

                        auto closest_condition = skin_model_part.GetCondition(Results[nearestNodeId]->Id());
                        best_node = closest_condition.GetGeometry()[0];

                        std::ofstream outputFile("txt_files/Projection_Coordinates.txt", std::ios::app);
                        outputFile << best_node[0] << " " << best_node[1] << " "  << gaussPoint.X() << " " << gaussPoint.Y() <<"\n";
                        outputFile.close();
                        // if (is_inner)
                        //     exit(0);

                        if (obtainedResults == 0) {
                             KRATOS_WATCH("0 POINTS FOUND: EXIT")
                             KRATOS_WATCH(pointToSearch)
                             exit(0);}

                        
                        listIdClosestCondition[j] = Results[nearestNodeId]->Id();                          
                    }
                    if (is_inner) {
                        this->CreateConditions(
                        geometries.ptr_begin(), geometries.ptr_end(),
                        rModelPart, skin_sub_model_part_in, listIdClosestCondition, id, PropertiesPointerType(), is_inner, meshSizes_uv);
                    }
                    else{
                        this->CreateConditions(
                        geometries.ptr_begin(), geometries.ptr_end(),
                        rModelPart, skin_sub_model_part_out, listIdClosestCondition, id, PropertiesPointerType(), is_inner, meshSizes_uv);
                    }
                } else if (isCoincidentToExternalParameterSpace){
                    
                    if (body_fitted_condition_name == "")
                        KRATOS_ERROR << "::[IgaModeler]:: The body_fitted_condition_name has not been defined " 
                                     << "for an outer loop that touches the background domain boundary." << std::endl;
                    this->CreateConditions(
                        geometries.ptr_begin(), geometries.ptr_end(),
                        rModelPart, body_fitted_condition_name, id, PropertiesPointerType());
                }
                else { 
                    this->CreateConditions(
                        geometries.ptr_begin(), geometries.ptr_end(),
                        rModelPart, name, id, PropertiesPointerType());
                }
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

        static int starting_brep_ids;

        // get surrogate model part name
        std::string surrogate_model_part_name;
        if (!mParameters.Has("surrogate_model_part_name")) surrogate_model_part_name = "surrogate_model_part";
        else {
            surrogate_model_part_name = mParameters["surrogate_model_part_name"].GetString();
        }

        if (rParameters.Has("brep_id")) {
            rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_id"].GetInt()));
        }
        if (rParameters.Has("brep_ids")) {
            for (SizeType i = 0; i < rParameters["brep_ids"].size(); ++i) {
                rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_ids"][i].GetInt()));
            }
            // int lastIndex
            starting_brep_ids = rParameters["brep_ids"][rParameters["brep_ids"].size()-1].GetInt() + 1;

            // OUTER
            std::string conditionName = rParameters["iga_model_part"].GetString();
            if (conditionName.rfind("SBM", 0) == 0) { 
                ModelPart& surrogateModelPart_outer = rModelPart.GetSubModelPart("surrogate_outer");
                // ModelPart& surrogateModelPart_outer = mpModel->GetModelPart(surrogate_model_part_name +"_outer");
                if (surrogateModelPart_outer.Conditions().size() > 0) {
                    // 2D
                    if ((*surrogateModelPart_outer.ConditionsBegin()).GetGeometry().size() == 2) {
                        int sizeSurrogateLoop_outer = surrogateModelPart_outer.Nodes().size();
                        for (int j = 0; j < (sizeSurrogateLoop_outer-1); ++j) {
                            // Add the brep_ids of the internal boundary for SBMLaplacianCondition
                            rGeometryList.push_back(rModelPart.pGetGeometry(starting_brep_ids));
                            starting_brep_ids++;
                        }
                        
                    } else { // 3D
                        int sizeSurrogateLoop = surrogateModelPart_outer.Conditions().size(); //lastSurrogateNode.Id() - firstSurrogateNode.Id() + 1 ;
                        for (SizeType j = 0; j < sizeSurrogateLoop-1; ++j) {
                            // Add the brep_ids of the internal boundary for SBMLaplacianCondition
                            rGeometryList.push_back(rModelPart.pGetGeometry(starting_brep_ids));
                            starting_brep_ids++;
                        }
                    }     
                }
            }
        }
        else {
            // INNER   
            ModelPart& surrogateModelPart_inner = rModelPart.GetSubModelPart("surrogate_inner");
            // ModelPart& surrogateModelPart_inner = mpModel->GetModelPart(surrogate_model_part_name + "_inner");
            if (surrogateModelPart_inner.Elements().size() > 0) {

                int sizeSurrogateLoop = surrogateModelPart_inner.Nodes().size();

                for (int iel = 1; iel < surrogateModelPart_inner.Elements().size()+1; iel++) {
                    // Each element in the surrogate_model_part represents a surrogate boundary loop. First "node" is the initial ID of the first surrogate node and
                    // the second "node" is the last surrogate node of that loop. (We have done this in the case we have multiple surrogate boundaries and 1 model part)
                    Node& firstSurrogateNode = surrogateModelPart_inner.pGetElement(iel)->GetGeometry()[0]; // Element 1 because is the only surrogate loop
                    Node& lastSurrogateNode = surrogateModelPart_inner.pGetElement(iel)->GetGeometry()[1];  // Element 1 because is the only surrogate loop
                    int sizeSurrogateLoop = lastSurrogateNode.Id() - firstSurrogateNode.Id() + 1 ;

                    for (SizeType j = 0; j < sizeSurrogateLoop; ++j) {
                        // Add the brep_ids of the internal boundary for SBMLaplacianCondition
                        rGeometryList.push_back(rModelPart.pGetGeometry(starting_brep_ids));
                        starting_brep_ids++;
                    }
                }
            } else { // 3D
                int sizeSurrogateLoop = surrogateModelPart_inner.Conditions().size(); //lastSurrogateNode.Id() - firstSurrogateNode.Id() + 1 ;
                for (SizeType j = 0; j < sizeSurrogateLoop; ++j) {
                    // Add the brep_ids of the internal boundary for SBMLaplacianCondition
                    rGeometryList.push_back(rModelPart.pGetGeometry(starting_brep_ids));
                    starting_brep_ids++;
                }
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
                // rModelPart.AddNode((*it)->pGetPoint(i));
                rModelPart.Nodes().push_back((*it)->pGetPoint(i));
            }
            rIdCounter++;
        }

        rModelPart.AddElements(new_element_list.begin(), new_element_list.end());
    }

    void IgaModeler::CreateConditions(
        typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
        typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
        ModelPart& rModelPart,
        ModelPart& rSkinModelPart,
        std::vector<int>& listIdClosestCondition,
        SizeType& rIdCounter,
        PropertiesPointerType pProperties,
        bool isInner,
        Vector mesh_size) const
    {

        ModelPart::ConditionsContainerType new_condition_list;

        int countListClosestCondition = 0;
        bool is2D = true;
        // if(rSkinModelPart.GetCondition(listIdClosestCondition[0]).GetGeometry().size() > 2) { is2D = false;}

        

        KRATOS_INFO_IF("CreateConditions", mEchoLevel > 2)
            << "Creating conditions of type " << "SBM Condition"
            << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;


        // correct conditions by projections
        
        // count how many conditions for each type
        std::vector<std::string> bc_on_skin_projection_type;
        std::vector<int> count_bc_on_skin_projection_type;
        for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it) {
            std::string rConditionName = rSkinModelPart.GetCondition(listIdClosestCondition[countListClosestCondition]).GetValue(CONDITION_NAME);

            // Find the condition name in bc_on_skin_projection_type
            auto it_name = std::find(bc_on_skin_projection_type.begin(), bc_on_skin_projection_type.end(), rConditionName);
            if (it_name != bc_on_skin_projection_type.end()) {
                // Increment the count for the existing condition name
                size_t index = std::distance(bc_on_skin_projection_type.begin(), it_name);
                count_bc_on_skin_projection_type[index]++;
            } else {
                // Add new condition name and initialize its count
                bc_on_skin_projection_type.push_back(rConditionName);
                count_bc_on_skin_projection_type.push_back(1); // Initialize count to 1
            }
            countListClosestCondition++;
        }
        //--------------------
        std::string max_condition_name;
        std::string layer_name;
        int max_count = 0; 
        for (size_t i = 0; i < count_bc_on_skin_projection_type.size(); ++i) {
            if (count_bc_on_skin_projection_type[i] > max_count) {
                max_count = count_bc_on_skin_projection_type[i];
                max_condition_name = bc_on_skin_projection_type[i];
            }
        }
        
        // create a pool of conditions of the right type
        std::vector<int> list_id_closest_condition_of_correct_bc;
        for (IndexType i = 0; i < listIdClosestCondition.size(); i++) {
            int condId = listIdClosestCondition[i];
            std::string rConditionName = rSkinModelPart.GetCondition(listIdClosestCondition[i]).GetValue(CONDITION_NAME);
            if (rConditionName == max_condition_name) 
            {
                list_id_closest_condition_of_correct_bc.push_back(condId);
                layer_name = rSkinModelPart.GetCondition(listIdClosestCondition[i]).GetValue(LAYER_NAME);
            }
        }
        KRATOS_ERROR_IF(list_id_closest_condition_of_correct_bc.size() != max_count) << "ERROR in list_id_closest_condition_of_correct_bc" << std::endl;

        // correct the projections
        countListClosestCondition = 0;
        for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it) {
            // int condId = listIdClosestCondition[countListClosestCondition];
            std::string rConditionName = rSkinModelPart.GetCondition(listIdClosestCondition[countListClosestCondition]).GetValue(CONDITION_NAME);

            auto gp_coord = (*it)->Center();
            int best_cond_id = -1;
            if (rConditionName != max_condition_name) 
            {
                double best_distance = 1e16;
                for (IndexType i = 0; i < max_count; i++)
                {
                    int condId = list_id_closest_condition_of_correct_bc[i];
                    auto cond_center = (&rSkinModelPart.GetCondition(condId))->GetGeometry().Center();
                    double curr_distance = norm_2(cond_center-gp_coord);

                    if (curr_distance < best_distance) 
                    {
                        best_distance = curr_distance;
                        best_cond_id = condId;
                    }
                }
                
                listIdClosestCondition[countListClosestCondition] = best_cond_id;
            }
            countListClosestCondition++;
        }
        
        if (max_condition_name == "SBMContact2DCondition") {
            ModelPart& analysis_model_part = rModelPart.GetParentModelPart();
            // KRATOS_WATCH(analysis_model_part)
            countListClosestCondition = 0;
            for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it) {
                int condId = listIdClosestCondition[countListClosestCondition];
                Condition::Pointer cond = rSkinModelPart.pGetCondition(condId);


                NodePointerVector empty_vector;
                empty_vector.push_back(cond->GetGeometry()(0)); // Just it_node-plane neighbours
                (*it)->SetValue(NEIGHBOUR_NODES, empty_vector);

                countListClosestCondition++;
            }

            ModelPart& layerModelPart = analysis_model_part.HasSubModelPart(layer_name) ? 
                            analysis_model_part.GetSubModelPart(layer_name) : 
                            analysis_model_part.CreateSubModelPart(layer_name);

            // retrieve the id of the brep
            IndexType id = (*rGeometriesBegin)->GetGeometryParent(0).Id();

            // add the brep to the SubModelPart
            layerModelPart.AddGeometry(analysis_model_part.pGetGeometry(id));

            return;
        }

        //--------------------------
        const Condition& rReferenceCondition = KratosComponents<Condition>::Get(max_condition_name);
        countListClosestCondition = 0;
        // in the Contact case the Condition is not created here

        ModelPart& layerModelPart = rModelPart.HasSubModelPart(layer_name) ? 
                            rModelPart.GetSubModelPart(layer_name) : 
                            rModelPart.CreateSubModelPart(layer_name);
        if (is2D) {
            for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it) {
                new_condition_list.push_back(
                    rReferenceCondition.Create(rIdCounter, (*it), pProperties));

                int condId = listIdClosestCondition[countListClosestCondition];

                Condition::Pointer cond1 = &rSkinModelPart.GetCondition(condId);
                int condId2;  
                if (condId == rSkinModelPart.ConditionsBegin()->Id()) {
                    condId2 = (rSkinModelPart.ConditionsEnd()-1)->Id();
                }
                else condId2 = condId-1;
                Condition::Pointer cond2 = &rSkinModelPart.GetCondition(condId2);

                // Add closest projection node
                NodePointerVector empty_vector;
                PointTypePointer projection_node = cond1->GetGeometry()(0);
                if (norm_2(projection_node->GetValue(NORMAL)) > 1e-13)
                {
                    empty_vector.push_back(cond1->GetGeometry()(0)); // Just it_node-plane neighbours
                    (*it)->SetValue(NEIGHBOUR_NODES, empty_vector);
                }
                //---------------------------------------------

                new_condition_list.GetContainer()[countListClosestCondition]->SetValue(NEIGHBOUR_CONDITIONS, GlobalPointersVector<Condition>({cond1,cond2}));
                if (isInner) {
                    new_condition_list.GetContainer()[countListClosestCondition]->SetValue(IDENTIFIER, "inner");
                } else {
                    new_condition_list.GetContainer()[countListClosestCondition]->SetValue(IDENTIFIER, "outer");
                }
                
                new_condition_list.GetContainer()[countListClosestCondition]->SetValue(MARKER_MESHES, mesh_size);
                                
                for (SizeType i = 0; i < (*it)->size(); ++i) {
                    // These are the control points associated with the basis functions involved in the condition we are creating
                    // rModelPart.AddNode((*it)->pGetPoint(i));
                    rModelPart.Nodes().push_back((*it)->pGetPoint(i));
                }
                rIdCounter++;
                countListClosestCondition++;
            }
        } else {
            // // 3D
            // for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it) {
            //     new_condition_list.push_back(
            //         rReferenceCondition.Create(rIdCounter, (*it), pProperties));

            //     int condId = listIdClosestCondition[countListClosestCondition];
            //     Condition::Pointer cond = &rSkinModelPart.GetCondition(condId);
            //     new_condition_list.GetContainer()[countListClosestCondition]->SetValue(NEIGHBOUR_CONDITIONS, GlobalPointersVector<Condition>({cond}));
            //     if (isInner) {
            //         new_condition_list.GetContainer()[countListClosestCondition]->SetValue(IDENTIFIER, "inner");
            //     } else {
            //         new_condition_list.GetContainer()[countListClosestCondition]->SetValue(IDENTIFIER, "outer");
            //     }
            //     for (SizeType i = 0; i < (*it)->size(); ++i) {
            //         // These are the control points associated with the basis functions involved in the condition we are creating
            //         // rModelPart.AddNode((*it)->pGetPoint(i));
            //         rModelPart.Nodes().push_back((*it)->pGetPoint(i));
            //     }
            //     rIdCounter++;
            //     countListClosestCondition++;
            // }
        }
        layerModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());
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
                // rModelPart.AddNode((*it)->pGetPoint(i));
                rModelPart.Nodes().push_back((*it)->pGetPoint(i));
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
                KRATOS_DEBUG_ERROR_IF(geom.LocalSpaceDimension() != 1) << "Geometry #" << geom.Id()
                    << " needs to have a dimension of 1 for type GeometryCurveNodes. LocalSpaceDimension: " << geom.LocalSpaceDimension()
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
                KRATOS_DEBUG_ERROR_IF(geom.LocalSpaceDimension() != 1) << "Geometry #" << geom.Id()
                    << " needs to have a dimension of 1 for type GeometryCurveVariationNodes. LocalSpaceDimension: " << geom.LocalSpaceDimension()
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
                KRATOS_DEBUG_ERROR_IF(geom.LocalSpaceDimension() != 2) << "Geometry #" << geom.Id()
                    << " needs to have a dimension of 2 for type " << rGeometryType << ". LocalSpaceDimension: "
                    << geom.LocalSpaceDimension() << ". Geometry" << geom << std::endl;

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

    bool IgaModeler::CheckIsOnExternalParameterSpace(Point point, Vector parameters_external_coordinates) const 
    {
        double tolerance = 1e-12;
        double u_initial = parameters_external_coordinates[0];
        double v_initial = parameters_external_coordinates[1];
        double u_final   = parameters_external_coordinates[2];
        double v_final   = parameters_external_coordinates[3];

        if (parameters_external_coordinates.size() > 4) {
            double w_initial = parameters_external_coordinates[4];
            double w_final   = parameters_external_coordinates[5];

            if (std::abs(point[0]-w_initial) < tolerance || std::abs(point[1]-w_final) < tolerance) {
                return true;
            }
        }
        
        if (std::abs(point[0]-u_initial) < tolerance || std::abs(point[1]-v_initial) < tolerance || 
            std::abs(point[0]-u_final  ) < tolerance || std::abs(point[1]-v_final  ) < tolerance) 
        {
            return true;
        }
        else {
            return false;
        }
    
    }

    void IgaModeler::prepareIntegrationOnTrueBoundary(ModelPart& analysis_model_part) const {
        // create the test bins containing all the boundary integration point (ALL!)
        
        PointVector points;
        std::string skin_model_part_name;
        if (!mParameters.Has("skin_model_part_name")) skin_model_part_name = "skin_model_part";
        else {
            skin_model_part_name = mParameters["skin_model_part_name"].GetString();
        }
        ModelPart& skin_model_part_in  = mpModel->GetModelPart(skin_model_part_name + "_in");
        ModelPart& skin_model_part_out = mpModel->GetModelPart(skin_model_part_name + "_out");

   
        

        if (!(skin_model_part_in.Nodes().size() > 0 || skin_model_part_out.Nodes().size() > 0) ) {
            KRATOS_ERROR << "Trying to integrate on true boundary when no skin boundary is defined" << std::endl;
        }

        for (auto i_cond : analysis_model_part.Conditions()) {
            points.push_back(PointTypePointer(new PointType(i_cond.Id(), i_cond.GetGeometry().Center().X(), i_cond.GetGeometry().Center().Y(), i_cond.GetGeometry().Center().Z())));
        }

        int order = 5;
        if (mParameters.Has("precision_order_on_integration")) {
            order = mParameters["precision_order_on_integration"].GetInt();
        } 
        int num_points = order+1;

        DynamicBins testBins(points.begin(), points.end());

        std::string surrogate_model_part_name;
        if (!mParameters.Has("surrogate_model_part_name")) surrogate_model_part_name = "surrogate_model_part";
        else {
            surrogate_model_part_name = mParameters["surrogate_model_part_name"].GetString();
        }

        Vector meshSizes_uv_inner = analysis_model_part.GetSubModelPart("surrogate_inner").GetValue(MARKER_MESHES);
        Vector meshSizes_uv_outer = analysis_model_part.GetSubModelPart("surrogate_outer").GetValue(MARKER_MESHES);
        
        // KRATOS_WATCH(meshSizes_uv_inner)
        // KRATOS_WATCH(meshSizes_uv_outer)

        const double meshSize= std::max(std::max(std::max(meshSizes_uv_inner[0], meshSizes_uv_inner[1]), meshSizes_uv_outer[0]), meshSizes_uv_outer[1]);
        const double radius = 2*sqrt(2)*(meshSize);

        const int numberOfResults = 1e6; 
        ModelPart::NodesContainerType::ContainerType Results(numberOfResults);
        std::vector<double> list_of_distances(numberOfResults);
        
        if (skin_model_part_in.Nodes().size() > 0) {

            const std::vector<std::array<double, 2>>& integration_point_list_u = IntegrationPointUtilities::s_gauss_legendre[num_points - 1];

            for (auto &i_cond : skin_model_part_in.Conditions()) {
                    // First and second point of the condition
                    const CoordinateVector U_0 = i_cond.GetGeometry()[0]; 
                    const CoordinateVector U_1 = i_cond.GetGeometry()[1];
                    CoordinateVector distance_u = U_1 - U_0;
                    const double length_u = norm_2(distance_u);
                    // Compue the integration points on this true segment
                    for (SizeType u = 0; u < num_points; ++u)
                    {
                        const CoordinateVector curr_integration_point = U_0 + distance_u * integration_point_list_u[u][0];
                        
                        const double curr_integration_weight = integration_point_list_u[u][1] * length_u;

                        PointerType pointToSearch = PointerType(new PointType(10000, curr_integration_point));

                        int obtainedResults = testBins.SearchInRadius(*pointToSearch, radius, Results.begin(), list_of_distances.begin(), numberOfResults);
                    
                    
                        double minimum_distance=1e10;
                        int nearestNodeId;
                        for (int i_distance = 0; i_distance < obtainedResults; i_distance++) {
                            double new_distance = list_of_distances[i_distance];   
                            if (new_distance < minimum_distance) { 
                                minimum_distance = new_distance;
                                nearestNodeId = i_distance;
                                }
                        }
                        if (obtainedResults == 0) {
                             KRATOS_WATCH("0 POINTS FOUND: EXIT")
                             KRATOS_WATCH(pointToSearch)
                             exit(0);}

                        
                        IndexType idCond = Results[nearestNodeId]->Id();
                        
                        std::vector<Vector> integration_point_list = analysis_model_part.GetCondition(idCond).GetValue(INTEGRATION_POINTS);
                        integration_point_list.push_back(curr_integration_point) ;
                        analysis_model_part.GetCondition(idCond).SetValue(INTEGRATION_POINTS, integration_point_list);

                        std::vector<double> integration_weight_list = analysis_model_part.GetCondition(idCond).GetValue(INTEGRATION_WEIGHTS);
                        integration_weight_list.push_back(curr_integration_weight) ;
                        analysis_model_part.GetCondition(idCond).SetValue(INTEGRATION_WEIGHTS, integration_weight_list);
                    }
            }

            



            // auto [points, weights] = boost::math::quadrature::gauss<double, num_points>::quadrature_points_and_weights();
            // for (auto &i_cond : skin_model_part_in.Conditions()) {

                    

            //         // Trasforma i punti e i pesi per l'intervallo [A, B]
            //         std::vector<double> gauss_points(num_points);
            //         for (int i = 0; i < num_points; ++i) {
            //             gauss_points[i] = 0.5 * ((B - A) * points[i] + (A + B));
            //         }

     
            //     }
            // }
        }
        

        

        // if (name.substr(0, 3) == "SBM") {
        //     is_inner = rParameters["is_inner"].GetBool();
        //     double meshSize;
        //     if (is_inner) { // INNER
        //         for (auto &i_cond : skin_model_part_in.Conditions()) {
        //             points.push_back(PointTypePointer(new PointType(i_cond.Id(), i_cond.GetGeometry().Center().X(), i_cond.GetGeometry().Center().Y(), i_cond.GetGeometry().Center().Z())));
        //         }
        //         meshSizes_uv = mpModel->GetModelPart(surrogate_model_part_name + "_inner").GetProcessInfo().GetValue(MARKER_MESHES);
        //         parameterExternalCoordinates = mpModel->GetModelPart(surrogate_model_part_name + "_inner").GetProcessInfo().GetValue(LOAD_MESHES);
        //         meshSize = meshSizes_uv[0];
        //         if (meshSizes_uv[1] > meshSize) {meshSize = meshSizes_uv[1];}
        //         if (meshSizes_uv.size() > 2) {if (meshSizes_uv[2] > meshSize) {meshSize = meshSizes_uv[2];}}
        //     } 
        //     else { // OUTER
        //         for (auto &i_cond : skin_model_part_out.Conditions()) {
        //             points.push_back(PointTypePointer(new PointType(i_cond.Id(), i_cond.GetGeometry().Center().X(), i_cond.GetGeometry().Center().Y(), i_cond.GetGeometry().Center().Z())));
        //         }

        //         meshSizes_uv = mpModel->GetModelPart(surrogate_model_part_name + "_outer").GetProcessInfo().GetValue(MARKER_MESHES);
        //         parameterExternalCoordinates = mpModel->GetModelPart(surrogate_model_part_name + "_outer").GetProcessInfo().GetValue(LOAD_MESHES);
        //         meshSize = meshSizes_uv[0];
        //         if (meshSizes_uv[1] > meshSize) {meshSize = meshSizes_uv[1];}
        //         if (meshSizes_uv.size() > 2) {if (meshSizes_uv[2] > meshSize) {meshSize = meshSizes_uv[2];}}
        //     }
            
        //     is_SBM = true;

        //     radius = sqrt(3)*(meshSize); 
        //     // radius = 30*(meshSize);
        // }
    }

    ///@}
}
