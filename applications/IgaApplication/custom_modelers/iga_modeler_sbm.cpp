//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolo' Antonelli
//                   Andrea Gorgi
//


// Project includes
#include "iga_modeler_sbm.h"
#include "integration/integration_point_utilities.h"
#include "iga_application_variables.h"


namespace Kratos
{
///@name Stages
///@{

void IgaModelerSbm::SetupModelPart()
{
    KRATOS_ERROR_IF_NOT(mParameters.Has("analysis_model_part_name"))
        << "Missing \"analysis_model_part_name\" section" << std::endl;
    ModelPart& analysis_model_part = mpModel->GetModelPart(mParameters["analysis_model_part_name"].GetString());

    KRATOS_ERROR_IF_NOT(mParameters.Has("element_condition_list"))
        << "Missing \"element_condition_list\" section" << std::endl;

    const Parameters iga_physics_parameters = mParameters["element_condition_list"];

    CreateIntegrationDomain(
        analysis_model_part,
        iga_physics_parameters);

}

///@}

void IgaModelerSbm::CreateIntegrationDomain(
    ModelPart& rModelPart,
    const Parameters rPhysicsParameters) const
{

    KRATOS_ERROR_IF_NOT(rPhysicsParameters.IsArray())
        << "\"element_condition_list\" needs to be an array." << std::endl;

    for (SizeType i = 0; i < rPhysicsParameters.size(); ++i)
    {
        CreateIntegrationDomainPerUnit(
            rModelPart,
            rPhysicsParameters[i]);
    }
}

void IgaModelerSbm::CreateIntegrationDomainPerUnit(
    ModelPart& rModelPart,
    const Parameters rPhysicsParameters) const
{
    KRATOS_ERROR_IF_NOT(rPhysicsParameters.Has("iga_model_part"))
        << "\"iga_model_part\" need to be specified." << std::endl;

    std::string sub_model_part_name = rPhysicsParameters["iga_model_part"].GetString();

    ModelPart& sub_model_part = rModelPart.HasSubModelPart(sub_model_part_name)
        ? rModelPart.GetSubModelPart(sub_model_part_name)
        : rModelPart.CreateSubModelPart(sub_model_part_name);

    // Generate the list of geometries, which are needed, here.
    GeometriesArrayType geometry_list;
    GetGeometryList(geometry_list, rModelPart, rPhysicsParameters);
    
    KRATOS_ERROR_IF_NOT(rPhysicsParameters.Has("geometry_type")) << 
         "::[IgaModelerSbm]:: Missing \"geometry_type\" parameter." << rPhysicsParameters << std::endl;
                
    std::string geometry_type = rPhysicsParameters["geometry_type"].GetString();

    if (!rPhysicsParameters.Has("sbm_parameters"))
        CreateQuadraturePointGeometries(
            geometry_list, sub_model_part, rPhysicsParameters, geometry_type); 
    else 
        // it must be a condition
        CreateQuadraturePointGeometriesSbm(
            geometry_list, sub_model_part, rPhysicsParameters, geometry_type);

    KRATOS_INFO_IF("CreateIntegrationDomainElementCondition", mEchoLevel > 3)
        << "Creation of elements/ conditions finished in: " << sub_model_part << std::endl;
}




///@}
///@name CAD functionalities
///@{

void IgaModelerSbm::GetGeometryList(
    GeometriesArrayType& rGeometryList,
    ModelPart& rModelPart,
    const Parameters rPhysicsParameters) const
{
    /* we have three cases:
        1) background geometry: i.e. the surface or volume. they are type = "element"
        2) outer_loop -> is_inner = false : type = "condition" (2D brep curves & 3D brep surfaces)
        2) inner_loop -> is_inner = true  : type = "condition" (2D brep curves & 3D brep surfaces)
    */ 

    const std::string type = rPhysicsParameters["type"].GetString();

    if (type == "element")
    {
        int surface_brep_id = 1; 
        rGeometryList.push_back(rModelPart.pGetGeometry(surface_brep_id));
    } 
    else if (type == "condition")
    {
        if (rPhysicsParameters.Has("sbm_parameters"))
        {
            KRATOS_ERROR_IF_NOT(rPhysicsParameters["sbm_parameters"].Has("is_inner")) << 
                "::[IgaModelerSbm]:: Missing \"is_inner\" parameter in the \"sbm_parameters\"." << rPhysicsParameters["sbm_parameters"] 
                << std::endl;
                        
            if (rPhysicsParameters["sbm_parameters"]["is_inner"].GetBool()) // inner loop
            {
                // Surface id is 1
                int inner_brep_id = 2;
                ModelPart& surrogate_model_part_outer = rModelPart.GetSubModelPart("surrogate_outer");
                if (surrogate_model_part_outer.NumberOfConditions() == 0)
                    // 2D case
                    inner_brep_id += 4; // if there is no outer we use the 4 sides of the rectangle
                    // 3D -> //TODO:
                else 
                    // if outer loop is present take the number of conditions
                    inner_brep_id += surrogate_model_part_outer.NumberOfConditions();

                // INNER   
                ModelPart& surrogate_model_part_inner = rModelPart.GetSubModelPart("surrogate_inner");

                KRATOS_ERROR_IF(surrogate_model_part_inner.NumberOfElements() == 0) 
                    << "::[IgaModelerSbm]:: The surrogate_model_part_inner has zero elements (no inner loop/boundary defined)."
                    << "Something might be missing in the NurbsModelerSbm." << std::endl;

                for (IndexType iel = 1; iel < surrogate_model_part_inner.NumberOfElements()+1; iel++) { //loop over the number of inner loops
                    /*
                    Each element in the surrogate_model_part_inner represents a surrogate boundary loop. First "node.Id()" is the id of the first condition and
                        the second "node.Id()" is the last condition of that loop. (Essential for multiple inner loops)
                    */
                    IndexType first_condition_id = surrogate_model_part_inner.pGetElement(iel)->GetGeometry()[0].Id();
                    IndexType last_condition_id = surrogate_model_part_inner.pGetElement(iel)->GetGeometry()[1].Id();

                    SizeType size_surrogate_loop = last_condition_id - first_condition_id + 1;

                    for (SizeType j = 0; j < size_surrogate_loop; ++j) {
                        rGeometryList.push_back(rModelPart.pGetGeometry(inner_brep_id));
                        inner_brep_id++;
                    }
                }
            } else // outer loop
            {
                int outer_brep_id = 2;
                // OUTER
                ModelPart& surrogate_model_part_outer = rModelPart.GetSubModelPart("surrogate_outer");
                
                if (surrogate_model_part_outer.NumberOfConditions() > 0) {
                    // 2D
                    if ((surrogate_model_part_outer.ConditionsBegin())->GetGeometry().size() == 2) {
                        const int size_surrogate_loop_outer = surrogate_model_part_outer.NumberOfConditions();
                        for (int j = 0; j < size_surrogate_loop_outer; ++j) {
                            rGeometryList.push_back(rModelPart.pGetGeometry(outer_brep_id));
                            outer_brep_id++;
                        }
                    } 
                }
            }
        }
        // else -> body-fitted case
        else {
            KRATOS_ERROR_IF_NOT(rPhysicsParameters.Has("brep_ids")) << 
                "::[IgaModelerSbm]:: Missing \"brep_ids\" parameter in body-fitted boundary condition\"." << rPhysicsParameters 
                << std::endl;
            if (rPhysicsParameters.Has("brep_ids")) {
                for (SizeType i = 0; i < rPhysicsParameters["brep_ids"].size(); ++i) {
                    rGeometryList.push_back(rModelPart.pGetGeometry(rPhysicsParameters["brep_ids"][i].GetInt()));
                }
            }
        }
    } else 
    {
        KRATOS_ERROR << "::[IgaModelerSbm]:: type " << type << " not defined in the IgaModelerSbm."
                     << " Available types are: elements, conditions." << std::endl;
    }

    KRATOS_ERROR_IF(rGeometryList.size() == 0)
        << "::[IgaModelerSbm]:: Empty geometry list in GetGeometryList. Physics parameter is: " << rPhysicsParameters << std::endl;
}




void IgaModelerSbm::CreateQuadraturePointGeometries(
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
        KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 3)
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
            std::vector<int> list_id_closest_condition(geometries.size());
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

void IgaModelerSbm::CreateQuadraturePointGeometriesSbm(
    GeometriesArrayType& rGeometryList,
    ModelPart& rModelPart,
    const Parameters rParameters,
    std::string GeometryType) const
{
    KRATOS_ERROR_IF_NOT(rParameters.Has("type"))
        << "\"type\" needs to be specified." << std::endl;
    KRATOS_ERROR_IF_NOT(rParameters.Has("name"))
        << "\"name\" needs to be specified." << std::endl;

    // Only conditions should call CreateQuadraturePointGeometriesSbm
    std::string type = rParameters["type"].GetString();
    bool check_input_type = (type == "condition" || type == "Condition");
    KRATOS_ERROR_IF_NOT(check_input_type) << ":::[IgaModelerSbm]::: type != \"condition\" in CreateQuadraturePointGeometriesSbm. "
                                          << "It must be a condition to apply the sbm operators. type: " << type << std::endl;
                                          
    std::string name = rParameters["name"].GetString();

    SizeType shape_function_derivatives_order = 1;
    if (rParameters.Has("shape_function_derivatives_order")) {
        shape_function_derivatives_order = rParameters["shape_function_derivatives_order"].GetInt();
    }
    else {
        KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 1)
            << "shape_function_derivatives_order is not provided and thus being considered as 1. " << std::endl;
    }

    std::string quadrature_method = rParameters.Has("quadrature_method")
        ? rParameters["integration_rule"].GetString()
        : "GAUSS";

    KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 0)
        << "Creating " << name << "s of type: " << type
        << " for " << rGeometryList.size() << " geometries"
        << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

    // Check if the sbm projection operation is needed (there is no need for background domain and for body-fitted boundary conditions)
    PointVector points;   
    const std::string skin_model_part_name = mParameters.Has("skin_model_part_name")
            ? mParameters["skin_model_part_name"].GetString()
            : "skin_model_part";

    ModelPart& skin_model_part = mpModel->HasModelPart(skin_model_part_name)
            ? mpModel->GetModelPart(skin_model_part_name)
            : KRATOS_ERROR << "::[CreateQuadraturePointGeometriesSbm]::: Sbm case -> skin_model_part has not been defined before. "
                            << "Maybe you are not calling the nurbs_modeler_sbm" << std::endl;

    // inner & outer are defaulf sub model part names
    ModelPart& skin_sub_model_part_in = skin_model_part.GetSubModelPart("inner");
    ModelPart& skin_sub_model_part_out = skin_model_part.GetSubModelPart("outer");

    const bool is_inner = rParameters["sbm_parameters"]["is_inner"].GetBool();
    
    const std::string surrogate_sub_model_part_name = is_inner ? "surrogate_inner" : "surrogate_outer";
    
    ModelPart& surrogate_sub_model_part = rModelPart.GetParentModelPart().HasSubModelPart(surrogate_sub_model_part_name)
            ? rModelPart.GetParentModelPart().GetSubModelPart(surrogate_sub_model_part_name)
            : KRATOS_ERROR << "::[CreateQuadraturePointGeometriesSbm]::: Sbm case -> surrogate_sub_model_part has not been defined before."
                            << "Maybe you are not calling the nurbs_modeler_sbm" << std::endl;

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
    
    // Get the mesh sizes from the surrogate model part
    const Vector knot_span_sizes = surrogate_sub_model_part.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);
    // Get the parameter space corners from the surrogate model part
    const std::vector<Vector> parameter_space_corners = surrogate_sub_model_part.GetParentModelPart().GetValue(PARAMETER_SPACE_CORNERS);

    double knot_span_reference_size = knot_span_sizes[0];
    if (knot_span_sizes[1] > knot_span_reference_size) {knot_span_reference_size = knot_span_sizes[1];}
    if (knot_span_sizes.size() > 2) {if (knot_span_sizes[2] > knot_span_reference_size) {knot_span_reference_size = knot_span_sizes[2];}}

    // 2D case
    const double search_radius = sqrt(2)*(knot_span_reference_size); 
    // 3D case //TODO:
    DynamicBins testBins(points.begin(), points.end());
    
    // Maximum number of results to be found in the search in radius
    const int number_of_results = 1e6; 

    ModelPart::NodesContainerType::ContainerType results(number_of_results);
    std::vector<double> list_of_distances(number_of_results);
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
    
        rGeometryList[i].CreateQuadraturePointGeometries(geometries, shape_function_derivatives_order, integration_info);

        KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 1)
            << geometries.size() << " quadrature point geometries have been created." << std::endl;

        SizeType id = 1;
        if (rModelPart.GetRootModelPart().NumberOfConditions() > 0)
            id = rModelPart.GetRootModelPart().Conditions().back().Id() + 1;
        
        std::vector<int> list_id_closest_condition(geometries.size());

        for (SizeType j= 0; j < geometries.size() ; j++) {  

            const Point integration_point = geometries[j].Center(); 
            PointerType p_integration_point = PointerType(new PointType(1, integration_point.X(), integration_point.Y(), integration_point.Z()));

            // Use the search in radius to find the closest point
            SizeType obtained_results = testBins.SearchInRadius(*p_integration_point, search_radius, results.begin(), list_of_distances.begin(), number_of_results);

            double minimum_distance = 1e14;

            // Find the nearest node
            IndexType nearest_node_id;
            for (IndexType k = 0; k < obtained_results; k++) {
                double current_distance = list_of_distances[k];   
                if (current_distance < minimum_distance) { 
                    minimum_distance = current_distance;
                    nearest_node_id = k;
                }
            }
            KRATOS_ERROR_IF(obtained_results == 0) << "::[IgaModelerSbm]:: Zero points found in serch for projection of point: " <<
                p_integration_point << std::endl;
            
            // store id closest condition
            list_id_closest_condition[j] = results[nearest_node_id]->Id();
        }

        if (is_inner) {
            // pass the skin_sub_model_part_in
            this->CreateConditions( geometries.ptr_begin(), geometries.ptr_end(),
                rModelPart, skin_sub_model_part_in, list_id_closest_condition, name, id, PropertiesPointerType(), is_inner, knot_span_sizes);
        }
        else{
            // pass the skin_sub_model_part_out
            this->CreateConditions(geometries.ptr_begin(), geometries.ptr_end(),
                rModelPart, skin_sub_model_part_out, list_id_closest_condition, name, id, PropertiesPointerType(), is_inner, knot_span_sizes);
        }
    }
}


///@}
///@name Generate Elements and Conditions
///@{

void IgaModelerSbm::CreateElements(
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

    SizeType num_elements = std::distance(rGeometriesBegin, rGeometriesEnd);
    new_element_list.reserve(num_elements);

    for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it)
    {
        new_element_list.push_back(
            rReferenceElement.Create(rIdCounter, (*it), pProperties));
        for (SizeType i = 0; i < (*it)->size(); ++i) {
            rModelPart.Nodes().push_back((*it)->pGetPoint(i));
        }
        rIdCounter++;
    }

    rModelPart.AddElements(new_element_list.begin(), new_element_list.end());
}


void IgaModelerSbm::CreateConditions(
    typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
    typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
    ModelPart& rModelPart,
    std::string& rConditionName,
    SizeType& rIdCounter,
    PropertiesPointerType pProperties) const
{
    const Condition& reference_condition = KratosComponents<Condition>::Get(rConditionName);

    ModelPart::ConditionsContainerType new_condition_list;

    KRATOS_INFO_IF("CreateConditions", mEchoLevel > 2)
        << "Creating conditions of type " << rConditionName
        << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

    for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it)
    {
        new_condition_list.push_back(
            reference_condition.Create(rIdCounter, (*it), pProperties));
        for (SizeType i = 0; i < (*it)->size(); ++i) {
            rModelPart.Nodes().push_back((*it)->pGetPoint(i));
        }
        rIdCounter++;
    }

    rModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());
}


void IgaModelerSbm::CreateConditions(
    typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
    typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
    ModelPart& rModelPart,
    ModelPart& rSkinModelPart,
    std::vector<int>& listIdClosestCondition,
    std::string& rConditionName,
    SizeType& rIdCounter,
    PropertiesPointerType pProperties,
    bool IsInner,
    Vector KnotSpanSizes) const
{
    const Condition& reference_condition = KratosComponents<Condition>::Get(rConditionName);

    ModelPart::ConditionsContainerType new_condition_list;

    KRATOS_INFO_IF("CreateConditions", mEchoLevel > 2)
        << "Creating conditions of type " << rConditionName
        << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

    int countListClosestCondition = 0;

    // 2D case
    if (rSkinModelPart.ConditionsBegin()->GetGeometry().size() == 2) {

        for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it) {
            new_condition_list.push_back(reference_condition.Create(rIdCounter, (*it), pProperties));

            IndexType condId = listIdClosestCondition[countListClosestCondition];

            Condition::Pointer cond1 = &rSkinModelPart.GetCondition(condId);

            IndexType condId2;
            if (condId == rSkinModelPart.ConditionsBegin()->Id()) condId2 = (rSkinModelPart.ConditionsEnd()-1)->Id();
            else condId2 = condId-1;

            Condition::Pointer cond2 = &rSkinModelPart.GetCondition(condId2);

            new_condition_list.GetContainer()[countListClosestCondition]->SetValue(NEIGHBOUR_CONDITIONS, GlobalPointersVector<Condition>({cond1,cond2}));
            if (IsInner) {
                new_condition_list.GetContainer()[countListClosestCondition]->SetValue(IDENTIFIER, "inner");
            } else {
                new_condition_list.GetContainer()[countListClosestCondition]->SetValue(IDENTIFIER, "outer");
            }
            new_condition_list.GetContainer()[countListClosestCondition]->SetValue(KNOT_SPAN_SIZES, KnotSpanSizes);
                        
            for (SizeType i = 0; i < (*it)->size(); ++i) {
                // These are the control points associated with the basis functions involved in the condition we are creating
                rModelPart.Nodes().push_back((*it)->pGetPoint(i));
            }
            rIdCounter++;
            countListClosestCondition++;
        }
    } else {
        // TODO: 3D case
        KRATOS_ERROR << "CreateConditions: 3D case not implemented yet." << std::endl;
    }
    
    rModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());
}

///@}
}