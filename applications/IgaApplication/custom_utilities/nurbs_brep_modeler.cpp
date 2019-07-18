//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Thomas Oberbichler
//

// System includes

// External includes

// Project includes
#include "nurbs_brep_modeler.h"


namespace Kratos
{
    void NurbsBrepModeler::ImportGeometry(BrepJsonIO& rBrepJsonIO, Parameters& rNurbsBrepGeometryJson)
    {
        std::cout << "ImportGeometry in of geometry" << std::endl;
        std::vector<BrepModel> brep_model_vector = rBrepJsonIO.ImportNurbsBrepGeometry(m_model_part, rNurbsBrepGeometryJson);
        for (auto brep_model = brep_model_vector.begin(); brep_model != brep_model_vector.end(); ++brep_model)
        {
            m_brep_model_vector.push_back(*brep_model);
        }
    }

    const BrepEdge& NurbsBrepModeler::GetBrepEdge(int& rBrepId) const
    {
        for (int i = 0; i < m_brep_model_vector.size(); ++i)
        {
            const std::vector<BrepEdge>& edge_vector = m_brep_model_vector[i].GetEdgeVector();

            for (int j = 0; j < edge_vector.size(); ++j)
            {
                if (edge_vector[j].Id() == rBrepId)
                    return edge_vector[j];
            }
        }
        KRATOS_ERROR << "Brep Id: " << rBrepId << " is not of type edge." << std::endl;
    }

    const BrepFace& NurbsBrepModeler::GetBrepFace(int& rBrepId) const
    {
        for (int i = 0; i < m_brep_model_vector.size(); ++i)
        {
            const std::vector<BrepFace>& face_vector = m_brep_model_vector[i].GetFaceVector();

            for (int j = 0; j < face_vector.size(); ++j)
            {
                if (face_vector[j].Id() == rBrepId)
                    return face_vector[j];
            }
        }
        KRATOS_ERROR << "Brep Id: " << rBrepId << " is not of type face." << std::endl;
    }

    const BrepVertex& NurbsBrepModeler::GetBrepVertex(int& rBrepId) const
    {
        for (int i = 0; i < m_brep_model_vector.size(); ++i)
        {
            const std::vector<BrepVertex>& vertex_vector = m_brep_model_vector[i].GetVertexVector();

            for (int j = 0; j < vertex_vector.size(); ++j)
            {
                if (vertex_vector[j].Id() == rBrepId)
                    return vertex_vector[j];
            }
        }
        KRATOS_ERROR << "Brep Id: " << rBrepId << " is not of type vertex." << std::endl;
    }

    void NurbsBrepModeler::ImportModelPart(ModelPart& model_part, Parameters& rModelPartParameters)
    {
        for (int i = 0; i < rModelPartParameters["element_condition_list"].size(); ++i)
        {
            Parameters element_parameter = rModelPartParameters["element_condition_list"][i];

            std::string sub_model_part_name = element_parameter["iga_model_part"].GetString();
            ModelPart& sub_model_part = model_part.HasSubModelPart(sub_model_part_name)
                ? model_part.GetSubModelPart(sub_model_part_name)
                : model_part.CreateSubModelPart(sub_model_part_name);

            std::string geometry_type = element_parameter["geometry_type"].GetString();

            if (geometry_type == "GeometryCurveNodes")
            {
                for (int j = 0; j < element_parameter["brep_ids"].size(); ++j)
                {
                    int brep_id = element_parameter["brep_ids"][j].GetInt();
                    Vector parameter = element_parameter["parameters"]["local_parameters"].GetVector();
                    GetBrepEdge(brep_id).GetGeometryNodes(sub_model_part, parameter[0]);
                }
            }
            else if (geometry_type == "GeometryCurveVariationNodes")
            {
                for (int j = 0; j < element_parameter["brep_ids"].size(); ++j)
                {
                    int brep_id = element_parameter["brep_ids"][j].GetInt();
                    Vector parameter = element_parameter["parameters"]["local_parameters"].GetVector();
                    GetBrepEdge(brep_id).GetGeometryVariationNodes(sub_model_part, parameter[0]);
                }
            }
            if (geometry_type == "GeometrySurfaceNodes")
            {
                for (int j = 0; j < element_parameter["brep_ids"].size(); ++j)
                {
                    int brep_id = element_parameter["brep_ids"][j].GetInt();
                    Vector parameter = element_parameter["parameters"]["local_parameters"].GetVector();
                    GetBrepFace(brep_id).GetGeometryNodes(sub_model_part, parameter[0], parameter[1]);
                }
            }
            else if (geometry_type == "GeometrySurfaceVariationNodes")
            {
                for (int j = 0; j < element_parameter["brep_ids"].size(); ++j)
                {
                    int brep_id = element_parameter["brep_ids"][j].GetInt();
                    Vector parameter = element_parameter["parameters"]["local_parameters"].GetVector();
                    GetBrepFace(brep_id).GetGeometryVariationNodes(sub_model_part, parameter[0], parameter[1]);
                }
            }
            else
            {
                std::string type = element_parameter["parameters"]["type"].GetString();
                std::string name = element_parameter["parameters"]["name"].GetString();
                int shape_function_derivatives_order = element_parameter["parameters"]["shape_function_derivatives_order"].GetInt();

                std::vector<std::string> variable_list;
                for (int j = 0; j < element_parameter["parameters"]["variables"].size(); ++j)
                {
                    variable_list.push_back(element_parameter["parameters"]["variables"][j].GetString());
                }

                //std::vector<shared_ptr<Element>> element_vector;

                if (element_parameter.Has("brep_ids"))
                {
                    for (int j = 0; j < element_parameter["brep_ids"].size(); ++j)
                    {
                        int brep_id = element_parameter["brep_ids"][j].GetInt();
                        if (geometry_type == "GeometrySurface")
                        {
                            auto brep_face = GetBrepFace(brep_id);

                            auto element_vector = IgaIntegrationUtilities::GetIntegrationDomainGeometrySurface(
                                brep_face.GetSurface(),
                                brep_face.GetSurfaceClipper(0.001, 0.00001),
                                shape_function_derivatives_order);

                            for(auto element = element_vector.begin(); element != element_vector.end(); ++element)
                                (*element)->SetValue(BREP_ID, brep_id);

                            if (type == "element")
                            {
                                int id = 1;
                                if (sub_model_part.GetRootModelPart().Elements().size() > 0)
                                    id = sub_model_part.GetRootModelPart().Elements().back().Id() + 1;

                                IgaIntegrationUtilities::ChangeElementType(
                                    element_vector, sub_model_part, name, id);
                            }

                            if (type == "condition")
                            {
                                int id = 0;
                                if (sub_model_part.GetRootModelPart().Conditions().size() > 0)
                                    id = sub_model_part.GetRootModelPart().Conditions().back().Id() + 1;

                                IgaIntegrationUtilities::ChangeConditionType(
                                    element_vector, sub_model_part, name, id);
                            }
                        }
                        if (geometry_type == "SurfaceEdge")
                        {
                            const auto edge = GetBrepEdge(brep_id);
                            auto topology = edge.GetEdgeTopology(0);

                            const auto brep_face = GetBrepFace(topology.brep_id);
                            const auto surface = brep_face.GetSurface();
                            const auto trimming_curve = brep_face.GetTrimCurve(topology.trim_index);

                            auto geometry_vector = IgaIntegrationUtilities::GetIntegrationDomainSurfaceEdge(
                                surface,
                                trimming_curve,
                                shape_function_derivatives_order);

                            if (type == "condition")
                            {
                                int id = 1;
                                if (sub_model_part.GetRootModelPart().Conditions().size() > 0)
                                    id = sub_model_part.GetRootModelPart().Conditions().back().Id() + 1;

                                IgaIntegrationUtilities::CreateConditions(
                                    geometry_vector, sub_model_part, name, id);
                            }
                        }
                        if (geometry_type == "SurfacePoint")
                        {
                            GetBrepEdge(brep_id);
                            bool success = m_brep_model_vector[0].GetIntegrationDomainBrep(
                                sub_model_part, brep_id, type, name,
                                shape_function_derivatives_order, variable_list);
                        }
                        else if (geometry_type == "GeometryCurve")
                        {
                            GetBrepEdge(brep_id).GetIntegrationGeometry(
                                sub_model_part, type, name,
                                shape_function_derivatives_order,
                                variable_list);
                        }
                        else if (geometry_type == "CurveEdge")
                        {
                            GetBrepEdge(brep_id).GetIntegrationGeometry(
                                sub_model_part, type, name,
                                shape_function_derivatives_order,
                                variable_list);
                        }
                        else if (geometry_type == "CurvePoint")
                        {
                            GetBrepEdge(brep_id).GetIntegrationGeometry(
                                sub_model_part, type, name,
                                shape_function_derivatives_order,
                                variable_list);
                        }
                        if (geometry_type == "SurfaceEdgeSurfaceEdge")
                        {
                            const auto coupling_edge = GetBrepEdge(brep_id);
                            auto master_topology = coupling_edge.GetEdgeTopology(0);

                            const auto brep_face_1 = GetBrepFace(master_topology.brep_id);
                            const auto surface_1 = brep_face_1.GetSurface();
                            const auto trimming_curve_1 = brep_face_1.GetTrimCurve(master_topology.trim_index);

                            int number_of_couplings = coupling_edge.GetNumberOfEdgeTopologies();
                            if (number_of_couplings > 1)
                            {
                                for (int i = 1; i < number_of_couplings; ++i)
                                {
                                    auto slave_topology = coupling_edge.GetEdgeTopology(i);

                                    const auto brep_face_2 = GetBrepFace(slave_topology.brep_id);
                                    const auto surface_2 = brep_face_2.GetSurface();
                                    const auto trimming_curve_2 = brep_face_2.GetTrimCurve(slave_topology.trim_index);

                                    auto geometry_vector = IgaIntegrationUtilities::GetIntegrationDomainSurfaceEdgeSurfaceEdge(
                                        surface_1,
                                        trimming_curve_1,
                                        surface_2,
                                        trimming_curve_2,
                                        shape_function_derivatives_order);

                                    if (type == "element")
                                    {
                                        int id = 1;
                                        if (sub_model_part.GetRootModelPart().Elements().size() > 0)
                                            id = sub_model_part.GetRootModelPart().Elements().back().Id() + 1;

                                        //IgaIntegrationUtilities::ChangeElementType(
                                        //    element_vector, sub_model_part, name, id);
                                    }

                                    if (type == "condition")
                                    {
                                        int id = 1;
                                        if (sub_model_part.GetRootModelPart().Conditions().size() > 0)
                                            id = sub_model_part.GetRootModelPart().Conditions().back().Id() + 1;

                                        IgaIntegrationUtilities::CreateConditions(
                                            geometry_vector, sub_model_part, name, id);
                                    }
                                }
                            }
                        }
                        if (geometry_type == "SurfaceEdgeCurveEdge")
                        {
                            GetBrepEdge(brep_id);
                            bool success = m_brep_model_vector[0].GetIntegrationDomainBrepCoupling(
                                sub_model_part, brep_id, type, name,
                                shape_function_derivatives_order, variable_list);
                        }
                        if (geometry_type == "SurfacePointSurfacePoint")
                        {
                            GetBrepEdge(brep_id);
                            bool success = m_brep_model_vector[0].GetIntegrationDomainBrepCoupling(
                                sub_model_part, brep_id, type, name,
                                shape_function_derivatives_order, variable_list);
                        }
                        if (geometry_type == "SurfacePointCurvePoint")
                        {
                            GetBrepEdge(brep_id);
                            bool success = m_brep_model_vector[0].GetIntegrationDomainBrepCoupling(
                                sub_model_part, brep_id, type, name,
                                shape_function_derivatives_order, variable_list);
                        }
                        if (geometry_type == "CurveEdgeCurveEdge")
                        {
                            GetBrepEdge(brep_id);
                            bool success = m_brep_model_vector[0].GetIntegrationDomainBrepCoupling(
                                sub_model_part, brep_id, type, name,
                                shape_function_derivatives_order, variable_list);
                        }
                        if (geometry_type == "CurvePointCurvePoint")
                        {
                            GetBrepEdge(brep_id);
                            bool success = m_brep_model_vector[0].GetIntegrationDomainBrepCoupling(
                                sub_model_part, brep_id, type, name,
                                shape_function_derivatives_order, variable_list);
                        }
                    }
                }
/*
                KRATOS_WATCH(element_vector.size())

                if (type == "element")
                {
                    int id = 1;
                    if (sub_model_part.GetRootModelPart().Elements().size() > 0)
                        id = sub_model_part.GetRootModelPart().Elements().back().Id() + 1;

                    IgaIntegrationUtilities::ChangeElementType(
                        element_vector, sub_model_part, name, id);
                }

                if (type == "condition")
                {
                    int id = 0;
                    if (sub_model_part.GetRootModelPart().Conditions().size() > 0)
                        id = sub_model_part.GetRootModelPart().Conditions().back().Id() + 1;

                    IgaIntegrationUtilities::ChangeConditionType(
                        element_vector, sub_model_part, name, id);
                }

*/
                KRATOS_WATCH(model_part)
            }
        }
    }

    void NurbsBrepModeler::GetInterfaceConditions(
        ModelPart& rExternalModelPart,
        ModelPart& rIgaModelPart,
        ModelPart& rInterfaceConditionsModelPart,
        const std::string& rConditionName)
    {
        //const Condition& ReferenceCondition = KratosComponents<Condition>::Get(rConditionName);
        //ModelPart::ConditionsContainerType new_condition_list;

        //BinsIgaConfigure::ContainerType BinsIgaObjects;
        //BinsIgaObjects.resize(rIgaModelPart.NumberOfElements());
        //const auto elements_begin = rIgaModelPart.Elements().ptr_begin();

        //for (int i = 0; i < rIgaModelPart.NumberOfElements(); ++i)
        //{
        //    auto it_elem = elements_begin + i;
        //    KRATOS_WATCH(*it_elem)
        //    BinsIgaObjects[i] = (Kratos::make_shared<BinsIgaObject>(*it_elem));
        //}
        //auto iga_bins_structure = BinsObjectDynamic<BinsIgaConfigure>(BinsIgaObjects.begin(), BinsIgaObjects.end());

        //int num_interface_obj_bin = BinsIgaObjects.size();

        //BinsIgaConfigure::ResultContainerType neighbor_results(num_interface_obj_bin);
        //std::vector<double> neighbor_distances(num_interface_obj_bin);

        ////Instance of external interface nodes. Updated for each node.
        //auto external_object = Kratos::make_shared<BinsIgaObject>(array_1d<double, 3>(0.0));

        //for (auto external_node_ptr = rExternalModelPart.NodesBegin();
        //    external_node_ptr != rExternalModelPart.NodesEnd();
        //    ++external_node_ptr)
        //{
        //    auto results_itr = neighbor_results.begin();
        //    auto distance_itr = neighbor_distances.begin();

        //    //update coordinates of extern
        //    external_object->UpdateCoordinates(
        //        external_node_ptr->Coordinates());

        //    const std::size_t number_of_results = iga_bins_structure.SearchObjectsInRadius(
        //        external_object,
        //        rSearchRadius,
        //        results_itr,
        //        distance_itr,
        //        num_interface_obj_bin);

        //    // todo: neglect bad solutions

        //    if (number_of_results > 0)
        //    {
        //        Element::Pointer p_closest_element;
        //        double distance_to_node = 1e10;
        //        for (int i = 0; i < number_of_results; ++i)
        //        {
        //            int brep_id = neighbor_results[i]->pGetBaseElement()->GetValue(BREP_ID);

        //            auto surface = GetBrepFace(brep_id).GetSurface();

        //            Element::Pointer p_projected_element =  IgaIntegrationUtilities::ProjectNodeOnSurface(
        //                surface,
        //                neighbor_results[i]->pGetBaseElement()->GetValue(LOCAL_COORDINATES),
        //                *external_node_ptr,
        //                rShapeFunctionDerivativesOrder,
        //                rAccuray,
        //                rNumberOfIterations);

        //            double new_distance = IgaGeometryUtilities::CalculateDistanceNodeElement(
        //                p_projected_element,
        //                *external_node_ptr);

        //            if (new_distance < distance_to_node)
        //            {
        //                p_closest_element = p_projected_element;
        //                distance_to_node = new_distance;
        //            }
        //        }


        //        //Create new tolerances
        //        if (rTolerance > distance_to_node)
        //        {
        //            auto p_condition = ReferenceCondition.Create(
        //                rIdCounter,
        //                p_closest_element->pGetGeometry(),
        //                p_closest_element->pGetProperties());

        //            // Deep copy elemental data and flags
        //            p_condition->Data() = p_closest_element->Data();
        //            p_condition->Set(Flags(*p_closest_element));

        //            new_condition_list.push_back(p_condition);

        //            rIdCounter++;
        //            continue;
        //        }
        //        KRATOS_ERROR << "No projected point found within tolerance!" << std::endl;
        //    }
        //    else
        //    {
        //        KRATOS_ERROR << "No Iga element in search radius found!" << std::endl;
        //    }
        //    rInterfaceConditionsModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());
        //}
    }

    void NurbsBrepModeler::GetInterfaceConditionsDEM(
        ModelPart& rExternalModelPart,
        ModelPart& rIgaModelPart,
        ModelPart& rInterfaceConditionsModelPart,
        const std::string& rConditionName,
        const int ShapeFunctionDerivativesOrder,
        const double SearchRadius,
        const double Accuracy,
        const double Tolerance,
        const double NumberOfIterations)
    {
        const Condition& ReferenceCondition = KratosComponents<Condition>::Get(rConditionName);

        BinsIgaConfigure::ContainerType BinsIgaObjects;
        BinsIgaObjects.resize(rIgaModelPart.NumberOfElements());
        const auto elements_begin = rIgaModelPart.Elements().ptr_begin();

        for (int i = 0; i < rIgaModelPart.NumberOfElements(); ++i)
        {
            auto it_elem = elements_begin + i;
            BinsIgaObjects[i] = (Kratos::make_shared<BinsIgaObject>(*it_elem));
        }
        auto iga_bins_structure = BinsObjectDynamic<BinsIgaConfigure>(BinsIgaObjects.begin(), BinsIgaObjects.end());

        int num_interface_obj_bin = BinsIgaObjects.size();

        BinsIgaConfigure::ResultContainerType neighbor_results(num_interface_obj_bin);
        std::vector<double> neighbor_distances(num_interface_obj_bin);

        //Instance of external interface nodes. Updated for each node.
        auto particle_object = Kratos::make_shared<BinsIgaObject>(array_1d<double, 3>(0.0));


        std::vector<int> checked_faces_trimming_curves;
        for (auto particle_element_ptr = rExternalModelPart.ElementsBegin();
            particle_element_ptr != rExternalModelPart.ElementsEnd();
            ++particle_element_ptr)
        {
            ModelPart::ConditionsContainerType new_condition_list;

            std::vector<Condition*> new_conditions;
            std::vector<array_1d<double, 3>> new_elastic_forces;
            std::vector<array_1d<double, 3>> new_total_forces;

            double radius = particle_element_ptr->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);

            auto results_itr = neighbor_results.begin();
            auto distance_itr = neighbor_distances.begin();

            //update coordinates of extern
            array_1d<double, 3> coords = particle_element_ptr->GetGeometry()[0].Coordinates();
            particle_object->UpdateCoordinates(coords);

            const std::size_t number_of_results = iga_bins_structure.SearchObjectsInRadius(
                particle_object,
                SearchRadius,
                results_itr,
                distance_itr,
                num_interface_obj_bin);

            // todo: neglect bad solutions
            if (number_of_results > 0)
            {
                Element::Pointer p_closest_element;
                double distance_to_node = 1e10;
                std::vector<Geometry<Node<3>>::Pointer> new_contact_geometries;
                for (int i = 0; i < number_of_results; ++i)
                {
                    int brep_id = neighbor_results[i]->pGetBaseElement()->GetValue(BREP_ID);

                    auto brep_face = GetBrepFace(brep_id);
                    auto surface = brep_face.GetSurface();

                    array_1d<double, 2> new_location(0.0);

                    if (IgaSurfaceUtilities::ProjectNodeOnSurface(
                        surface,
                        neighbor_results[i]->pGetBaseElement()->GetValue(LOCAL_COORDINATES),
                        particle_element_ptr->GetGeometry()[0].Coordinates(),
                        Accuracy,
                        Tolerance,
                        NumberOfIterations,
                        new_location))
                    {
                        if (!brep_face.IsInside(new_location))
                        {
                            for (int f = 0; f < checked_faces_trimming_curves.size(); ++f)
                            {
                                if (checked_faces_trimming_curves[f] == brep_face.Id())
                                {
                                    goto cnt;
                                }
                            }
                            checked_faces_trimming_curves.push_back(brep_face.Id());
                            auto projections = brep_face.ProjectOnFaceCurves(
                                particle_element_ptr->GetGeometry()[0].Coordinates(),
                                radius,
                                NumberOfIterations,
                                Accuracy);

                            KRATOS_WATCH(projections.size())
                            for (int j = 0; j < projections.size(); ++j)
                            {
                                array_1d<double, 3> point_3d = surface->PointAt(projections[j][0], projections[j][1]);

                                bool new_condition = true;
                                // check if new location is already existant
                                for (int locations_itr = 0; locations_itr < new_contact_geometries.size(); ++locations_itr)
                                {
                                    double distance_to_existant_point = norm_2(point_3d - new_contact_geometries[locations_itr]->Center());
                                    if (distance_to_existant_point < Tolerance)
                                    {
                                        new_condition = false;
                                    }
                                }

                                if (new_condition)
                                {
                                    auto new_geometry = IgaIntegrationPointUtilities::GetIntegrationPointSurface(
                                        surface,
                                        projections[j],
                                        ShapeFunctionDerivativesOrder
                                    );
                                    new_contact_geometries.push_back(new_geometry);
                                }
                            }
                            for (int j = 0; j < new_contact_geometries.size(); ++j)
                            {
                                KRATOS_WATCH(new_contact_geometries[j]->Center())
                            }
                        }
                        else
                        {
                            array_1d<double, 3> point = surface->PointAt(new_location[0], new_location[1]);

                            // check if new location is already existant
                            for (int locations_itr = 0; locations_itr < new_contact_geometries.size(); ++locations_itr)
                            {
                                double distance_to_existant_point = norm_2(point - new_contact_geometries[locations_itr]->Center());
                                if (distance_to_existant_point < Tolerance)
                                {
                                    goto cnt;
                                }
                            }

                            double distance_particle_to_contact = norm_2(point - particle_element_ptr->GetGeometry()[0].Coordinates());
                            if (distance_particle_to_contact <= radius)
                            {
                                auto new_geometry = IgaIntegrationPointUtilities::GetIntegrationPointSurface(
                                    surface,
                                    new_location,
                                    ShapeFunctionDerivativesOrder
                                );
                                new_contact_geometries.push_back(new_geometry);
                            }
                        }
                    }
                    cnt:;
                }
                Vector delete_condition = ZeroVector(particle_element_ptr->GetValue(WALL_POINT_CONDITION_POINTERS).size());

                for (auto p_contact_geometry : new_contact_geometries)
                {
                    int i = 0;
                    bool found_old_condition = false;
                    for ( auto interface_condition = rInterfaceConditionsModelPart.ConditionsBegin();
                        interface_condition != rInterfaceConditionsModelPart.ConditionsEnd();
                        ++interface_condition)
                    {
                        array_1d<double, 3> distance_vector = p_contact_geometry->Center() - (*&interface_condition)->GetGeometry().Center();

                        if (norm_2(distance_vector) < 220*Tolerance)
                        {
                            interface_condition->pGetGeometry() = p_contact_geometry;

                            delete_condition[i] = 1;
                            new_elastic_forces.push_back(particle_element_ptr->GetValue(WALL_POINT_CONDITION_ELASTIC_FORCES)[i]);
                            new_total_forces.push_back(particle_element_ptr->GetValue(WALL_POINT_CONDITION_TOTAL_FORCES)[i]);

                            new_conditions.push_back(&*interface_condition);

                            found_old_condition = true;
                            break;
                        }
                        else
                        {
                            KRATOS_WATCH(norm_2(distance_vector))
                        }
                        i++;
                    }
                    if (!found_old_condition)
                    {
                        int id = 1;
                        if (rInterfaceConditionsModelPart.GetRootModelPart().Conditions().size() > 0)
                            id = rInterfaceConditionsModelPart.GetRootModelPart().Conditions().back().Id() + 1;

                        auto p_condition = ReferenceCondition.Create(
                            id,
                            p_contact_geometry,
                            particle_element_ptr->pGetProperties());

                        new_elastic_forces.push_back(ZeroVector(3));
                        new_total_forces.push_back(ZeroVector(3));

                        new_condition_list.push_back(p_condition);
                        new_conditions.push_back(&*p_condition);
                    }
                }

                int condition_length = particle_element_ptr->GetValue(WALL_POINT_CONDITION_POINTERS).size();

                for (int i = 0; i < condition_length; ++i)
                {
                    if (!delete_condition[i] > 0.1)
                    {
                        rInterfaceConditionsModelPart.RemoveConditionFromAllLevels(particle_element_ptr->GetValue(WALL_POINT_CONDITION_POINTERS)[i]->Id());
                    }
                }

                rInterfaceConditionsModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());

                particle_element_ptr->SetValue(WALL_POINT_CONDITION_ELASTIC_FORCES, new_elastic_forces);
                particle_element_ptr->SetValue(WALL_POINT_CONDITION_TOTAL_FORCES, new_total_forces);
                particle_element_ptr->SetValue(WALL_POINT_CONDITION_POINTERS, new_conditions);
                KRATOS_WATCH(particle_element_ptr->GetValue(WALL_POINT_CONDITION_POINTERS).size())
            }
        }
    }

    void NurbsBrepModeler::ExportGeometry()
    {
        std::cout << "ExportGeometry to FILE..." << std::endl;
        BrepJsonIO a; 
        //a.ExportNurbsGeometry(m_brep_model_vector);         
    }

    NurbsBrepModeler::NurbsBrepModeler(ModelPart& rModelPart)
        : m_model_part(rModelPart)
    {
    }
}  // namespace Kratos.


