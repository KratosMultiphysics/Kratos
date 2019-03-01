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

    //void NurbsBrepModeler::ImportModelPart(ModelPart& model_part, Parameters& rModelPartParameters)
    //{
    //    for (int i = 0; i < rModelPartParameters["element_condition_list"].size(); ++i)
    //    {
    //        Parameters element_parameter = rModelPartParameters["element_condition_list"][i];

    //        std::string sub_model_part_name = element_parameter["iga_model_part"].GetString();
    //        ModelPart& sub_model_part = model_part.HasSubModelPart(sub_model_part_name)
    //            ? model_part.GetSubModelPart(sub_model_part_name)
    //            : model_part.CreateSubModelPart(sub_model_part_name);

    //        std::string geometry_type = element_parameter["geometry_type"].GetString();

    //        if (geometry_type == "Geometry3DStrong")
    //        {
    //            bool success = false;
    //            for (int j = 0; j < element_parameter["brep_ids"].size(); ++j)
    //            {
    //                int brep_id = element_parameter["brep_ids"][j].GetInt();
    //                Vector parameter = element_parameter["parameters"]["local_parameters"].GetVector();
    //                success = m_brep_model_vector[0].GetNodesGeometry(
    //                    sub_model_part, brep_id, parameter[0], parameter[1]);
    //            }
    //            continue;
    //        }
    //        else
    //        {
    //            std::string type = element_parameter["parameters"]["type"].GetString();
    //            std::string name = element_parameter["parameters"]["name"].GetString();
    //            int property_id = element_parameter["parameters"]["properties_id"].GetInt();
    //            int shape_function_derivatives_order = element_parameter["parameters"]["shape_function_derivatives_order"].GetInt();


    //            std::vector<std::string> variable_list;
    //            for (int j = 0; j < element_parameter["parameters"]["variables"].size(); ++j)
    //            {
    //                variable_list.push_back(element_parameter["parameters"]["variables"][j].GetString());
    //            }

    //            if (element_parameter.Has("brep_ids"))
    //            {
    //                for (int j = 0; j < element_parameter["brep_ids"].size(); ++j)
    //                {
    //                    int brep_id = element_parameter["brep_ids"][j].GetInt();
    //                    if (geometry_type == "Geometry3D")
    //                        bool success = m_brep_model_vector[0].GetIntegrationDomainGeometry(
    //                            sub_model_part, brep_id, type, name,
    //                            property_id, shape_function_derivatives_order, variable_list);
    //                    if (geometry_type == "BrepCoupling")
    //                        bool success = m_brep_model_vector[0].GetIntegrationDomainBrepCoupling(
    //                            sub_model_part, brep_id, type, name,
    //                            property_id, shape_function_derivatives_order, variable_list);
    //                    if (geometry_type == "Brep")
    //                        bool success = m_brep_model_vector[0].GetIntegrationDomainBrep(
    //                            sub_model_part, brep_id, type, name,
    //                            property_id, shape_function_derivatives_order, variable_list);
    //                }
    //            }
    //            else
    //            {
    //                if (geometry_type == "BrepCoupling")
    //                    bool success = m_brep_model_vector[0].GetIntegrationDomainBrepCoupling(
    //                        sub_model_part, type, name,
    //                        property_id, shape_function_derivatives_order, variable_list);
    //            }
    //        }
    //    }
    //}


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

                if (element_parameter.Has("brep_ids"))
                {
                    for (int j = 0; j < element_parameter["brep_ids"].size(); ++j)
                    {
                        int brep_id = element_parameter["brep_ids"][j].GetInt();
                        if (geometry_type == "GeometryFace")
                        {
                            GetBrepFace(brep_id).GetGeometryIntegrationTrimmed(
                                sub_model_part, type, name,
                                shape_function_derivatives_order,
                                variable_list);
                        }
                        if (geometry_type == "SurfaceEdge")
                        {
                            GetBrepEdge(brep_id);
                            bool success = m_brep_model_vector[0].GetIntegrationDomainBrep(
                                sub_model_part, brep_id, type, name,
                                shape_function_derivatives_order, variable_list);
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

                                    auto element_vector = IgaIntegrationUtilities::GetIntegrationDomainSurfaceEdgeSurfaceEdge(
                                        surface_1,
                                        trimming_curve_1,
                                        surface_2,
                                        trimming_curve_2,
                                        shape_function_derivatives_order);

                                    if (type == "element")
                                    {
                                        int id = 0;
                                        if (sub_model_part.GetRootModelPart().Elements().size() > 0)
                                            int id = sub_model_part.GetRootModelPart().Elements().back().Id() + 1;

                                        IgaIntegrationUtilities::ChangeElementType(
                                            element_vector, sub_model_part, name, id);
                                    }

                                    if (type == "condition")
                                    {
                                        int id = 0;
                                        if (sub_model_part.GetRootModelPart().Conditions().size() > 0)
                                            int id = sub_model_part.GetRootModelPart().Conditions().back().Id() + 1;

                                        IgaIntegrationUtilities::ChangeConditionType(
                                            element_vector, sub_model_part, name, id);
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
            }
        }
    }

    void NurbsBrepModeler::ExportGeometry()
    {
        std::cout << "ExportGeometry to FILE..." << std::endl;
        BrepJsonIO a; 
        //a.ExportNurbsGeometry(m_brep_model_vector);         
    }


    // void NurbsBrepModeler::PrintBrepNodes()
    // {   
    //     for (int i = 0; i < m_brep_model_vector.size(); ++i)
    //     {
    //         m_brep_model_vector[i].GetModelNodes();
    //     }
    // }

    // void NurbsBrepModeler::PrintEdgePolygon()
    // {   
    //     for (int i = 0; i < m_brep_model_vector.size(); ++i)
    //     {
    //         m_brep_model_vector[i].PrintEdgeNodes();
    //     }
    // }

    // void NurbsBrepModeler::PrintTrimmingPolygon()
    // {
    //     for (int i = 0; i < m_brep_model_vector.size(); ++i)
    //     {
    //         m_brep_model_vector[i].PrintTrimmingNodes();
    //     }
    // }

    NurbsBrepModeler::NurbsBrepModeler(ModelPart& rModelPart)
        : m_model_part(rModelPart)
    {
    }
}  // namespace Kratos.


