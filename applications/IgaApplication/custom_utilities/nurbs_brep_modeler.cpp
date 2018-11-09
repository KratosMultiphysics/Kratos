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
        std::vector<BrepModel> brep_model_vector = rBrepJsonIO.ImportNurbsBrepGeometry(m_model_part, rNurbsBrepGeometryJson);
        for (auto brep_model = brep_model_vector.begin(); brep_model != brep_model_vector.end(); ++brep_model)
        {
            m_brep_model_vector.push_back(*brep_model);
        }
    }

    void NurbsBrepModeler::ImportModelPart(ModelPart& model_part, Parameters& rModelPartParameters)
    {
        for (int i = 0; i < rModelPartParameters["element_condition_list"].size(); ++i)
        {
            Parameters element_parameter = rModelPartParameters["element_condition_list"][i];
            std::string type = element_parameter["type"].GetString();
            std::string name = element_parameter["name"].GetString();
            int property_id = element_parameter["parameters"]["properties_id"].GetInt();
            int shape_function_derivatives_order = element_parameter["parameters"]["shape_function_derivatives_order"].GetInt();

            std::string iga_model_part = element_parameter["iga_model_part"].GetString();
            std::string geometry_type = element_parameter["geometry_type"].GetString();

            std::vector<std::string> variable_list;
            for (int j = 0; j < element_parameter["parameters"]["variables"].size(); ++j)
            {
                variable_list.push_back(element_parameter["parameters"]["variables"][j].GetString());
            }

            std::string sub_model_part_name = element_parameter["iga_model_part"].GetString();

            if (!model_part.HasSubModelPart(sub_model_part_name))
            {
                model_part.CreateSubModelPart(sub_model_part_name);
            }
            ModelPart& sub_model_part = model_part.GetSubModelPart(sub_model_part_name);
            for (int j = 0; j < element_parameter["parameters"]["brep_ids"].size(); ++j)
            {
                int brep_id = element_parameter["parameters"]["brep_ids"][j].GetInt();
                bool success = m_brep_model_vector[0].GetIntegrationDomain(
                    sub_model_part, brep_id, type, name, property_id, shape_function_derivatives_order, variable_list);
            }
        }
    }

    NurbsBrepModeler::NurbsBrepModeler(ModelPart& rModelPart)
        : m_model_part(rModelPart)
    {
    }
}  // namespace Kratos.


