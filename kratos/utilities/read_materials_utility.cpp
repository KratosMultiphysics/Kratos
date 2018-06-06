//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Marcelo Raschi


// System includes

// External includes

// Project includes
#include "includes/properties.h"
#include "utilities/read_materials_utility.hpp"

namespace Kratos
{
    ReadMaterialsUtility::ReadMaterialsUtility(Parameters parameters, Model &rModel)
            : mrModel(rModel)
    {
        Parameters default_parameters(R"(
            {
                "materials_filename" : "please specify the file to be opened"
            }  )"
        );

        parameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        // read json string in materials file, create Parameters
        std::string materials_filename = parameters["materials_filename"].GetString();
        std::ifstream infile(materials_filename);
        std::stringstream buffer;
        buffer << infile.rdbuf();
        Parameters materials(buffer.str());

        GetPropertyBlock(materials);
    }

    ReadMaterialsUtility::ReadMaterialsUtility(std::string parameters_str, Model &rModel)
            : mrModel(rModel)
    {
        // receive json string with materials properties, create Parameters
        Parameters materials(parameters_str);

        GetPropertyBlock(materials);
    }

    void ReadMaterialsUtility::GetPropertyBlock(Parameters materials)
    {
        KRATOS_INFO("Read materials") << "Started" << std::endl;
        for (auto i = 0; i < materials["properties"].size(); ++i)
        {
            Parameters material = materials["properties"].GetArrayItem(i);
            AssignPropertyBlock(material);
        }
        KRATOS_INFO("Read materials") << "Finished" << std::endl;
    }

    void ReadMaterialsUtility::AssignPropertyBlock(Parameters data)
    {
        // Get the properties for the specified model part.
        ModelPart& model_part = mrModel.GetModelPart(data["model_part_name"].GetString());
        IndexType property_id = data["properties_id"].GetInt();
        IndexType mesh_id = 0;
        Properties::Pointer prop = model_part.pGetProperties(property_id, mesh_id);

        //TODO(marcelo): Implement the "keys()" part? Not sure if this check is necessary
        //if (data["Material"]["Variables"].end() - data["Material"]["Variables"].begin())
        //    KRATOS_INFO("::[Reading materials process DEBUG]::")
        //            << "Property " << property_id << " is not empty." << std::endl;
        if (prop->HasVariables())
            KRATOS_INFO("Read materials")
                << "Property " << property_id << " already has variables." << std::endl;
        //if (len(data["Material"]["Tables"].keys()) > 0 && prop.HasTables())
        if (prop->HasTables())
            KRATOS_INFO("Read materials")
                << "Property " << property_id << " already has tables." << std::endl;

        // Assign the properties to the model part's elements and conditions.
        for (auto elem = model_part.ElementsBegin(); elem != model_part.ElementsEnd(); elem++)
            elem->SetProperties(prop);

        for (auto cond = model_part.ConditionsBegin(); cond != model_part.ConditionsEnd(); cond++)
            cond->SetProperties(prop);

        //Set the CONSTITUTIVE_LAW for the current properties.
        std::string constitutive_law_name = data["Material"]["constitutive_law"]["name"].GetString();
        auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get(constitutive_law_name).Clone();
        prop->SetValue(CONSTITUTIVE_LAW, p_constitutive_law);

        // Add / override the values of material parameters in the properties
        auto variables = data["Material"]["Variables"];
        for(auto iter = variables.begin(); iter != variables.end(); iter++)
        {
            auto value = variables.GetValue(iter.name());
            if (value.IsDouble()){
                auto variable = KratosComponents<Variable<double>>().Get(iter.name());
                prop->SetValue(variable, value.GetDouble());
            }
            else if (value.IsInt()){
                auto variable = KratosComponents<Variable<int>>().Get(iter.name());
                prop->SetValue(variable, value.GetInt());
            }
            else if (value.IsBool()){
                auto variable = KratosComponents<Variable<bool>>().Get(iter.name());
                prop->SetValue(variable, value.GetBool());
            }
            else if (value.IsString()){
                auto variable = KratosComponents<Variable<std::string>>().Get(iter.name());
                prop->SetValue(variable, value.GetString());
            }
            else if (value.IsVector()){
                auto variable = KratosComponents<Variable<Vector>>().Get(iter.name());
                prop->SetValue(variable, value.GetVector());
            }
            else if (value.IsMatrix()){
                auto variable = KratosComponents<Variable<Matrix>>().Get(iter.name());
                prop->SetValue(variable, value.GetMatrix());
            }
                //TODO(marcelo): Add error here
            else {
                KRATOS_INFO("Read materials")
                        << "Type of value is not available" << std::endl;
            }
        }

        // Add / override tables in the properties
        auto tables = data["Material"]["Tables"];
        for(auto iter = tables.begin(); iter != tables.end(); iter++)
        {
            auto table_param = tables.GetValue(iter.name());
            // Case table is double, double. How is it defined? How to check?
                Table<double> table;
                auto input_var = KratosComponents<Variable<double>>().Get(
                    table_param["input_variable"].GetString());
                auto output_var = KratosComponents<Variable<double>>().Get(
                    table_param["output_variable"].GetString());
                for (auto i = 0; i < table_param["data"].size(); i++) {
                    table.insert(table_param["data"][i][0].GetDouble(),
                                 table_param["data"][i][1].GetDouble());
                }
                prop->SetTable(input_var, output_var, table);
    //            KRATOS_WATCH(table);
        }
    }

}  // namespace Kratos.
