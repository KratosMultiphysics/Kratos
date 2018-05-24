//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//

// System includes


// External includes


// Project includes
#include "custom_utilities/read_materials.hpp"

namespace Kratos
{
    ReadMaterialProcess::ReadMaterialProcess(ModelPart &rModelPart,
                                             Parameters parameters)
            : mrModelPart(rModelPart)
    {
        Parameters default_parameters(R"(
            {
                "materials_filename" : "please specify the file to be opened"
            }  )"
        );

        parameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        KRATOS_INFO("::[Reading materials process DEBUG]::") << "Started" << std::endl;

        // read json string in materials file, create Parameters
        std::string materials_filename = parameters["materials_filename"].GetString();
        std::ifstream infile(materials_filename);
        std::stringstream buffer;
        buffer << infile.rdbuf();
        Parameters materials(buffer.str());

        for (auto i = 0; i < materials["properties"].size(); ++i)
        {
            Parameters material = materials["properties"].GetArrayItem(i);
            AssignPropertyBlock(material);
        }

        KRATOS_INFO("::[Reading materials process DEBUG]::") << "Finished" << std::endl;
    }

    void ReadMaterialProcess::AssignPropertyBlock(Parameters data)
    {
        // Get the properties for the specified model part.
        IndexType property_id = data["properties_id"].GetInt();
        IndexType mesh_id = 0;
        Properties::Pointer prop = mrModelPart.pGetProperties(property_id, mesh_id);

        //TODO(marcelo): Implement the "keys()" part
        //if (data["Material"]["Variables"].end() - data["Material"]["Variables"].begin())
        //    KRATOS_INFO("::[Reading materials process DEBUG]::")
        //            << "Property " << property_id << " is not empty." << std::endl;
        if (prop->HasVariables())
            KRATOS_INFO("::[Reading materials process DEBUG]::")
                << "Property " << property_id << " already has variables." << std::endl;
        //if (len(data["Material"]["Tables"].keys()) > 0 && prop.HasTables())
        if (prop->HasTables())
            KRATOS_INFO("::[Reading materials process DEBUG]::")
                << "Property " << property_id << " already has tables." << std::endl;

        // Assign the properties to the model part's elements and conditions.
        for (auto elem = mrModelPart.ElementsBegin(); elem != mrModelPart.ElementsEnd(); elem++)
            elem->SetProperties(prop);

        for (auto cond = mrModelPart.ConditionsBegin(); cond != mrModelPart.ConditionsEnd(); cond++)
            cond->SetProperties(prop);

        //Set the CONSTITUTIVE_LAW for the current properties.
        std::string constitutive_law_name = data["Material"]["constitutive_law"]["name"].GetString();
        auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get(constitutive_law_name).Clone();
        prop->SetValue(CONSTITUTIVE_LAW, p_constitutive_law);

        // Add / override the values of material parameters in the properties
        for(auto iter = data["Material"]["Variables"].begin();
            iter != data["Material"]["Variables"].end(); iter++){
            auto value = data["Material"]["Variables"].GetValue(iter.name());
            auto variable = KratosComponents<VariableData>().Get(iter.name());
            if (value.IsDouble()){
                //prop->SetValue(variable, value.GetDouble());
            }
            /*
            else if (value.IsInt())
                prop->SetValue(iter, value.GetInt());
            else if (value.IsBool())
                prop->SetValue(iter, value.GetBool());
            else if (value.IsString())
                prop->SetValue(iter, value.GetString());
            else if (value.IsMatrix())
                prop->SetValue(iter, value.GetMatrix());
            else if (value.IsVector())
                prop->SetValue(iter, value.GetVector());
            else{
                KRATOS_INFO("Type of value is not available");
            }
                */
        }
/*
        # Add / override tables in the properties
        for key, table in mat["Tables"].items():
            table_name = key

            input_var = self._GetVariable(table["input_variable"].GetString())
            output_var = self._GetVariable(table["output_variable"].GetString())

            new_table = KratosMultiphysics.PiecewiseLinearTable()

            for i in range(table["data"].size()):
                new_table.AddRow(table["data"][i][0].GetDouble(), table["data"][i][1].GetDouble())

            prop.SetTable(input_var,output_var,new_table)
        */
    }
}  // namespace Kratos.
