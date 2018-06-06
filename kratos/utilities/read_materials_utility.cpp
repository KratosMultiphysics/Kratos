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
//

// System includes

// External includes

// Project includes
#include "includes/properties.h"
#include "utilities/read_materials_utility.hpp"

namespace Kratos
{

ReadMaterialsUtility::ReadMaterialsUtility(
    Parameters& rParameters,
    Model& rModel
    ) : mrModel(rModel)
{
    Parameters default_parameters(R"(
    {
        "materials_filename" : "please specify the file to be opened"
    }  )"
    );

    rParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // read json string in materials file, create Parameters
    const std::string& materials_filename = rParameters["materials_filename"].GetString();
    std::ifstream infile(materials_filename);
    std::stringstream buffer;
    buffer << infile.rdbuf();
    Parameters materials(buffer.str());

    GetPropertyBlock(materials);
}

/***********************************************************************************/
/***********************************************************************************/

ReadMaterialsUtility::ReadMaterialsUtility(
    const std::string& rParametersStr,
    Model& rModel
    ) : mrModel(rModel)
{
    // Receive json string with materials properties, create Parameters
    Parameters materials(rParametersStr);

    GetPropertyBlock(materials);
}

/***********************************************************************************/
/***********************************************************************************/

void ReadMaterialsUtility::GetPropertyBlock(Parameters& materials)
{
    KRATOS_INFO("Read materials") << "Started" << std::endl;
    for (auto i = 0; i < materials["properties"].size(); ++i) {
        Parameters material = materials["properties"].GetArrayItem(i);
        AssignPropertyBlock(material);
    }
    KRATOS_INFO("Read materials") << "Finished" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

std::string CleanVariableName(std::string);
std::string CleanVariableName(std::string line){
    std::stringstream ss(line);
    std::string variable_name;
    while (std::getline(ss, variable_name, '.')){}
    return variable_name;
}

/***********************************************************************************/
/***********************************************************************************/

void ReadMaterialsUtility::AssignPropertyBlock(Parameters& data)
{
    // Get the properties for the specified model part.
    ModelPart& model_part = mrModel.GetModelPart(data["model_part_name"].GetString());
    const IndexType property_id = data["properties_id"].GetInt();
    const IndexType mesh_id = 0;
    Properties::Pointer p_prop = model_part.pGetProperties(property_id, mesh_id);

    /*
    //TODO(marcelo): Implement the "keys()" part? Not sure if this check is necessary
    //if (data["Material"]["Variables"].end() - data["Material"]["Variables"].begin())
    //    KRATOS_INFO("::[Reading materials process DEBUG]::")
    //            << "Property " << property_id << " is not empty." << std::endl;
    if (p_prop->HasVariables())
        KRATOS_INFO("Read materials")
            << "Property " << property_id << " already has variables." << std::endl;
    //if (len(data["Material"]["Tables"].keys()) > 0 && p_prop.HasTables())
    if (p_prop->HasTables())
        KRATOS_INFO("Read materials")
            << "Property " << property_id << " already has tables." << std::endl;
    */

    // Assign the p_properties to the model part's elements and conditions.
    auto& elements_array = model_part.Elements();
    auto& conditions_array = model_part.Conditions();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(elements_array.size()); ++i) {
        auto it_elem = elements_array.begin() + i;
        it_elem->SetProperties(p_prop);
    }

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i) {
        auto it_cond = conditions_array.begin() + i;
        it_cond->SetProperties(p_prop);
    }

    //Set the CONSTITUTIVE_LAW for the current p_properties.
    if (data["Material"].Has("constitutive_law")) {
        std::string constitutive_law_name = data["Material"]["constitutive_law"]["name"].GetString();

        // Remove application info from consitutive law name.
        // Ex: KratosMultiphysics.StructuralMechanicsApplication.LinearElastic3D -> LinearElastic3D
        constitutive_law_name = CleanVariableName(constitutive_law_name);

        auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get(constitutive_law_name).Clone();
        p_prop->SetValue(CONSTITUTIVE_LAW, p_constitutive_law);
    } else {
        KRATOS_INFO("Read materials") << "Not consitutive law defined for material ID: " << property_id << std::endl;
    }

    // Add / override the values of material parameters in the p_properties
    auto variables = data["Material"]["Variables"];
    for(auto iter = variables.begin(); iter != variables.end(); iter++) {
        const auto& value = variables.GetValue(iter.name());

        // Remove application info from variable name.
        // Ex: KratosMultiphysics.YOUNG_MODULUS -> YOUNG_MODULUS
        const std::string& variable_name = CleanVariableName(iter.name());

        if (value.IsDouble()){
            const auto variable = KratosComponents<Variable<double>>().Get(variable_name);
            p_prop->SetValue(variable, value.GetDouble());
        } else if (value.IsInt()){
            const auto variable = KratosComponents<Variable<int>>().Get(variable_name);
            p_prop->SetValue(variable, value.GetInt());
        } else if (value.IsBool()){
            const auto variable = KratosComponents<Variable<bool>>().Get(variable_name);
            p_prop->SetValue(variable, value.GetBool());
        } else if (value.IsString()){
            const auto variable = KratosComponents<Variable<std::string>>().Get(variable_name);
            p_prop->SetValue(variable, value.GetString());
        } else if (value.IsVector()){
            const auto variable = KratosComponents<Variable<Vector>>().Get(variable_name);
            p_prop->SetValue(variable, value.GetVector());
        } else if (value.IsMatrix()){
            const auto variable = KratosComponents<Variable<Matrix>>().Get(variable_name);
            p_prop->SetValue(variable, value.GetMatrix());
        }
        else {
            KRATOS_ERROR << "Type of value is not available";
        }
    }

    // Add / override tables in the p_properties
    auto tables = data["Material"]["Tables"];
    for(auto iter = tables.begin(); iter != tables.end(); iter++) {
        auto table_param = tables.GetValue(iter.name());
        // Case table is double, double. How is it defined? How to check?
        Table<double> table;

        // Remove application info from variable name.
        // Ex: KratosMultiphysics.YOUNG_MODULUS -> YOUNG_MODULUS
        const std::string& input_var_name = CleanVariableName(table_param["input_variable"].GetString());
        const std::string& output_var_name = CleanVariableName(table_param["output_variable"].GetString());

        const auto input_var = KratosComponents<Variable<double>>().Get(input_var_name);
        const auto output_var = KratosComponents<Variable<double>>().Get(output_var_name);
        for (auto i = 0; i < table_param["data"].size(); i++) {
            table.insert(table_param["data"][i][0].GetDouble(),
                            table_param["data"][i][1].GetDouble());
        }
        p_prop->SetTable(input_var, output_var, table);
    }
}

}  // namespace Kratos.
