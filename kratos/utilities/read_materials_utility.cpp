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
    Parameters rParameters,
    Model& rModel
    ) : mrModel(rModel)
{
    Parameters default_parameters(R"(
    {
        "Parameters" : {
            "materials_filename" : "please specify the file to be opened"
        }
    }  )"
    );

    rParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // read json string in materials file, create Parameters
    const std::string& materials_filename = rParameters["Parameters"]["materials_filename"].GetString();
    std::ifstream infile(materials_filename);
    KRATOS_ERROR_IF_NOT(infile.good()) << "The masterials file cannot be found: " << materials_filename << std::endl;
    std::stringstream buffer;
    buffer << infile.rdbuf();
    Parameters materials(buffer.str());

    GetPropertyBlock(materials);
}

/***********************************************************************************/
/***********************************************************************************/

ReadMaterialsUtility::ReadMaterialsUtility(
    const std::string& rParametersName,
    Model& rModel
    ) : mrModel(rModel)
{
    // Receive json string with materials properties, create Parameters
    Parameters materials(rParametersName);

    GetPropertyBlock(materials);
}

/***********************************************************************************/
/***********************************************************************************/

void ReadMaterialsUtility::GetPropertyBlock(Parameters materials)
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

void ReadMaterialsUtility::AssignPropertyBlock(Parameters data)
{
    // Get the properties for the specified model part.
    ModelPart& model_part = mrModel.GetModelPart(data["model_part_name"].GetString());
    const IndexType property_id = data["properties_id"].GetInt();
    const IndexType mesh_id = 0;
    Properties::Pointer p_prop = model_part.pGetProperties(property_id, mesh_id);

    KRATOS_INFO_IF("::[Reading materials process]:: Property", data["Material"]["Variables"].size() > 0 && p_prop->HasVariables()) << std::to_string(property_id) << "already has variables." << std::endl;
    KRATOS_INFO_IF("::[Reading materials process]:: Property", data["Material"]["Tables"].size() > 0 && p_prop->HasTables()) << std::to_string(property_id) << "already has tables." << std::endl;

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
    Parameters variables = data["Material"]["Variables"];
    for(auto iter = variables.begin(); iter != variables.end(); iter++) {
        const Parameters value = variables.GetValue(iter.name());

        // Remove application info from variable name.
        // Ex: KratosMultiphysics.YOUNG_MODULUS -> YOUNG_MODULUS
        const std::string& variable_name = CleanVariableName(iter.name());

        // We don't just copy the values, we do some tyransformation depending of the destination variable
        if(KratosComponents<Variable<double> >::Has(variable_name)) {
            const Variable<double>& variable = KratosComponents<Variable<double>>().Get(variable_name);
            if (value.IsDouble()) {
                p_prop->SetValue(variable, value.GetDouble());
            } else if (value.IsInt()) {
                p_prop->SetValue(variable, static_cast<double>(value.GetInt()));
            } else {
                KRATOS_ERROR << "Check you write the value in a correct format: " << value << std::endl;
            }
        } else if(KratosComponents<Variable<bool> >::Has(variable_name)) {
            const Variable<bool>& variable = KratosComponents<Variable<bool>>().Get(variable_name);
            if (value.IsBool()) {
                p_prop->SetValue(variable, value.GetBool());
            } else if (value.IsInt()) {
                p_prop->SetValue(variable, static_cast<bool>(value.GetInt()));
            } else if (value.IsDouble()) {
                p_prop->SetValue(variable, static_cast<bool>(value.GetDouble()));
            } else {
                KRATOS_ERROR << "Check you write the value in a correct format: " << value << std::endl;
            }
        } else if(KratosComponents<Variable<int> >::Has(variable_name)) {
            const Variable<int>& variable = KratosComponents<Variable<int>>().Get(variable_name);
            if (value.IsInt()) {
                p_prop->SetValue(variable, value.GetInt());
            } else if (value.IsDouble()) {
                p_prop->SetValue(variable, static_cast<int>(value.GetDouble()));
            } else {
                KRATOS_ERROR << "Check you write the value in a correct format: " << value << std::endl;
            }
        } else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name)) {
            const Variable<array_1d<double, 3>>& variable = KratosComponents<Variable<array_1d<double, 3>>>().Get(variable_name);
            if (value.IsVector()) {
                array_1d<double, 3> temp(3, 0.0);
                const Vector& value_variable = value.GetVector();
                const std::size_t iter_number = (3 < value_variable.size()) ? 3 : value_variable.size();
                for (std::size_t index = 0; index < iter_number; index++)
                    temp[index] = value_variable[index];
                p_prop->SetValue(variable, temp);
            } else {
                KRATOS_ERROR << "Check you write the value in a correct format: " << value << std::endl;
            }
        } else if(KratosComponents<Variable<array_1d<double, 6> > >::Has(variable_name)) {
            const Variable<array_1d<double, 6>>& variable = KratosComponents<Variable<array_1d<double, 6>>>().Get(variable_name);
            if (value.IsVector()) {
                array_1d<double, 6> temp(6, 0.0);
                const Vector& value_variable = value.GetVector();
                const std::size_t iter_number = (6 < value_variable.size()) ? 6 : value_variable.size();
                for (std::size_t index = 0; index < iter_number; index++)
                    temp[index] = value_variable[index];
                p_prop->SetValue(variable, temp);
            } else {
                KRATOS_ERROR << "Check you write the value in a correct format: " << value << std::endl;
            }
        } else if(KratosComponents<Variable<Vector > >::Has(variable_name)) {
            const Variable<Vector>& variable = KratosComponents<Variable<Vector>>().Get(variable_name);
            if (value.IsVector()) {
                p_prop->SetValue(variable, value.GetVector());
            } else if (value.IsMatrix()) {
                Vector temp;
                const Matrix& value_variable = value.GetMatrix();
                for (std::size_t index = 0; index < value_variable.size1(); index++)
                    temp[index] = value_variable(index, 0);
                p_prop->SetValue(variable, temp);
            } else {
                KRATOS_ERROR << "Check you write the value in a correct format: " << value << std::endl;
            }
        } else if(KratosComponents<Variable<Matrix> >::Has(variable_name)) {
            const Variable<Matrix>& variable = KratosComponents<Variable<Matrix>>().Get(variable_name);
            if (value.IsMatrix()) {
                p_prop->SetValue(variable, value.GetMatrix());
            } else if (value.IsVector()) {
                Matrix temp;
                const Vector& value_variable = value.GetVector();
                for (std::size_t index = 0; index < value_variable.size(); index++)
                    temp(index, 0) = value_variable[index];
                p_prop->SetValue(variable, temp);
            } else {
                KRATOS_ERROR << "Check you write the value in a correct format: " << value << std::endl;
            }
        } else if(KratosComponents<Variable<std::string> >::Has(variable_name)) {
            const Variable<std::string>& variable = KratosComponents<Variable<std::string>>().Get(variable_name);
            if (value.IsString()) {
                p_prop->SetValue(variable, value.GetString());
            } else {
                KRATOS_ERROR << "Check you write the value in a correct format: " << value << std::endl;
            }
        } else {
            KRATOS_ERROR << "Type of value is not available";
        }
    }

    // Add / override tables in the p_properties
    Parameters tables = data["Material"]["Tables"];
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
