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
//                   Vicente Mataix Ferrandiz
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
    KRATOS_ERROR_IF_NOT(infile.good()) << "Materials file: " << materials_filename << " cannot be found" << std::endl;
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

// Remove application info from variable name.
// Ex: KratosMultiphysics.YOUNG_MODULUS -> YOUNG_MODULUS
// Ex: KratosMultiphysics.StructuralMechanicsApplication.LinearElastic3D -> LinearElastic3D

void TrimComponentName(std::string&);
void TrimComponentName(std::string& line){
    std::stringstream ss(line);
    std::size_t counter = 0;
    while (std::getline(ss, line, '.')){counter++;}
    if (counter > 1)
        KRATOS_WARNING("Read materials") << "Ignoring module information for component " << line << std::endl;
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

    // Compute the size using the iterators
    std::size_t variables_size = 0;
    for(auto it=data["Material"]["Variables"].begin(); it!=data["Material"]["Variables"].end(); ++it)
        variables_size++;
    
    std::size_t tables_size = 0;
    for(auto it=data["Material"]["Tables"].begin(); it!=data["Material"]["Tables"].end(); ++it)
        tables_size++;
    
    KRATOS_WARNING_IF("Read materials", variables_size > 0 && p_prop->HasVariables())
        << "Property " << std::to_string(property_id) << " already has variables." << std::endl;
    KRATOS_WARNING_IF("Read materials", tables_size > 0 && p_prop->HasTables())
        << "Property " << std::to_string(property_id) << " already has tables." << std::endl;

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
        TrimComponentName(constitutive_law_name);

        auto p_constitutive_law = KratosComponents<ConstitutiveLaw>().Get(constitutive_law_name).Clone();
        p_prop->SetValue(CONSTITUTIVE_LAW, p_constitutive_law);
    } else {
        KRATOS_INFO("Read materials") << "No constitutive law defined for material ID: " << property_id << std::endl;
    }

    // Add / override the values of material parameters in the p_properties
    Parameters variables = data["Material"]["Variables"];
    for(auto iter = variables.begin(); iter != variables.end(); iter++) {
        const Parameters value = variables.GetValue(iter.name());

        std::string variable_name = iter.name();
        TrimComponentName(variable_name);

        // We don't just copy the values, we do some tyransformation depending of the destination variable
        if(KratosComponents<Variable<double> >::Has(variable_name)) {
            const Variable<double>& variable = KratosComponents<Variable<double>>().Get(variable_name);
            if (value.IsDouble()) {
                p_prop->SetValue(variable, value.GetDouble());
            } else if (value.IsInt()) {
                p_prop->SetValue(variable, static_cast<double>(value.GetInt()));
            } else {
                KRATOS_ERROR << "Check the value: " << value << " is in the correct format" << std::endl;
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
                KRATOS_ERROR << "Check the value: " << value << " is in the correct format" << std::endl;
            }
        } else if(KratosComponents<Variable<int> >::Has(variable_name)) {
            const Variable<int>& variable = KratosComponents<Variable<int>>().Get(variable_name);
            if (value.IsInt()) {
                p_prop->SetValue(variable, value.GetInt());
            } else if (value.IsDouble()) {
                p_prop->SetValue(variable, static_cast<int>(value.GetDouble()));
            } else {
                KRATOS_ERROR << "Check the value: " << value << " is in the correct format" << std::endl;
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
                KRATOS_ERROR << "Check the value: " << value << " is in the correct format" << std::endl;
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
                KRATOS_ERROR << "Check the value: " << value << " is in the correct format" << std::endl;
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
                KRATOS_ERROR << "Check the value: " << value << " is in the correct format" << std::endl;
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
                KRATOS_ERROR << "Check the value: " << value << " is in the correct format" << std::endl;
            }
        } else if(KratosComponents<Variable<std::string> >::Has(variable_name)) {
            const Variable<std::string>& variable = KratosComponents<Variable<std::string>>().Get(variable_name);
            if (value.IsString()) {
                p_prop->SetValue(variable, value.GetString());
            } else {
                KRATOS_ERROR << "Check the value: " << value << " is in the correct format" << std::endl;
            }
        } else {
            KRATOS_ERROR << "Value type not defined";
        }
    }

    // Add / override tables in the p_properties
    Parameters tables = data["Material"]["Tables"];
    for(auto iter = tables.begin(); iter != tables.end(); iter++) {
        auto table_param = tables.GetValue(iter.name());
        // Case table is double, double. TODO(marandra): Does it make sense to consider other cases?
        Table<double> table;

        std::string input_var_name = table_param["input_variable"].GetString();
        TrimComponentName(input_var_name);
        std::string output_var_name = table_param["output_variable"].GetString();
        TrimComponentName(output_var_name);

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
