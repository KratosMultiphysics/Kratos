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
#include "utilities/read_materials_utility.h"

namespace Kratos {
namespace {

template <class TValueType>
void CheckIfOverwritingValue(const Properties& rProps,
                             const Variable<TValueType>& rVariable,
                             const TValueType& rValue)
{
    KRATOS_WARNING_IF("ReadMaterialsUtility", rProps.Has(rVariable)) << "The properties ID: "
        << rProps.Id() << " already has " << rVariable.Name() << "\nOverwriting "
        << rProps[rVariable] << " with " << rValue << std::endl;
}

}

/***********************************************************************************/
/***********************************************************************************/

ReadMaterialsUtility::ReadMaterialsUtility(
    Parameters Params,
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

    Params.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Read json string in materials file, create Parameters
    const std::string& materials_filename = Params["Parameters"]["materials_filename"].GetString();
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

void ReadMaterialsUtility::ReadMaterials(Parameters MaterialData)
{
    GetPropertyBlock(MaterialData);
}

/***********************************************************************************/
/***********************************************************************************/

void ReadMaterialsUtility::GetPropertyBlock(Parameters Materials)
{
    KRATOS_INFO("Read materials") << "Started" << std::endl;

    CheckUniqueMaterialAssignment(Materials);

    for (IndexType i = 0; i < Materials["properties"].size(); ++i) {
        Parameters material = Materials["properties"].GetArrayItem(i);
        AssignPropertyBlock(material);
    }
    KRATOS_INFO("Read materials") << "Finished" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void ReadMaterialsUtility::TrimComponentName(std::string& rLine){
    std::stringstream ss(rLine);
    std::size_t counter = 0;
    while (std::getline(ss, rLine, '.')){++counter;}
    KRATOS_WARNING_IF("Read materials", counter > 1) << "Ignoring module information for component " << rLine << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void ReadMaterialsUtility::AssignPropertyBlock(Parameters Data)
{
    // Get the properties for the specified model part.
    ModelPart& r_model_part = mrModel.GetModelPart(Data["model_part_name"].GetString());
    const IndexType property_id = Data["properties_id"].GetInt();
    const IndexType mesh_id = 0;
    Properties::Pointer p_prop;
    if (r_model_part.RecursivelyHasProperties(property_id, mesh_id)) {
        KRATOS_WARNING("ReadMaterialsUtility") << "WARNING:: The properties ID: " << property_id
            << " in mesh ID: " << mesh_id << " is already defined. "
            << "This will overwrite the existing values" << std::endl;
        p_prop = r_model_part.pGetProperties(property_id, mesh_id);

        // Compute the size using the iterators
        std::size_t variables_size = 0;
        if (Data["Material"].Has("Variables")) {
            for(auto it=Data["Material"]["Variables"].begin(); it!=Data["Material"]["Variables"].end(); ++it) {
                ++variables_size;
            }
        }
        std::size_t tables_size = 0;
        if (Data["Material"].Has("Tables")) {
            for(auto it=Data["Material"]["Tables"].begin(); it!=Data["Material"]["Tables"].end(); ++it) {
                ++tables_size;
            }
        }

        KRATOS_WARNING_IF("ReadMaterialsUtility", variables_size > 0 && p_prop->HasVariables())
            << "WARNING:: The properties ID: " << property_id << " already has variables." << std::endl;
        KRATOS_WARNING_IF("ReadMaterialsUtility", tables_size > 0 && p_prop->HasTables())
            << "WARNING:: The properties ID: " << property_id << " already has tables." << std::endl;
    } else {
        p_prop = r_model_part.CreateNewProperties(property_id, mesh_id);
    }

    // Assign the p_properties to the model part's elements and conditions.
    auto& r_elements_array = r_model_part.Elements();
    auto& r_conditions_array = r_model_part.Conditions();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i) {
        auto it_elem = r_elements_array.begin() + i;
        it_elem->SetProperties(p_prop);
    }

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
        auto it_cond = r_conditions_array.begin() + i;
        it_cond->SetProperties(p_prop);
    }

    // Set the CONSTITUTIVE_LAW for the current p_properties.
    if (Data["Material"].Has("constitutive_law")) {
        Parameters cl_parameters = Data["Material"]["constitutive_law"];
        std::string constitutive_law_name = cl_parameters["name"].GetString();
        TrimComponentName(constitutive_law_name);
        cl_parameters["name"].SetString(constitutive_law_name);

        KRATOS_ERROR_IF_NOT(KratosComponents<ConstitutiveLaw>::Has(constitutive_law_name))
            << "Kratos components missing \"" << constitutive_law_name << "\"" << std::endl;
        auto p_constitutive_law =
            KratosComponents<ConstitutiveLaw>::Get(constitutive_law_name).Create(cl_parameters);
        p_prop->SetValue(CONSTITUTIVE_LAW, p_constitutive_law);
    } else {
        KRATOS_INFO("Read materials") << "No constitutive law defined for material ID: " << property_id << std::endl;
    }

    // Add / override the values of material parameters in the p_properties
    if (Data["Material"].Has("Variables")) {
        Parameters variables = Data["Material"]["Variables"];
        for (auto iter = variables.begin(); iter != variables.end(); ++iter) {
            const Parameters value = variables.GetValue(iter.name());

            std::string variable_name = iter.name();
            TrimComponentName(variable_name);

            // We don't just copy the values, we do some tyransformation depending of the destination variable
            if (KratosComponents<Variable<double> >::Has(variable_name)) {
                const Variable<double>& r_variable = KratosComponents<Variable<double>>().Get(variable_name);
                CheckIfOverwritingValue(*p_prop, r_variable, value.GetDouble());
                p_prop->SetValue(r_variable, value.GetDouble());
            } else if(KratosComponents<Variable<bool> >::Has(variable_name)) {
                const Variable<bool>& r_variable = KratosComponents<Variable<bool>>().Get(variable_name);
                CheckIfOverwritingValue(*p_prop, r_variable, value.GetBool());
                p_prop->SetValue(r_variable, value.GetBool());
            } else if(KratosComponents<Variable<int> >::Has(variable_name)) {
                const Variable<int>& r_variable = KratosComponents<Variable<int>>().Get(variable_name);
                CheckIfOverwritingValue(*p_prop, r_variable, value.GetInt());
                p_prop->SetValue(r_variable, value.GetInt());
            } else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name)) {
                const Variable<array_1d<double, 3>>& r_variable = KratosComponents<Variable<array_1d<double, 3>>>().Get(variable_name);
                array_1d<double, 3> temp = ZeroVector(3);
                const Vector& r_value_variable = value.GetVector();
                KRATOS_ERROR_IF(r_value_variable.size() != 3) << "The vector of variable " << variable_name << " has size " << r_value_variable.size() << " and it is supposed to be 3" << std::endl;
                for (IndexType index = 0; index < 3; ++index)
                    temp[index] = r_value_variable[index];
                CheckIfOverwritingValue(*p_prop, r_variable, temp);
                p_prop->SetValue(r_variable, temp);
            } else if(KratosComponents<Variable<array_1d<double, 6> > >::Has(variable_name)) {
                const Variable<array_1d<double, 6>>& r_variable = KratosComponents<Variable<array_1d<double, 6>>>().Get(variable_name);
                array_1d<double, 6> temp(6, 0.0);
                const Vector& r_value_variable = value.GetVector();
                KRATOS_ERROR_IF(r_value_variable.size() != 6) << "The vector of variable " << variable_name << " has size " << r_value_variable.size() << " and it is supposed to be 6" << std::endl;
                for (IndexType index = 0; index < 6; ++index)
                    temp[index] = r_value_variable[index];
                CheckIfOverwritingValue(*p_prop, r_variable, temp);
                p_prop->SetValue(r_variable, temp);
            } else if(KratosComponents<Variable<Vector > >::Has(variable_name)) {
                const Variable<Vector>& r_variable = KratosComponents<Variable<Vector>>().Get(variable_name);
                CheckIfOverwritingValue(*p_prop, r_variable, value.GetVector());
                p_prop->SetValue(r_variable, value.GetVector());
            } else if(KratosComponents<Variable<Matrix> >::Has(variable_name)) {
                const Variable<Matrix>& r_variable = KratosComponents<Variable<Matrix>>().Get(variable_name);
                CheckIfOverwritingValue(*p_prop, r_variable, value.GetMatrix());
                p_prop->SetValue(r_variable, value.GetMatrix());
            } else if(KratosComponents<Variable<std::string> >::Has(variable_name)) {
                const Variable<std::string>& r_variable = KratosComponents<Variable<std::string>>().Get(variable_name);
                CheckIfOverwritingValue(*p_prop, r_variable, value.GetString());
                p_prop->SetValue(r_variable, value.GetString());
            } else {
                KRATOS_ERROR << "Value type for \"" << variable_name << "\" not defined";
            }
        }
    } else {
        KRATOS_INFO("Read materials") << "No variables defined for material ID: " << property_id << std::endl;
    }

    // Add / override tables in the p_properties
    if (Data["Material"].Has("Tables")) {
        Parameters tables = Data["Material"]["Tables"];
        for (auto iter = tables.begin(); iter != tables.end(); ++iter) {
            auto table_param = tables.GetValue(iter.name());
            // Case table is double, double. TODO(marandra): Does it make sense to consider other cases?
            Table<double> table;

            std::string input_var_name = table_param["input_variable"].GetString();
            TrimComponentName(input_var_name);
            std::string output_var_name = table_param["output_variable"].GetString();
            TrimComponentName(output_var_name);

            const auto& r_input_var = KratosComponents<Variable<double>>().Get(input_var_name);
            const auto& r_output_var = KratosComponents<Variable<double>>().Get(output_var_name);
            for (IndexType i = 0; i < table_param["data"].size(); ++i) {
                table.insert(table_param["data"][i][0].GetDouble(),
                            table_param["data"][i][1].GetDouble());
            }
            p_prop->SetTable(r_input_var, r_output_var, table);
        }
    } else {
        KRATOS_INFO("Read materials") << "No tables defined for material ID: " << property_id << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ReadMaterialsUtility::CheckUniqueMaterialAssignment(Parameters Materials)
{
    const std::size_t num_props = Materials["properties"].size();

    // save all ModelPartNames in a vector
    std::vector<std::string> model_part_names(num_props);
    for (IndexType i = 0; i < num_props; ++i) {
        model_part_names[i] = Materials["properties"].GetArrayItem(i)["model_part_name"].GetString();
    }

    // sort the names
    std::sort(model_part_names.begin(), model_part_names.end());

    // check if the same name exists multiple times (this requires the sorting)
    const auto it = std::unique( model_part_names.begin(), model_part_names.end() );
    KRATOS_ERROR_IF_NOT(it == model_part_names.end()) << "Materials for ModelPart \""
        << *it << "\" are specified multiple times!" << std::endl;

    // checking if a parent also has a materials definition, i.e. if the assignment is unique
    std::string parent_model_part_name;
    for (IndexType i = 0; i < num_props; ++i) {
        parent_model_part_name = model_part_names[i];

        // removing the submodelpart-names one-by-one
        while (parent_model_part_name.find(".") != std::string::npos) {
            std::size_t found_pos = parent_model_part_name.find_last_of(".");
            parent_model_part_name = parent_model_part_name.substr(0, found_pos);

            // check if the parent-modelpart-name also has a materials definition
            const bool parent_has_materials = std::find(model_part_names.begin(), model_part_names.end(),
                parent_model_part_name) != model_part_names.end();

            KRATOS_ERROR_IF(parent_has_materials) << "Materials for ModelPart \""
                << model_part_names[i] << "\" are specified multiple times!\n"
                << "Overdefined due to also specifying the materials for Parent-ModelPart \""
                << parent_model_part_name << "\"!" << std::endl;
        }
    }
}

}  // namespace Kratos.
