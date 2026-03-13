//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Andrea Gorgi, 
//                  Polytimi Zisimoupoulu

// Project includes
#include "geometries/nurbs_volume_geometry.h"

// Application includes
#include "assign_iga_external_conditions_process.h"
#include "iga_application_variables.h"

namespace Kratos
{

AssignIgaExternalConditionsProcess::AssignIgaExternalConditionsProcess(
    Model& rModel,
    Parameters ThisParameters)
    : mpModel(&rModel)
    , mParameters(ThisParameters)
{
    ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());

    KRATOS_ERROR_IF_NOT(mParameters.Has("model_part_name"))
        << "Missing \"model_part_name\" section" << std::endl;

    mElementConditionList = mParameters["element_condition_list"];

    ExpandElementConditionList();
}

void AssignIgaExternalConditionsProcess::ExecuteInitialize(){

    // Number of element_condition_list (Elements or Conditions), in the end are all the Gauss Points
    const SizeType number_of_entitites = mElementConditionList.size();
    
    for (IndexType i = 0; i < number_of_entitites; i++) {
        std::string model_part_name = mElementConditionList[i]["iga_model_part"].GetString();
        // Get the sub_model_part. Example "SBM_Support_outer"
        ModelPart& r_model_part = mpModel->GetModelPart(model_part_name);

        const ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        double t = r_process_info[TIME];

        if (mElementConditionList[i].Has("variables")) {
            const auto& parameters = mElementConditionList[i]["variables"];
            
            for (IndexType i_var = 0; i_var < parameters.size(); i_var++){
                
                // set the variable as static by deafault
                parameters[i_var].AddEmptyValue("is_static").SetBool(true);

                // external conditions by direction
                //--------------------------------------------------------------------------------------------------------------------------------------
                bool assign_by_direction = parameters[i_var].Has("assign_by_direction");
                if (assign_by_direction) assign_by_direction = parameters[i_var]["assign_by_direction"].GetBool();

                if (assign_by_direction)
                {
                    // check input
                    if (!(parameters[i_var].Has("direction") && parameters[i_var].Has("modulus"))) 
                        KRATOS_ERROR << "ERROR in assign_iga_external_conditions_process. ""direction"" or ""modulus"" not defined in assign_by_direction. \n";

                    KRATOS_ERROR_IF_NOT(parameters[i_var].Has("variable_name"))
                        << "ERROR in assign_iga_external_conditions_process. \"variable_name\" not defined in assign_by_direction. \n";

                    std::string variable_name = parameters[i_var]["variable_name"].GetString();

                    std::string direction_x_string = parameters[i_var]["direction"][0].GetString();
                    std::string direction_y_string = parameters[i_var]["direction"][1].GetString();
                    std::string direction_z_string = parameters[i_var]["direction"][2].GetString();
                    std::string module_string = parameters[i_var]["modulus"].GetString();

                    if (direction_x_string.find("t") != std::string::npos || direction_y_string.find("t") != std::string::npos ||
                        direction_z_string.find("t") != std::string::npos || module_string.find("t") != std::string::npos) {
                        // If the direction or the module depends on time, then the variable is not static
                        parameters[i_var]["is_static"].SetBool(false);
                    }

                    Kratos::unique_ptr<GenericFunctionUtility> func_dir_x = Kratos::make_unique<GenericFunctionUtility>(direction_x_string);
                    Kratos::unique_ptr<GenericFunctionUtility> func_dir_y = Kratos::make_unique<GenericFunctionUtility>(direction_y_string);
                    Kratos::unique_ptr<GenericFunctionUtility> func_dir_z = Kratos::make_unique<GenericFunctionUtility>(direction_z_string);
                    Kratos::unique_ptr<GenericFunctionUtility> func_module = Kratos::make_unique<GenericFunctionUtility>(module_string);
                    
                    if (r_model_part.Conditions().size() > 0) {
                        for (auto i_cond = r_model_part.ConditionsBegin(); i_cond != r_model_part.ConditionsEnd(); i_cond++) {
                            const double x = i_cond->GetGeometry().Center().X();
                            const double y = i_cond->GetGeometry().Center().Y();
                            const double z = i_cond->GetGeometry().Center().Z();

                            const double dir_x_value = func_dir_x->CallFunction(x, y, z, t);
                            const double dir_y_value = func_dir_y->CallFunction(x, y, z, t);
                            const double dir_z_value = func_dir_z->CallFunction(x, y, z, t);
                            const double module_value = func_module->CallFunction(x, y, z, t);

                            Vector dirichlet_vector(3);
                            dirichlet_vector[0] = module_value*dir_x_value;
                            dirichlet_vector[1] = module_value*dir_y_value;
                            dirichlet_vector[2] = module_value*dir_z_value;
                            
                            i_cond->SetValue(DIRECTION_X, dir_x_value);
                            i_cond->SetValue(DIRECTION_Y, dir_y_value);
                            i_cond->SetValue(DIRECTION_Z, dir_z_value);

                            Condition::Pointer p_condition = *i_cond.base();

                            std::vector<std::string> components_array = {"_X","_Y","_Z"};
                            for (IndexType i_component = 0; i_component < 3; i_component++ ) {
                                std::string component_variable_name = variable_name+components_array[i_component];
                                
                                // Set the value to the condition
                                SetVariableValueToCondition(component_variable_name, dirichlet_vector[i_component], p_condition);
                            }
                        } 
                    } else {KRATOS_ERROR << "AssignIgaExternalConditionsProcess : No Conditions defined for AssignByDirection" ;}
                }
                //-------------------------------------------------------------------------------------------------------------------------------------
                else // standard external condition (in all directions)
                {

                    // Get the variable. Example "HEAT_FLUX"
                    std::string variable_name = parameters[i_var]["variable_name"].GetString();
                    std::vector<std::string> variable_component_name;
                    std::vector<std::string> values_string;
                    std::vector<std::string> components_array = {"_X","_Y","_Z"};
                    bool isScalarValue = parameters[i_var]["value"].IsString();
                    if (isScalarValue == 1 ) {
                        values_string.push_back(parameters[i_var]["value"].GetString()); 
                        variable_component_name.push_back(variable_name) ;
                    } else {
                        SizeType number_of_components = parameters[i_var]["value"].size();
                        for (IndexType i_component = 0; i_component < number_of_components; i_component++ ) {
                            values_string.push_back(parameters[i_var]["value"][i_component].GetString());
                            std::string component_variable_name = variable_name+components_array[i_component];
                            variable_component_name.push_back(component_variable_name) ;
                        }  
                    }  

                    // check if the variable is static or not
                    for (const auto& str : values_string) {
                        if (str.find("t") != std::string::npos) {
                            // the variable is not static
                            parameters[i_var]["is_static"].SetBool(false);
                            break;  
                        }
                    }

                    SetExternalConditionToElementsAndConditions(r_model_part, variable_name, variable_component_name, values_string, 
                                                                r_process_info, t);
                }
            }
        }
    }
}

void AssignIgaExternalConditionsProcess::ExecuteInitializeSolutionStep(){

    const SizeType number_of_entitites = mElementConditionList.size();  
    
    for (IndexType i = 0; i < number_of_entitites; i++) {
        std::string model_part_name = mElementConditionList[i]["iga_model_part"].GetString();
        // Get the sub_model_part. Example "SBM_Support_outer"
        ModelPart& r_model_part = mpModel->GetModelPart(model_part_name);

        const ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        double t = r_process_info[TIME];

        if (mElementConditionList[i].Has("variables")) {
            auto parameters = mElementConditionList[i]["variables"];
            for (IndexType i_var = 0; i_var < parameters.size(); i_var++){

                if (parameters[i_var].Has("is_static")) 
                    if(parameters[i_var]["is_static"].GetBool()) {
                        // Skip if the variable is static. It has already been assigned to the model part
                        continue;
                    }

                // external conditions by direction
                //--------------------------------------------------------------------------------------------------------------------------------------
                bool assign_by_direction = parameters[i_var].Has("assign_by_direction");
                if (assign_by_direction) assign_by_direction = parameters[i_var]["assign_by_direction"].GetBool();

                if (assign_by_direction)
                {
                    // check input
                    if (!(parameters[i_var].Has("direction") && parameters[i_var].Has("modulus"))) 
                        KRATOS_ERROR << "ERROR in assign_iga_external_conditions_process. ""direction"" or ""modulus"" not defined in assign_by_direction. \n";

                    KRATOS_ERROR_IF_NOT(parameters[i_var].Has("variable_name"))
                        << "ERROR in assign_iga_external_conditions_process. \"variable_name\" not defined in assign_by_direction. \n";

                    std::string variable_name = parameters[i_var]["variable_name"].GetString();

                    std::string direction_x_string = parameters[i_var]["direction"][0].GetString();
                    std::string direction_y_string = parameters[i_var]["direction"][1].GetString();
                    std::string direction_z_string = parameters[i_var]["direction"][2].GetString();
                    std::string module_string = parameters[i_var]["modulus"].GetString();

                    Kratos::unique_ptr<GenericFunctionUtility> func_dir_x = Kratos::make_unique<GenericFunctionUtility>(direction_x_string);
                    Kratos::unique_ptr<GenericFunctionUtility> func_dir_y = Kratos::make_unique<GenericFunctionUtility>(direction_y_string);
                    Kratos::unique_ptr<GenericFunctionUtility> func_dir_z = Kratos::make_unique<GenericFunctionUtility>(direction_z_string);
                    Kratos::unique_ptr<GenericFunctionUtility> func_module = Kratos::make_unique<GenericFunctionUtility>(module_string);
                    
                    if (r_model_part.Conditions().size() > 0) {
                        for (auto i_cond = r_model_part.ConditionsBegin(); i_cond != r_model_part.ConditionsEnd(); i_cond++) {
                            const double x = i_cond->GetGeometry().Center().X();
                            const double y = i_cond->GetGeometry().Center().Y();
                            const double z = i_cond->GetGeometry().Center().Z();

                            const double dir_x_value = func_dir_x->CallFunction(x, y, z, t);
                            const double dir_y_value = func_dir_y->CallFunction(x, y, z, t);
                            const double dir_z_value = func_dir_z->CallFunction(x, y, z, t);
                            const double module_value = func_module->CallFunction(x, y, z, t);
                            
                            Vector dirichlet_vector(3);
                            dirichlet_vector[0] = module_value*dir_x_value;
                            dirichlet_vector[1] = module_value*dir_y_value;
                            dirichlet_vector[2] = module_value*dir_z_value;
                            
                            i_cond->SetValue(DIRECTION_X, dir_x_value);
                            i_cond->SetValue(DIRECTION_Y, dir_y_value);
                            i_cond->SetValue(DIRECTION_Z, dir_z_value);

                            Condition::Pointer p_condition = *i_cond.base();

                            std::vector<std::string> components_array = {"_X","_Y","_Z"};
                            for (IndexType i_component = 0; i_component < 3; i_component++ ) {
                                std::string component_variable_name = variable_name+components_array[i_component];
                                
                                // Set the value to the condition
                                SetVariableValueToCondition(component_variable_name, dirichlet_vector[i_component], p_condition);
                            }
                        } 
                    } else {KRATOS_ERROR << "AssignIgaExternalConditionsProcess : No Conditions defined for AssignByDirection" ;}
                }
                //-------------------------------------------------------------------------------------------------------------------------------------
                else // standard external condition (in all directions)
                {

                    // Get the variable. Example "HEAT_FLUX"
                    std::string variable_name = parameters[i_var]["variable_name"].GetString();
                    std::vector<std::string> variable_component_name;
                    std::vector<std::string> values_string;
                    std::vector<std::string> components_array = {"_X","_Y","_Z"};
                    bool isScalarValue = parameters[i_var]["value"].IsString();
                    if (isScalarValue == 1 ) {
                        values_string.push_back(parameters[i_var]["value"].GetString()); 
                        variable_component_name.push_back(variable_name) ;
                    } else {
                        SizeType number_of_components = parameters[i_var]["value"].size();
                        for (IndexType i_component = 0; i_component < number_of_components; i_component++ ) {
                            values_string.push_back(parameters[i_var]["value"][i_component].GetString());
                            std::string component_variable_name = variable_name+components_array[i_component];
                            variable_component_name.push_back(component_variable_name) ;
                        }  
                    }  
                    SetExternalConditionToElementsAndConditions(r_model_part, variable_name, variable_component_name, values_string, 
                                                                r_process_info, t);
                }
            }
        }
    }
}


void AssignIgaExternalConditionsProcess::SetExternalConditionToElementsAndConditions(
    ModelPart& rModelPart,
    const std::string& rVariableName,
    const std::vector<std::string>& rVariableComponentName,
    const std::vector<std::string>& rValuesString,
    const ProcessInfo& rProcessInfo,
    const double Time)
{
    for (IndexType i_component = 0; i_component < rValuesString.size(); i_component++) {   
        // Define the function to be parsed      
        Kratos::unique_ptr<GenericFunctionUtility> eval_function = Kratos::make_unique<GenericFunctionUtility>(rValuesString[i_component]);
        const std::string component_rVariableName = rVariableComponentName[i_component];
        if (rModelPart.Elements().size() > 0) {
            for (auto i_element = rModelPart.ElementsBegin(); i_element != rModelPart.ElementsEnd(); i_element++) {
                const double x = i_element->GetGeometry().Center().X();
                const double y = i_element->GetGeometry().Center().Y();
                const double z = i_element->GetGeometry().Center().Z();
                const double value = eval_function->CallFunction(x, y, z, Time);

                Element::Pointer p_element = *i_element.base();
                SetVariableValueToElement(component_rVariableName, value, p_element);
            }
        }
        else if (rModelPart.Conditions().size() > 0) {
            for (auto i_cond = rModelPart.ConditionsBegin(); i_cond != rModelPart.ConditionsEnd(); i_cond++) {
                const double x = i_cond->GetGeometry().Center().X();
                const double y = i_cond->GetGeometry().Center().Y();
                const double z = i_cond->GetGeometry().Center().Z();
                const double value = eval_function->CallFunction(x, y, z, Time);

                Condition::Pointer p_condition = *i_cond.base();
                SetVariableValueToCondition(component_rVariableName, value, p_condition);
            } 
        } 
        // else {KRATOS_ERROR << "AssignIgaExternalConditionsProcess : No Condition or Elements defined" ;}
    }
}

void AssignIgaExternalConditionsProcess::SetVariableValueToElement(
    const std::string& rVariableName,
    const double& value,
    Element::Pointer p_element)
{
    if (rVariableName == "HEAT_FLUX") {
        p_element->SetValue(HEAT_FLUX, value);
    } else if (rVariableName == "VELOCITY_X") {
        p_element->SetValue(VELOCITY_X, value);
    } else if (rVariableName == "VELOCITY_Y") {
        p_element->SetValue(VELOCITY_Y, value);
    } else if (rVariableName == "VELOCITY_Z") {
        p_element->SetValue(VELOCITY_Z, value);
    } else if (rVariableName == "BODY_FORCE_X") {
        p_element->SetValue(BODY_FORCE_X, value);
    } else if (rVariableName == "BODY_FORCE_Y") {
        p_element->SetValue(BODY_FORCE_Y, value);
    } else if (rVariableName == "BODY_FORCE_Z") {
        p_element->SetValue(BODY_FORCE_Z, value);
    } else {
        KRATOS_ERROR << "Variable with name " << rVariableName << " not found or defined "
        << "in the AssignIgaExternalConditionsProcess" << std::endl;
    }
} 

void AssignIgaExternalConditionsProcess::SetVariableValueToCondition(
    const std::string& rVariableName,
    const double& value,
    Condition::Pointer p_condition)
{
    if (rVariableName == "TEMPERATURE") {
        p_condition->SetValue(TEMPERATURE, value);
    } else if (rVariableName == "FACE_HEAT_FLUX") {
        p_condition->SetValue(FACE_HEAT_FLUX, value);
    } else if (rVariableName == "VELOCITY_X") {
        p_condition->SetValue(VELOCITY_X, value);
    } else if (rVariableName == "VELOCITY_Y") {
        p_condition->SetValue(VELOCITY_Y, value);
    } else if (rVariableName == "VELOCITY_Z") {
        p_condition->SetValue(VELOCITY_Z, value);
    } else if (rVariableName == "PRESSURE") {
        p_condition->SetValue(PRESSURE, value);
    } else if (rVariableName == "DISPLACEMENT_X") {
        p_condition->SetValue(DISPLACEMENT_X, value);
    } else if (rVariableName == "DISPLACEMENT_Y") {
        p_condition->SetValue(DISPLACEMENT_Y, value);
    } else if (rVariableName == "DISPLACEMENT_Z") {
        p_condition->SetValue(DISPLACEMENT_Z, value);
    } else if (rVariableName == "FORCE_X") {
        p_condition->SetValue(FORCE_X, value);
    } else if (rVariableName == "FORCE_Y") {
        p_condition->SetValue(FORCE_Y, value);
    } else if (rVariableName == "FORCE_Z") {
        p_condition->SetValue(FORCE_Z, value);
    } else {
        KRATOS_ERROR << "No name found" ;
    }
} 

void AssignIgaExternalConditionsProcess::ExpandElementConditionList()
{
    Parameters expanded("[]");

    const std::string root_name = mParameters["model_part_name"].GetString();
    ModelPart& r_root = mpModel->GetModelPart(root_name);

    for (IndexType i = 0; i < mElementConditionList.size(); ++i) {
        const Parameters item = mElementConditionList[i];

        const bool for_all = item.Has("apply_to_all_patches") && item["apply_to_all_patches"].GetBool();
        const std::string suffix = item.Has("iga_model_part_suffix") ? item["iga_model_part_suffix"].GetString() : "";
        const std::string patch_prefix = item.Has("patch_prefix") ? item["patch_prefix"].GetString() : "Patch";

        if (for_all) {
            KRATOS_ERROR_IF(suffix.empty())
                << "AssignIgaExternalConditionsProcess: 'apply_to_all_patches' requires 'iga_model_part_suffix'." << std::endl;

            for (auto& rPatch : r_root.SubModelParts()) {
                const std::string& patch_name = rPatch.Name();                 // <-- define it
                if (patch_name.rfind(patch_prefix, 0) != 0)                    // starts-with check
                    continue;

                const std::string full_target = rPatch.FullName() + "." + suffix;
                if (!mpModel->HasModelPart(full_target))
                    continue; // skip missing children

                Parameters clone = item.Clone();
                // materialize the concrete target and strip helper keys
                if (clone.Has("iga_model_part")) clone["iga_model_part"].SetString(full_target);
                else clone.AddEmptyValue("iga_model_part").SetString(full_target);

                if (clone.Has("apply_to_all_patches"))  clone.RemoveValue("apply_to_all_patches");
                if (clone.Has("iga_model_part_suffix")) clone.RemoveValue("iga_model_part_suffix");
                if (clone.Has("patch_prefix"))          clone.RemoveValue("patch_prefix");

                expanded.Append(clone);
            }
        } else {
            expanded.Append(item); // keep explicit entries
        }
    }

    mElementConditionList = expanded;
    // KRATOS_WATCH(mElementConditionList)
    // exit(0);
}


} // namespace Kratos
