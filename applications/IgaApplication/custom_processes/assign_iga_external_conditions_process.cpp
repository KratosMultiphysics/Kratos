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
    Model& rModel, Parameters ThisParameters) : mpModel(&rModel), mParameters(ThisParameters)
{
    ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());

    KRATOS_ERROR_IF_NOT(mParameters.Has("model_part_name"))
        << "Missing \"model_part_name\" section" << std::endl;

    mElementConditionList = mParameters["element_condition_list"];
}

void AssignIgaExternalConditionsProcess::ExecuteInitialize(){

    ModelPart& r_model_part = mpModel->GetModelPart(mParameters["model_part_name"].GetString());

    const std::vector<std::string> sub_model_part_names = r_model_part.GetSubModelPartNames();

    // Number of element_condition_list (Elements or Conditions), in the end are all the Gauss Points
    const SizeType number_of_entitites = mElementConditionList.size();  

    const ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    double t = r_process_info[TIME];
    
    for (IndexType i = 0; i < number_of_entitites; i++) {
        std::string sub_model_part_name = mElementConditionList[i]["iga_model_part"].GetString();
        // Get the sub_model_part. Example "SBM_Support_outer"
        ModelPart& sub_model_part = r_model_part.GetSubModelPart(sub_model_part_name);

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
                    
                    if (sub_model_part.Conditions().size() > 0) {
                        for (auto i_cond = sub_model_part.ConditionsBegin(); i_cond != sub_model_part.ConditionsEnd(); i_cond++) {
                            const double x = i_cond->GetGeometry().Center().X();
                            const double y = i_cond->GetGeometry().Center().Y();
                            const double z = i_cond->GetGeometry().Center().Z();

                            const double dir_x_value = func_dir_x->CallFunction(x, y, z, t);
                            const double dir_y_value = func_dir_y->CallFunction(x, y, z, t);
                            const double dir_z_value = func_dir_z->CallFunction(x, y, z, t);
                            const double module_value = func_module->CallFunction(x, y, z, t);
                            
                            i_cond->SetValue(DIRECTION_X, dir_x_value);
                            i_cond->SetValue(DIRECTION_Y, dir_y_value);
                            i_cond->SetValue(DIRECTION_Z, dir_z_value);
                            i_cond->SetValue(MODULUS, module_value);
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

                    SetExternalConditionToElementsAndConditions(sub_model_part, variable_name, variable_component_name, values_string, 
                                                                r_process_info, t);
                }
            }
        }
    }
}

void AssignIgaExternalConditionsProcess::ExecuteInitializeSolutionStep(){

    ModelPart& r_model_part = mpModel->GetModelPart(mParameters["model_part_name"].GetString());

    const std::vector<std::string> sub_model_part_names = r_model_part.GetSubModelPartNames();

    // Number of element_condition_list (Elements or Conditions), in the end are all the Gauss Points
    const SizeType number_of_entitites = mElementConditionList.size();  

    const ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    double t = r_process_info[TIME];
    
    for (IndexType i = 0; i < number_of_entitites; i++) {
        std::string sub_model_part_name = mElementConditionList[i]["iga_model_part"].GetString();
        // Get the sub_model_part. Example "SBM_Support_outer"
        ModelPart& sub_model_part = r_model_part.GetSubModelPart(sub_model_part_name);

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

                    std::string direction_x_string = parameters[i_var]["direction"][0].GetString();
                    std::string direction_y_string = parameters[i_var]["direction"][1].GetString();
                    std::string direction_z_string = parameters[i_var]["direction"][2].GetString();
                    std::string module_string = parameters[i_var]["modulus"].GetString();

                    Kratos::unique_ptr<GenericFunctionUtility> func_dir_x = Kratos::make_unique<GenericFunctionUtility>(direction_x_string);
                    Kratos::unique_ptr<GenericFunctionUtility> func_dir_y = Kratos::make_unique<GenericFunctionUtility>(direction_y_string);
                    Kratos::unique_ptr<GenericFunctionUtility> func_dir_z = Kratos::make_unique<GenericFunctionUtility>(direction_z_string);
                    Kratos::unique_ptr<GenericFunctionUtility> func_module = Kratos::make_unique<GenericFunctionUtility>(module_string);
                    
                    if (sub_model_part.Conditions().size() > 0) {
                        for (auto i_cond = sub_model_part.ConditionsBegin(); i_cond != sub_model_part.ConditionsEnd(); i_cond++) {
                            const double x = i_cond->GetGeometry().Center().X();
                            const double y = i_cond->GetGeometry().Center().Y();
                            const double z = i_cond->GetGeometry().Center().Z();

                            const double dir_x_value = func_dir_x->CallFunction(x, y, z, t);
                            const double dir_y_value = func_dir_y->CallFunction(x, y, z, t);
                            const double dir_z_value = func_dir_z->CallFunction(x, y, z, t);
                            const double module_value = func_module->CallFunction(x, y, z, t);
                            
                            i_cond->SetValue(DIRECTION_X, dir_x_value);
                            i_cond->SetValue(DIRECTION_Y, dir_y_value);
                            i_cond->SetValue(DIRECTION_Z, dir_z_value);
                            i_cond->SetValue(MODULUS, module_value);
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
                    SetExternalConditionToElementsAndConditions(sub_model_part, variable_name, variable_component_name, values_string, 
                                                                r_process_info, t);
                }
            }
        }
    }
}


void AssignIgaExternalConditionsProcess::SetExternalConditionToElementsAndConditions(
    ModelPart& rSubModelPart,
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
        if (rSubModelPart.Elements().size() > 0) {
            for (auto i_element = rSubModelPart.ElementsBegin(); i_element != rSubModelPart.ElementsEnd(); i_element++) {
                const double x = i_element->GetGeometry().Center().X();
                const double y = i_element->GetGeometry().Center().Y();
                const double z = i_element->GetGeometry().Center().Z();

                const double value = eval_function->CallFunction(x, y, z, Time);
                if (component_rVariableName == "HEAT_FLUX") {
                    i_element->SetValue(HEAT_FLUX, value);
                } else if (component_rVariableName == "VELOCITY_X") {
                    i_element->SetValue(VELOCITY_X, value);
                } else if (component_rVariableName == "VELOCITY_Y") {
                    i_element->SetValue(VELOCITY_Y, value);
                } else if (component_rVariableName == "VELOCITY_Z") {
                    i_element->SetValue(VELOCITY_Z, value);
                } else if (component_rVariableName == "BODY_FORCE_X") {
                    i_element->SetValue(BODY_FORCE_X, value);
                } else if (component_rVariableName == "BODY_FORCE_Y") {
                    i_element->SetValue(BODY_FORCE_Y, value);
                } else if (component_rVariableName == "BODY_FORCE_Z") {
                    i_element->SetValue(BODY_FORCE_Z, value);
                } else {
                    KRATOS_ERROR << "Variable with name " << rVariableName << " not found or defined "
                    << "in the AssignIgaExternalConditionsProcess" << std::endl;
                }
            }
        }
        else if (rSubModelPart.Conditions().size() > 0) {
            for (auto i_cond = rSubModelPart.ConditionsBegin(); i_cond != rSubModelPart.ConditionsEnd(); i_cond++) {
                const double x = i_cond->GetGeometry().Center().X();
                const double y = i_cond->GetGeometry().Center().Y();
                const double z = i_cond->GetGeometry().Center().Z();
                const double value = eval_function->CallFunction(x, y, z, Time);
                if (component_rVariableName == "TEMPERATURE") {
                    i_cond->SetValue(TEMPERATURE, value);
                } else if (component_rVariableName == "FACE_HEAT_FLUX") {
                    i_cond->SetValue(FACE_HEAT_FLUX, value);
                } else if (component_rVariableName == "VELOCITY_X") {
                    i_cond->SetValue(VELOCITY_X, value);
                } else if (component_rVariableName == "VELOCITY_Y") {
                    i_cond->SetValue(VELOCITY_Y, value);
                } else if (component_rVariableName == "VELOCITY_Z") {
                    i_cond->SetValue(VELOCITY_Z, value);
                } else if (component_rVariableName == "DISPLACEMENT_X") {
                    i_cond->SetValue(DISPLACEMENT_X, value);
                } else if (component_rVariableName == "DISPLACEMENT_Y") {
                    i_cond->SetValue(DISPLACEMENT_Y, value);
                } else if (component_rVariableName == "DISPLACEMENT_Z") {
                    i_cond->SetValue(DISPLACEMENT_Z, value);
                } else if (component_rVariableName == "FORCE_X") {
                    i_cond->SetValue(FORCE_X, value);
                } else if (component_rVariableName == "FORCE_Y") {
                    i_cond->SetValue(FORCE_Y, value);
                } else if (component_rVariableName == "FORCE_Z") {
                    i_cond->SetValue(FORCE_Z, value);
                } else {
                    KRATOS_ERROR << "No name found" ;
                }
            } 
        } else {KRATOS_ERROR << "AssignIgaExternalConditionsProcess : No Condition or Elements defined" ;}
    }
}
} // End namespace Kratos
