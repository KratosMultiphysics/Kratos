//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Andrea Gorgi, Polytimi Zisimoupoulu

// Project includes
#include "assign_iga_external_conditions_process.h"
#include "geometries/nurbs_volume_geometry.h"
#include "iga_application_variables.h"

namespace Kratos
{

AssignIgaExternalConditionsProcess::AssignIgaExternalConditionsProcess(
    Model& rModel, Parameters ThisParameters) : mpModel(&rModel), mParameters(ThisParameters)
{
    ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());

    mEchoLevel = mParameters["echo_level"].GetInt();

    KRATOS_ERROR_IF_NOT(mParameters.Has("analysis_model_part_name"))
        << "Missing \"analysis_model_part_name\" section" << std::endl;
    ModelPart& analysis_model_part = mpModel->GetModelPart(mParameters["analysis_model_part_name"].GetString());

    mIgaPhysicsParameters = mParameters["element_condition_list"];
}

void AssignIgaExternalConditionsProcess::Execute(){

    ModelPart& analysis_model_part = mpModel->GetModelPart(mParameters["analysis_model_part_name"].GetString());

    const std::vector<std::string> sub_model_part_names = analysis_model_part.GetSubModelPartNames();

    // Number of element_condition_list (Elements or Conditions), in the end are all the Gauss Points
    const int number_of_entitites = mIgaPhysicsParameters.size();  

    const ProcessInfo& r_process_info = analysis_model_part.GetProcessInfo();
    double t = r_process_info[TIME];
    
    for (IndexType i = 0; i < number_of_entitites; i++) {
        std::string sub_model_part_name = mIgaPhysicsParameters[i]["iga_model_part"].GetString();
        // Get the sub_model_part. Example "SBM_Support_outer"
        ModelPart& sub_model_part = analysis_model_part.GetSubModelPart(sub_model_part_name);

        if (mIgaPhysicsParameters[i].Has("variables")) {
            auto parameters = mIgaPhysicsParameters[i]["variables"];
            for (IndexType i_var = 0; i_var < parameters.size(); i_var++){
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
                            i_cond->SetValue(MODULE, module_value);
                        } 
                    } else {KRATOS_ERROR << "AssignIgaExternalConditionsProcess : No Conditions defined for AssignByDirection" ;}
                }
                //-------------------------------------------------------------------------------------------------------------------------------------
                else // standard external condition (in all directions)
                {

                    // Get the variable. Example "HEAT_FLUX"
                    std::string variable_name = parameters[i_var]["variable_name"].GetString();
                    std::vector<std::string> variable_name_array;
                    std::vector<std::string> values_string;
                    std::vector<std::string> components_array = {"_X","_Y","_Z"};
                    bool isScalarValue = parameters[i_var]["value"].IsString();
                    if (isScalarValue == 1 ) {
                        values_string.push_back(parameters[i_var]["value"].GetString()); 
                        variable_name_array.push_back(variable_name) ;
                    } else {
                        const int numberOfComponents = parameters[i_var]["value"].size();
                        for (IndexType i_component = 0; i_component < numberOfComponents; i_component++ ) {
                            values_string.push_back(parameters[i_var]["value"][i_component].GetString());
                            std::string component_variable_name = variable_name+components_array[i_component];
                            variable_name_array.push_back(component_variable_name) ;
                        }  
                    }  
                    for (IndexType i_component = 0; i_component < values_string.size(); i_component++) {   
                        // Define the function to be parsed      
                        Kratos::unique_ptr<GenericFunctionUtility> eval_function = Kratos::make_unique<GenericFunctionUtility>(values_string[i_component]);
                        const std::string component_variable_name = variable_name_array[i_component];
                        if (sub_model_part.Elements().size() > 0) {
                            for (auto i_element = sub_model_part.ElementsBegin(); i_element != sub_model_part.ElementsEnd(); i_element++) {
                                const double x = i_element->GetGeometry().Center().X();
                                const double y = i_element->GetGeometry().Center().Y();
                                const double z = i_element->GetGeometry().Center().Z();


                                const double value = eval_function->CallFunction(x, y, z, t);
                                if (component_variable_name == "HEAT_FLUX") {
                                    i_element->SetValue(HEAT_FLUX, value);
                                } else if (component_variable_name == "VELOCITY_X") {
                                    i_element->SetValue(VELOCITY_X, value);
                                } else if (component_variable_name == "VELOCITY_Y") {
                                    i_element->SetValue(VELOCITY_Y, value);
                                } else if (component_variable_name == "VELOCITY_Z") {
                                    i_element->SetValue(VELOCITY_Z, value);
                                } else if (component_variable_name == "BODY_FORCE_X") {
                                    i_element->SetValue(BODY_FORCE_X, value);
                                } else if (component_variable_name == "BODY_FORCE_Y") {
                                    i_element->SetValue(BODY_FORCE_Y, value);
                                } else if (component_variable_name == "BODY_FORCE_Z") {
                                    i_element->SetValue(BODY_FORCE_Z, value);
                                } else {
                                    KRATOS_ERROR << "No name found" ;
                                }
                            }
                        }
                        else if (sub_model_part.Conditions().size() > 0) {
                            for (auto i_cond = sub_model_part.ConditionsBegin(); i_cond != sub_model_part.ConditionsEnd(); i_cond++) {
                                const double x = i_cond->GetGeometry().Center().X();
                                const double y = i_cond->GetGeometry().Center().Y();
                                const double z = i_cond->GetGeometry().Center().Z();
                                const double value = eval_function->CallFunction(x, y, z, t);
                                if (component_variable_name == "TEMPERATURE") {
                                    i_cond->SetValue(TEMPERATURE, value);
                                } else if (component_variable_name == "FACE_HEAT_FLUX") {
                                    i_cond->SetValue(FACE_HEAT_FLUX, value);
                                } else if (component_variable_name == "VELOCITY_X") {
                                    i_cond->SetValue(VELOCITY_X, value);
                                } else if (component_variable_name == "VELOCITY_Y") {
                                    i_cond->SetValue(VELOCITY_Y, value);
                                } else if (component_variable_name == "VELOCITY_Z") {
                                    i_cond->SetValue(VELOCITY_Z, value);
                                } else if (component_variable_name == "DISPLACEMENT_X") {
                                    i_cond->SetValue(DISPLACEMENT_X, value);
                                } else if (component_variable_name == "DISPLACEMENT_Y") {
                                    i_cond->SetValue(DISPLACEMENT_Y, value);
                                } else if (component_variable_name == "DISPLACEMENT_Z") {
                                    i_cond->SetValue(DISPLACEMENT_Z, value);
                                } else if (component_variable_name == "FORCE_X") {
                                    i_cond->SetValue(FORCE_X, value);
                                } else if (component_variable_name == "FORCE_Y") {
                                    i_cond->SetValue(FORCE_Y, value);
                                } else if (component_variable_name == "FORCE_Z") {
                                    i_cond->SetValue(FORCE_Z, value);
                                } else {
                                    KRATOS_ERROR << "No name found" ;
                                }
                            } 
                        } else {KRATOS_ERROR << "AssignIgaExternalConditionsProcess : No Condition or Elements defined" ;}
                    }
                }
            }
        }
    }
}
} // End namespace Kratos
