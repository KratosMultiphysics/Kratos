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

namespace Kratos
{

    AssignIgaExternalConditionsProcess::AssignIgaExternalConditionsProcess(
        Model& rModel, Parameters ThisParameters) : mpModel(&rModel), mParameters(ThisParameters)
    {
        mEchoLevel = mParameters["echo_level"].GetInt();

        KRATOS_ERROR_IF_NOT(mParameters.Has("analysis_model_part_name"))
            << "Missing \"analysis_model_part_name\" section" << std::endl;
        ModelPart& analysis_model_part = mpModel->GetModelPart(mParameters["analysis_model_part_name"].GetString());

        const std::string& rDataFileName = mParameters.Has("physics_file_name")
            ? mParameters["physics_file_name"].GetString()
            : "physics.iga.json";

        iga_physics_parameters = ReadParamatersFile(rDataFileName);
   }

    void AssignIgaExternalConditionsProcess::ExecuteInitializeSolutionStep(){

        ModelPart& analysis_model_part = mpModel->GetModelPart(mParameters["analysis_model_part_name"].GetString());

        const std::vector<std::string> sub_model_part_names = analysis_model_part.GetSubModelPartNames();

        // Number of element_condition_list (Elements or Conditions), at the end are all the Gauss Points
        const int numberOfEntities = iga_physics_parameters["element_condition_list"].size();  

        const ProcessInfo& rProcessInfo = analysis_model_part.GetProcessInfo();
        double t = rProcessInfo[TIME];
        
        

        for (IndexType i = 0; i < numberOfEntities; i++) {
            std::string sub_model_part_name = iga_physics_parameters["element_condition_list"][i]["iga_model_part"].GetString();
            // Get the sub_model_part. Example "SBM_Support_outer"
            ModelPart& sub_model_part = analysis_model_part.GetSubModelPart(sub_model_part_name);

            if (iga_physics_parameters["element_condition_list"][i]["parameters"].Has("variables")) {
                for (IndexType i_var = 0; i_var < iga_physics_parameters["element_condition_list"][i]["parameters"]["variables"].size(); i_var++){
                    // Get the variable. Example "HEAT_FLUX"
                    std::string variable_name = iga_physics_parameters["element_condition_list"][i]["parameters"]["variables"][i_var]["variable_name"].GetString();
                    std::vector<std::string> variable_name_array;
                    std::vector<std::string> values_string;
                    std::vector<std::string> components_array = {"_X","_Y","_Z"};
                    bool isScalarValue = iga_physics_parameters["element_condition_list"][i]["parameters"]["variables"][i_var]["value"].IsString();
                    if (isScalarValue == 1 ) {
                        values_string.push_back(iga_physics_parameters["element_condition_list"][i]["parameters"]["variables"][i_var]["value"].GetString()); 
                        variable_name_array.push_back(variable_name) ;
                    } else {
                        const int numberOfComponents = iga_physics_parameters["element_condition_list"][i]["parameters"]["variables"][i_var]["value"].size();
                        for (IndexType i_component = 0; i_component < numberOfComponents; i_component++ ) {
                            values_string.push_back(iga_physics_parameters["element_condition_list"][i]["parameters"]["variables"][i_var]["value"][i_component].GetString());
                            std::string component_variable_name = variable_name+components_array[i_component];
                            variable_name_array.push_back(component_variable_name) ;
                        }  
                    }  
                    for (IndexType i_component = 0; i_component < values_string.size(); i_component++) {   
                        // Define the function to be parsed      
                        Kratos::unique_ptr<GenericFunctionUtility> evalFunction = Kratos::make_unique<GenericFunctionUtility>(values_string[i_component]);
                        const std::string component_variable_name = variable_name_array[i_component];
                        if (sub_model_part.Elements().size() > 0) {
                            for (auto i_element = sub_model_part.ElementsBegin(); i_element != sub_model_part.ElementsEnd(); i_element++) {
                                const double x = i_element->GetGeometry().Center().X();
                                const double y = i_element->GetGeometry().Center().Y();
                                const double z = i_element->GetGeometry().Center().Z();


                                const double value = evalFunction->CallFunction(x, y, z, t);
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
                                const double value = evalFunction->CallFunction(x, y, z, t);
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
                                } else if (component_variable_name == "PRESSURE") {
                                    i_cond->SetValue(PRESSURE, value);
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


    void AssignIgaExternalConditionsProcess::ExecuteInitialize(){

        ModelPart& analysis_model_part = mpModel->GetModelPart(mParameters["analysis_model_part_name"].GetString());

        const ProcessInfo& rProcessInfo = analysis_model_part.GetProcessInfo();
        double t = rProcessInfo[TIME];
        if (mParameters.Has("initial_variables")) {
            const auto& initial_variables = mParameters["initial_variables"];

            for (IndexType i_var = 0; i_var < initial_variables.size(); i_var++) {
                std::string variable_name = initial_variables[i_var]["variable_name"].GetString();
                bool isScalarValue = initial_variables[i_var]["initial_value"].IsString();

                std::vector<std::string> components_array = {"_X", "_Y", "_Z"};
                std::vector<std::string> variable_name_array;
                std::vector<std::string> values_string;

                if (isScalarValue) {
                    values_string.push_back(initial_variables[i_var]["initial_value"].GetString());
                    variable_name_array.push_back(variable_name);
                } else {
                    const int numberOfComponents = initial_variables[i_var]["initial_value"].size();
                    for (IndexType i_component = 0; i_component < numberOfComponents; i_component++) {
                        values_string.push_back(initial_variables[i_var]["initial_value"][i_component].GetString());
                        std::string component_variable_name = variable_name + components_array[i_component];
                        variable_name_array.push_back(component_variable_name);
                    }
                }
                KRATOS_WATCH(variable_name_array)
                for (auto i_element = analysis_model_part.ElementsBegin(); i_element != analysis_model_part.ElementsEnd(); i_element++) {
                    const double x = i_element->GetGeometry().Center().X();
                    const double y = i_element->GetGeometry().Center().Y();
                    const double z = i_element->GetGeometry().Center().Z();

                    for (IndexType i_component = 0; i_component < values_string.size(); i_component++) {
                        Kratos::unique_ptr<GenericFunctionUtility> evalFunction = Kratos::make_unique<GenericFunctionUtility>(values_string[i_component]);
                        const double value = evalFunction->CallFunction(x, y, z, t);

                        if (variable_name_array[i_component] == "VELOCITY_X") {
                            i_element->SetValue(VELOCITY_X, value);
                        } else if (variable_name_array[i_component] == "VELOCITY_Y") {
                            i_element->SetValue(VELOCITY_Y, value);
                        } else if (variable_name_array[i_component] == "VELOCITY_Z") {
                            i_element->SetValue(VELOCITY_Z, value);
                        } else if (variable_name_array[i_component] == "PRESSURE") {
                            i_element->SetValue(PRESSURE, value);
                        } else if (variable_name_array[i_component] == "FORCE_X") {
                            i_element->SetValue(FORCE_X, value);
                        } else if (variable_name_array[i_component] == "FORCE_Y") {
                            i_element->SetValue(FORCE_Y, value);
                        } else if (variable_name_array[i_component] == "FORCE_Z") {
                            i_element->SetValue(FORCE_Z, value);
                        } else {
                            KRATOS_ERROR << "Unsupported variable: " << variable_name_array[i_component];
                        }
                    }
                }
            }
        }
    }


    /// Reads in a json formatted file and returns its KratosParameters instance.
    Parameters AssignIgaExternalConditionsProcess::ReadParamatersFile(
        const std::string& rDataFileName) const
    {
        // Check if rDataFileName ends with ".cad.json" and add it if needed.
        const std::string data_file_name = (rDataFileName.compare(rDataFileName.size() - 9, 9, ".iga.json") != 0)
            ? rDataFileName + ".iga.json"
            : rDataFileName;

        std::ifstream infile(data_file_name);
        KRATOS_ERROR_IF_NOT(infile.good()) << "Physics fil: "
            << data_file_name << " cannot be found." << std::endl;
        KRATOS_INFO_IF("ReadParamatersFile", mEchoLevel > 3)
            << "Reading file: \"" << data_file_name << "\"" << std::endl;

        std::stringstream buffer;
        buffer << infile.rdbuf();

        return Parameters(buffer.str());
    }
} // End namespace Kratos