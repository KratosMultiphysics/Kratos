# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#
# ===============================================================================

import KratosMultiphysics as KM
from KratosMultiphysics import Parameters, Logger
from KratosMultiphysics.vtk_output_process import VtkOutputProcess

# Additional imports
import csv
# ==============================================================================
class OptimizationAlgorithm:
    def __init__(self,name,opt_settings,model,model_parts_controller,analyses_controller,responses_controller,controls_controller):

        self.name = name
        self.opt_settings =  opt_settings
        self.model=model
        self.model_parts_controller = model_parts_controller
        self.analyses_controller = analyses_controller
        self.responses_controller = responses_controller
        self.controls_controller = controls_controller

        # objectives
        self.objectives = self.opt_settings["objectives"].GetStringArray()
        self.objectives_improvements = self.opt_settings["objectives_improvements"].GetVector()

        # constraints
        self.constraints = self.opt_settings["constraints"].GetStringArray()
        self.constraints_types = self.opt_settings["constraints_types"].GetStringArray()
        self.constraints_ref_values = self.opt_settings["constraints_ref_values"].GetVector()
        # all responses
        self.responses = self.objectives + self.constraints

        # now extract analyses belong to responses
        self.analyses = self.responses_controller.GetResponsesAnalyses(self.responses)        

        # controls
        self.controls = opt_settings["controls"].GetStringArray()
        self.objectives_controls_weights = opt_settings["objectives_controls_weights"].GetVector()
        self.constraints_controls_weights = opt_settings["constraints_controls_weights"].GetVector()

        # algorithm_settings
        self.max_iterations = self.opt_settings["algorithm_settings"]["max_iterations"].GetInt()
        self.projection_step_size = self.opt_settings["algorithm_settings"]["projection_step_size"].GetDouble()
        self.correction_step_size = self.opt_settings["algorithm_settings"]["correction_step_size"].GetDouble()

        # compile settings for responses
        self.responses_controls = {}
        self.responses_types = {}
        self.responses_controlled_objects = {}
        self.responses_control_gradient_names = {}
        self.responses_control_types = {}
        self.responses_control_var_names = {}
        self.responses_var_names = {}
        self.controls_responses={}
        self.controls_response_var_names = {}
        self.controls_response_gradient_names = {}
        self.controls_response_control_gradient_names = {}
        self.root_model_part_data_field_names = {}
        self.analyses_responses = {}
        self.analysis_free_responses = []
        for response in self.responses:
            response_type = self.responses_controller.GetResponseType(response)
            response_analysis = self.responses_controller.GetResponseAnalysis(response)
            if response_analysis != None:
                if response_analysis in self.analyses_responses.keys():
                    self.analyses_responses[response_analysis].append(response)
                else:
                    self.analyses_responses[response_analysis] = [response]
            else:
                self.analysis_free_responses.append(response)

            self.responses_types[response] = response_type
            response_controlled_objects = self.responses_controller.GetResponseControlledObjects(response)
            response_control_types = self.responses_controller.GetResponseControlTypes(response)
            response_variable_name = self.responses_controller.GetResponseVariableName(response)
            self.responses_var_names[response] = response_variable_name
            for control in self.controls:
                control_type = self.controls_controller.GetControlType(control)
                control_variable_name = self.controls_controller.GetControlVariableName(control)
                control_update_name = self.controls_controller.GetControlUpdateName(control)
                control_output_names = self.controls_controller.GetControlOutputNames(control)
                control_controlling_objects = self.controls_controller.GetControlControllingObjects(control)
                for control_controlling_object in control_controlling_objects:
                    response_controlled_object_index = 0
                    for response_controlled_object in response_controlled_objects:
                        response_controlled_object_type = response_control_types[response_controlled_object_index]
                        if control_controlling_object == response_controlled_object and response_controlled_object_type == control_type:
                            response_gradient_name = self.responses_controller.GetResponseGradientVariableNameForType(response,control_type)
                            response_control_gradient_field = "D_"+response_variable_name+"_D_"+control_variable_name
                            if response in self.responses_controlled_objects.keys():
                                self.responses_controlled_objects[response].append(control_controlling_object)
                                self.responses_controls[response].append(control)
                                self.responses_control_types[response].append(control_type)
                                self.responses_control_var_names[response].append(control_variable_name)
                                self.responses_control_gradient_names[response].append(response_control_gradient_field)
                            else:
                                self.responses_controlled_objects[response]=[control_controlling_object]
                                self.responses_control_types[response]=[control_type]
                                self.responses_control_var_names[response]=[control_variable_name]
                                self.responses_control_gradient_names[response]=[response_control_gradient_field]
                                self.responses_controls[response]=[control]   

                            
                            control_controlling_root_model_part = self.model_parts_controller.GetRootModelPart(control_controlling_object)
                            extracted_root_model_part_name = control_controlling_object.split(".")[0]
                            
                            if extracted_root_model_part_name in self.root_model_part_data_field_names.keys():
                                if not response_control_gradient_field in self.root_model_part_data_field_names[extracted_root_model_part_name]:
                                    self.root_model_part_data_field_names[extracted_root_model_part_name].append(response_control_gradient_field)
                                    self.root_model_part_data_field_names[extracted_root_model_part_name].append(response_gradient_name)
                                    self.root_model_part_data_field_names[extracted_root_model_part_name].extend(control_output_names)
                                    control_controlling_root_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(response_control_gradient_field))
                            else:
                                self.root_model_part_data_field_names[extracted_root_model_part_name] = [response_control_gradient_field,response_gradient_name]
                                self.root_model_part_data_field_names[extracted_root_model_part_name].extend(control_output_names)
                                control_controlling_root_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(response_control_gradient_field))

                            if control in self.controls_responses.keys():
                                if not response in self.controls_responses[control]:
                                    self.controls_responses[control].append(response)
                                    self.controls_response_var_names[control].append(response_variable_name)
                                    self.controls_response_gradient_names[control].append(response_gradient_name)
                                    self.controls_response_control_gradient_names[control].append(response_control_gradient_field)                                   
                            else:
                                self.controls_responses[control] = [response]
                                self.controls_response_var_names[control] = [response_variable_name]
                                self.controls_response_gradient_names[control]= [response_gradient_name]
                                self.controls_response_control_gradient_names[control] = [response_control_gradient_field]  
                        
                        response_controlled_object_index +=1                          
        
        # compile settings for c++ optimizer
        self.opt_parameters = Parameters()
        self.opt_parameters.AddEmptyArray("objectives")
        self.opt_parameters.AddEmptyArray("constraints")
        self.opt_parameters.AddInt("opt_itr",0)
        self.opt_parameters.AddInt("max_opt_itr",self.max_iterations)
        self.opt_parameters.AddDouble("projection_step_size",self.projection_step_size)
        self.opt_parameters.AddDouble("correction_step_size",self.correction_step_size)
        self.opt_parameters.AddInt("num_active_consts",0)
        self.opt_parameters.AddDouble("sin_alpha",0)
        self.opt_parameters.AddBool("opt_converged",False)
        for response in self.responses_controlled_objects.keys():
            response_settings = Parameters()
            response_settings.AddString("name",response)
            response_settings.AddString("response_type",self.responses_types[response])
            response_settings.AddString("variable_name",self.responses_var_names[response])
            response_settings.AddDouble("value",0.0)
            response_settings.AddDouble("weight",1.0)
            response_settings.AddDouble("init_value",0.0)
            response_settings.AddDouble("prev_itr_value",0.0)
            response_settings.AddEmptyArray("controlled_objects")
            response_settings.AddEmptyArray("control_types")
            response_settings.AddEmptyArray("control_gradient_names")
            response_settings.AddEmptyArray("control_variable_names")
            response_settings.AddEmptyArray("controls")
                
            for control_obj in self.responses_controlled_objects[response]:
                response_settings["controlled_objects"].Append(control_obj)

            for control_type in self.responses_control_types[response]:
                response_settings["control_types"].Append(control_type)

            for control_variable in self.responses_control_var_names[response]:
                response_settings["control_variable_names"].Append(control_variable)

            for control_gradient_name in self.responses_control_gradient_names[response]:
                response_settings["control_gradient_names"].Append(control_gradient_name)

            for control in self.responses_controls[response]:
                response_settings["controls"].Append(control)

            if response in self.objectives:
                index = self.objectives.index(response)
                response_settings.AddDouble("objective_improvement",self.objectives_improvements[index])
                self.opt_parameters["objectives"].Append(response_settings)
            elif response in self.constraints:
                index = self.constraints.index(response)
                response_settings.AddString("type",self.constraints_types[index])
                response_settings.AddDouble("ref_value",self.constraints_ref_values[index])
                response_settings.AddBool("is_active",True)
                response_settings.AddBool("prev_itr_is_active",False)
                self.opt_parameters["constraints"].Append(response_settings)
            else:
                raise RuntimeError("OptimizationAlgorithm:__init__:error in compile settings for c++ optimizer")

        self.opt_parameters.AddEmptyArray("controls")
        for control in self.controls:
            control_index = self.controls.index(control)
            control_objectives_weight = self.objectives_controls_weights[control_index]
            control_constraints_weight = self.constraints_controls_weights[control_index]
            control_type = self.controls_controller.GetControlType(control)
            control_update_name = self.controls_controller.GetControlUpdateName(control)
            control_variable_name = self.controls_controller.GetControlVariableName(control)
            control_settings = Parameters()
            control_settings.AddString("name",control)
            control_settings.AddString("type",control_type)
            control_settings.AddString("update_name",control_update_name)
            control_settings.AddString("variable_name",control_variable_name)
            control_settings.AddDouble("objectives_weight",control_objectives_weight)
            control_settings.AddDouble("constraints_weight",control_constraints_weight)
            if control_type == "shape":
                control_settings.AddInt("size",3)
            elif control_type == "thickness":
                control_settings.AddInt("size",1)
            elif control_type == "material":
                control_settings.AddInt("size",1)
            else:
                raise RuntimeError("OptimizationAlgorithm:__init__:error in compile settings for c++ optimizer")

            control_settings.AddEmptyArray("controlling_objects")

            control_controlling_objects = self.controls_controller.GetControlControllingObjects(control)
            for control_controlling_object in control_controlling_objects:
                control_settings["controlling_objects"].Append(control_controlling_object) 

            self.opt_parameters["controls"].Append(control_settings) 

        Logger.PrintInfo("::[OptimizationAlgorithm]:: ", "Variables ADDED")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop( self ):

        # create vtkIOs
        self.root_model_parts_vtkIOs = {}
        for root_model_part,data_fields in self.root_model_part_data_field_names.items():
            vtk_parameters = Parameters()
            root_controlling_model_part = self.model_parts_controller.GetRootModelPart(root_model_part)
            vtk_parameters.AddString("model_part_name",root_controlling_model_part.Name)
            vtk_parameters.AddBool("write_ids",False)
            vtk_parameters.AddString("file_format","ascii")
            vtk_parameters.AddBool("output_sub_model_parts",False)
            vtk_parameters.AddString("output_path","Optimization_Results")
            vtk_parameters.AddEmptyArray("nodal_solution_step_data_variables")
            for nodal_result in data_fields:
                vtk_parameters["nodal_solution_step_data_variables"].Append(nodal_result)
            
            controlling_model_part_vtkIO = VtkOutputProcess(self.model, vtk_parameters)
            controlling_model_part_vtkIO.ExecuteInitialize()
            controlling_model_part_vtkIO.ExecuteBeforeSolutionLoop()
            self.root_model_parts_vtkIOs[root_model_part] = controlling_model_part_vtkIO

    # --------------------------------------------------------------------------
    def RunOptimizationLoop( self ):
        raise RuntimeError("Algorithm base class is called. Please check your implementation of the function >> RunOptimizationLoop << .")

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop( self ):
        for vtkIO in self.root_model_parts_vtkIOs.values():
            vtkIO.ExecuteFinalize() 

    # --------------------------------------------------------------------------
    def SetResponseValue( self,response,value ):

        for objective in self.opt_parameters["objectives"]:
            if objective["name"].GetString() == response:
                objective["prev_itr_value"].SetDouble(objective["value"].GetDouble())
                objective["value"].SetDouble(value)
                if self.opt_parameters["opt_itr"].GetInt()<2:
                   objective["init_value"].SetDouble(value)   
                
                return  True

        for constraint in self.opt_parameters["constraints"]:
            if constraint["name"].GetString() == response:
                constraint["prev_itr_value"].SetDouble(constraint["value"].GetDouble())
                constraint["value"].SetDouble(value) 
                if self.opt_parameters["opt_itr"].GetInt()<2:
                    constraint["init_value"].SetDouble(value)
                    valid_types = ["initial_value_equality","smaller_than_initial_value","bigger_than_initial_value"]
                    if constraint["type"].GetString() in valid_types:
                        constraint["ref_value"].SetDouble(value)

                ref_value = constraint["ref_value"].GetDouble()

                type = constraint["type"].GetString()
                is_active = False
                if type == "equality" or type == "initial_value_equality":
                    is_active = True
                elif (type == "smaller_than" or type == "smaller_than_initial_value") and value>ref_value:
                    is_active = True
                elif (type == "bigger_than" or type == "bigger_than_initial_value")  and value<ref_value:
                    is_active = True

                if is_active:
                    self.num_active_consts += 1

                constraint["prev_itr_is_active"].SetBool(constraint["is_active"].GetBool())
                constraint["is_active"].SetBool(is_active)

                return is_active       

    # --------------------------------------------------------------------------
    def _InitializeCSVLogger(self):
        self.complete_log_file_name = "optimization_log.csv"
        with open(self.complete_log_file_name, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>4s}".format("itr"))
            for objective in self.opt_parameters["objectives"]:
                row.append("{:>4s}".format("f: "+str(objective["name"].GetString())))
                row.append("{:>4s}".format("abs[%]"))
                row.append("{:>4s}".format("rel[%]"))
                row.append("{:>4s}".format("weight"))

            for constraint in self.opt_parameters["constraints"]:
                row.append("{:>4s}".format("g: "+str(constraint["name"].GetString())))
                row.append("{:>4s}".format("ref_val "))
                row.append("{:>4s}".format("ref_diff[%]"))
                row.append("{:>4s}".format("weight"))
                row.append("{:>4s}".format("is_active"))


            row.append("{:>4s}".format("projection_step_size"))
            if self.opt_parameters["constraints"].size() > 0:
                row.append("{:>4s}".format("correction_step_size"))
                row.append("{:>4s}".format("sin_alpha"))

            historyWriter.writerow(row)

    # --------------------------------------------------------------------------
    def _WriteCurrentOptItrToCSVFile( self ):
        
        with open(self.complete_log_file_name, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>4d}".format(self.optimization_iteration))

            Logger.Print("")

            for objective in self.opt_parameters["objectives"]:
                objectivs_current_value = objective["value"].GetDouble()
                objectivs_weight = objective["weight"].GetDouble()
                objective_initial_value = objectivs_current_value
                objective_previous_value = objectivs_current_value

                if self.optimization_iteration>1:                
                    objective_initial_value = objective["init_value"].GetDouble() 
                    objective_previous_value = objective["prev_itr_value"].GetDouble() 

                rel_change = 100 * (objectivs_current_value-objective_previous_value)/objective_previous_value
                abs_change = 100 * (objectivs_current_value-objective_initial_value)/objective_initial_value

                Logger.Print("  ===== objective: ",objective["name"].GetString())
                Logger.Print("                   current value: ",objectivs_current_value)
                Logger.Print("                   rel_change: ",rel_change)
                Logger.Print("                   abs_change: ",abs_change)
                Logger.Print("                   weight: ",objectivs_weight)
                row.append(" {:> .5E}".format(objectivs_current_value))
                row.append(" {:> .5E}".format(abs_change))
                row.append(" {:> .5E}".format(rel_change)) 
                row.append(" {:> .5E}".format(objectivs_weight))

            for constraint in self.opt_parameters["constraints"]:
                constraint_current_value = constraint["value"].GetDouble() 
                constraint_ref_val = constraint["ref_value"].GetDouble()
                constraint_weight = constraint["weight"].GetDouble()
                is_active = constraint["is_active"].GetBool()
                abs_change = 100 * abs(constraint_current_value-constraint_ref_val)/abs(constraint_ref_val)

                Logger.Print("  ===== constraint: ",constraint["name"].GetString())
                Logger.Print("                   current value: ",constraint_current_value)
                Logger.Print("                   ref value: ",constraint_ref_val)
                Logger.Print("                   change: ",abs_change)
                Logger.Print("                   weight: ",constraint_weight)
                Logger.Print("                   is_active: ",is_active)

                row.append(" {:> .5E}".format(constraint_current_value))
                row.append(" {:> .5E}".format(constraint_ref_val))
                row.append(" {:> .5E}".format(abs_change))
                row.append(" {:> .5E}".format(constraint_weight))
                row.append(" {:> .5E}".format(is_active))
                

            Logger.Print("  ===== projection_step_size: ",self.opt_parameters["projection_step_size"].GetDouble())            
            row.append(" {:> .5E}".format(self.opt_parameters["projection_step_size"].GetDouble()))

            if self.opt_parameters["constraints"].size() > 0 and self.opt_parameters["num_active_consts"].GetInt()>0:
                Logger.Print("  ===== correction_step_size: ",self.opt_parameters["correction_step_size"].GetDouble())
                row.append(" {:> .5E}".format(self.opt_parameters["correction_step_size"].GetDouble()))
                Logger.Print("  ===== sin alpha: ",self.opt_parameters["sin_alpha"].GetDouble())
                row.append(" {:> .5E}".format(self.opt_parameters["sin_alpha"].GetDouble()))

            historyWriter.writerow(row)

# ==============================================================================
