# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#
# ===============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication as KOA
from KratosMultiphysics import Parameters, Logger
from KratosMultiphysics.vtk_output_process import VtkOutputProcess

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
        self.objectives_weights = self.opt_settings["objectives_weights"].GetVector()


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
        self.controls_maximum_updates = opt_settings["controls_maximum_updates"].GetVector()

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
        for response in self.responses_controlled_objects.keys():
            response_settings = Parameters()
            response_settings.AddString("name",response)
            response_settings.AddString("response_type",self.responses_types[response])
            response_settings.AddString("variable_name",self.responses_var_names[response])
            response_settings.AddDouble("value",0.0)
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
                response_settings.AddDouble("objective_weight",self.objectives_weights[index])
                self.opt_parameters["objectives"].Append(response_settings)
            elif response in self.constraints:
                index = self.constraints.index(response)
                response_settings.AddString("type",self.constraints_types[index])
                response_settings.AddDouble("ref_value",self.constraints_ref_values[index])
                self.opt_parameters["constraints"].Append(response_settings)
            else:
                raise RuntimeError("OptimizationAlgorithm:__init__:error in compile settings for c++ optimizer")

        self.opt_parameters.AddEmptyArray("controls")
        for control in self.controls:
            control_index = self.controls.index(control)
            control_max_update = self.controls_maximum_updates[control_index]
            control_type = self.controls_controller.GetControlType(control)
            control_update_name = self.controls_controller.GetControlUpdateName(control)
            control_variable_name = self.controls_controller.GetControlVariableName(control)
            control_settings = Parameters()
            control_settings.AddString("name",control)
            control_settings.AddString("type",control_type)
            control_settings.AddString("update_name",control_update_name)
            control_settings.AddString("variable_name",control_variable_name)
            control_settings.AddDouble("max_update",control_max_update)
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

        self.constraints_hist = {}
        self.objectives_hist = {}

        for objective in self.objectives:
            self.objectives_hist[objective] = []

        for constraint in self.constraints:
            self.constraints_hist[constraint] = []

    # --------------------------------------------------------------------------
    def RunOptimizationLoop( self ):
        raise RuntimeError("Algorithm base class is called. Please check your implementation of the function >> RunOptimizationLoop << .")

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop( self ):
        for vtkIO in self.root_model_parts_vtkIOs.values():
            vtkIO.ExecuteFinalize() 

    # --------------------------------------------------------------------------
    def SetResponseValue( self,response,value ):

        for objective in self.objectives:
            if objective == response:
                self.objectives_hist[objective].append(value)
                break

        for constraint in self.constraints:
            if constraint == response:
                self.constraints_hist[constraint].append(value)
                break

        for objective in self.opt_parameters["objectives"]:
            if objective["name"].GetString() == response:
                objective["value"].SetDouble(value)  

        for constraint in self.opt_parameters["constraints"]:
            if constraint["name"].GetString() == response:
                constraint["value"].SetDouble(value) 

    # --------------------------------------------------------------------------
    def _InitializeCSVLogger(self):
        self.complete_log_file_name = "optimization_log.csv"
        with open(self.complete_log_file_name, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>4s}".format("itr"))
            for itr in range(len(self.objectives)):
                row.append("{:>4s}".format("f: "+str(self.objectives[itr])))
                row.append("{:>4s}".format("abs[%]"))
                row.append("{:>4s}".format("rel[%]"))

            for itr in range(len(self.constraints)):
                row.append("{:>4s}".format("c: "+str(self.constraints[itr])))
                row.append("{:>4s}".format("ref_val "))
                row.append("{:>4s}".format("ref_diff[%]"))

            historyWriter.writerow(row)

    # --------------------------------------------------------------------------
    def _WriteCurrentResponseValuesToCSVFile( self ):
        
        with open(self.complete_log_file_name, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>4d}".format(self.optimization_iteration))
            for objective,objective_hist_values in self.objectives_hist.items():
                objectivs_current_value = objective_hist_values[self.optimization_iteration-1]
                objective_initial_value = objectivs_current_value
                objective_previous_value = objectivs_current_value
                if self.optimization_iteration>1:                
                    objective_initial_value = objective_hist_values[0]
                    objective_previous_value = objective_hist_values[self.optimization_iteration-2]
                rel_change = 100 * (objectivs_current_value-objective_previous_value)/objective_previous_value
                abs_change = 100 * (objectivs_current_value-objective_initial_value)/objective_initial_value
                Logger.Print("  ===== objective: ",objective)
                Logger.Print("                   current value: ",objectivs_current_value)
                Logger.Print("                   rel_change: ",rel_change)
                Logger.Print("                   abs_change: ",abs_change)
                row.append(" {:> .5E}".format(objectivs_current_value))
                row.append(" {:> .5E}".format(abs_change))
                row.append(" {:> .5E}".format(rel_change)) 

            constraint_index = 0
            for constraint,constraint_hist_values in self.constraints_hist.items():
                constraint_current_value = constraint_hist_values[self.optimization_iteration-1]
                constraint_ref_val = self.constraints_ref_values[constraint_index]
                abs_change = 100 * abs(constraint_current_value-constraint_ref_val)/abs(constraint_ref_val)

                Logger.Print("  ===== constraint: ",constraint)
                Logger.Print("                   current value: ",constraint_current_value)
                Logger.Print("                   abs_change: ",abs_change)

                row.append(" {:> .5E}".format(constraint_current_value))
                row.append(" {:> .5E}".format(constraint_ref_val))
                row.append(" {:> .5E}".format(abs_change))
                constraint_index = constraint_index + 1

            historyWriter.writerow(row)

# ==============================================================================
