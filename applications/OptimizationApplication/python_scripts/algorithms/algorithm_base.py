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

        # all responses
        self.responses = self.objectives + self.constraints

        # now extract analyses belong to responses
        self.analyses = self.responses_controller.GetResponsesAnalyses(self.responses)        

        # controls
        self.controls = opt_settings["controls"].GetStringArray()
        # self.supported_control_types_variables_name = self.controls_controller.GetSupportedControlTypesVariablesName()
        # self.controls_type = {}
        # self.controls_responses_model_parts = {}
        # self.controlling_model_parts = []
        # for control in self.controls:
        #     control_type = self.controls_controller.GetControlType(control)
        #     self.controls_type[control] = control_type
        #     control_variable_name =  self.supported_control_types_variables_name[control_type]
        #     controls_controlling_parts = self.controls_controller.GetControlControllingObjects(control)
        #     self.controlling_model_parts.extend(controls_controlling_parts)
        #     responses_dict = self.responses_controller.GetResponsesForControl(control_type,controls_controlling_parts)
        #     # now we need to do the checks
        #     if not len(responses_dict)>0:
        #         raise RuntimeError("OptimizationAlgorithm: could not associate control {} to any response !.".format(control))
        #     # remove responses that are not in the response list
        #     all_found_responses = list(responses_dict.keys())
        #     for response in all_found_responses:
        #         if response not in self.responses:
        #             del responses_dict[response]

        #     all_found_objects = []
        #     for objects in responses_dict.values():
        #         all_found_objects.extend(objects)  
        #     all_found_objects = list(set(all_found_objects))
        #     for object in controls_controlling_parts:
        #         if not object in all_found_objects:
        #             raise RuntimeError("OptimizationAlgorithm: could not associate control object {} of control {} to any response !".format(object,control))
       
        #     self.controls_responses_model_parts[control] = responses_dict
        #     # add data fields here
        #     for response,controlled_objects in responses_dict.items():
        #         control_controlling_root_model_parts = self.model_parts_controller.GetRootModelParts(controlled_objects)
        #         response_variable_name = self.responses_controller.GetResponseVariableName(response)
        #         response_control_gradient_field = "D_"+response_variable_name+"_D_"+control_variable_name
        #         for root_model in control_controlling_root_model_parts:
        #             root_model.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(response_control_gradient_field))

        # # check that we could associate all responses to controls
        # all_found_responses = []
        # for control_associated_responses in self.controls_responses_model_parts.values():
        #     all_found_responses.extend(list(control_associated_responses.keys()))

        # if not set(all_found_responses) == set(self.responses):    
        #     raise RuntimeError("OptimizationAlgorithm: could not associate controls to any responses !")


        # compile settings for responses
        self.responses_controls = {}
        self.responses_controlled_objects = {}
        self.responses_control_types = {}
        self.responses_control_var_names = {}
        self.responses_var_names = {}
        self.controls_responses={}
        self.controls_response_var_names = {}
        self.controls_response_gradient_names = {}
        self.supported_control_types_variables_name = self.controls_controller.GetSupportedControlTypesVariablesName()
        for response in self.responses:
            response_type = self.responses_controller.GetResponseType(response)
            response_controlled_objects = self.responses_controller.GetResponseControlledObjects(response)
            response_control_types = self.responses_controller.GetResponseControlTypes(response)
            response_variable_name = self.responses_controller.GetResponseVariableName(response)
            self.responses_var_names[response] = response_variable_name
            for control in self.controls:
                control_type = self.controls_controller.GetControlType(control)
                control_variable_name =  self.supported_control_types_variables_name[control_type]
                control_controlling_objects = self.controls_controller.GetControlControllingObjects(control)
                for control_controlling_object in control_controlling_objects:
                    if control_controlling_object in response_controlled_objects:
                        index = response_controlled_objects.index(control_controlling_object)
                        if control_type == response_control_types[index]:
                            response_gradient_name = self.responses_controller.GetResponseGradientVariableNameForType(response,control_type)
                            if response in self.responses_controlled_objects.keys():
                                self.responses_controlled_objects[response].append(control_controlling_object)
                                self.responses_control_types[response].append(control_type)
                                self.responses_control_var_names[response].append(control_variable_name)
                            else:
                                self.responses_controlled_objects[response]=[control_controlling_object]
                                self.responses_control_types[response]=[control_type]
                                self.responses_control_var_names[response]=[control_variable_name]   

                            if response in self.responses_controls.keys():
                                if not control in self.responses_controls[response]:
                                    self.responses_controls[response].append(control)
                            else:
                                self.responses_controls[response]=[control]

                            if control in self.controls_responses.keys():
                                if not response in self.controls_responses[control]:
                                    self.controls_responses[control].append(response)
                                    self.controls_response_var_names[control].append(response_variable_name)
                                    self.controls_response_gradient_names[control].append(response_gradient_name)
                            else:
                                self.controls_responses[control] = [response]
                                self.controls_response_var_names[control] = [response_variable_name]
                                self.controls_response_gradient_names[control]= [response_gradient_name]





        Logger.PrintInfo("::[OptimizationAlgorithm]:: ", "Variables ADDED")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop( self ):

        # create vtkIO for every controlled object
        self.controlling_model_parts = list(set(self.controlling_model_parts)) # remove dependencies
        self.controlling_model_parts_vtkIOs = {}
        for controlling_model_part in self.controlling_model_parts:
            vtk_parameters = Parameters()
            root_controlling_model_part = self.model_parts_controller.GetRootModelPart(controlling_model_part)
            vtk_parameters.AddString("model_part_name",root_controlling_model_part.Name)
            vtk_parameters.AddBool("write_ids",False)
            vtk_parameters.AddString("file_format","ascii")
            vtk_parameters.AddBool("output_sub_model_parts",False)
            vtk_parameters.AddString("output_path","Optimization_Results")
            nodal_results = []
            for control,responses_model_parts_dict in self.controls_responses_model_parts.items():
                control_type = self.controls_type[control]
                control_variable_name =  self.supported_control_types_variables_name[control_type]
                for response,model_parts in responses_model_parts_dict.items():
                    reponse_gradient_variable_name = self.responses_controller.GetResponseGradientVariableNameForType(response,control_type,False)
                    response_variable_name = self.responses_controller.GetResponseVariableName(response)
                    response_control_gradient_variable_name = "D_"+response_variable_name+"_D_"+control_variable_name
                    
                    if controlling_model_part in model_parts:
                        nodal_results.append(reponse_gradient_variable_name)
                        nodal_results.append(response_control_gradient_variable_name)

            nodal_results = list(set(nodal_results))
            vtk_parameters.AddEmptyArray("nodal_solution_step_data_variables")
            for nodal_result in nodal_results:
                vtk_parameters["nodal_solution_step_data_variables"].Append(nodal_result)
            
            controlling_model_part_vtkIO = VtkOutputProcess(self.model, vtk_parameters)
            controlling_model_part_vtkIO.ExecuteInitialize()
            controlling_model_part_vtkIO.ExecuteBeforeSolutionLoop()
            self.controlling_model_parts_vtkIOs[controlling_model_part] = controlling_model_part_vtkIO

    # --------------------------------------------------------------------------
    def RunOptimizationLoop( self ):
        raise RuntimeError("Algorithm base class is called. Please check your implementation of the function >> RunOptimizationLoop << .")

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop( self ):
        for vtkIO in self.controlling_model_parts_vtkIOs.values():
            vtkIO.ExecuteFinalize()     

# ==============================================================================
