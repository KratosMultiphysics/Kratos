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
        self.supported_control_types_variables_name = self.controls_controller.GetSupportedControlTypesVariablesName()
        self.controls_type = {}
        self.controls_controling_objects = {}
        self.controls_controling_objects = {}
        self.controls_responses = {}
        self.controls_responses_model_parts = {}
        for control in self.controls:
            control_type = self.controls_controller.GetControlType(control)
            self.controls_type[control] = control_type
            control_variable_name =  self.supported_control_types_variables_name[control_type]
            controls_controlling_parts = self.controls_controller.GetControlControllingObjects(control)
            responses_dict = self.responses_controller.GetResponsesForControl(control_type,controls_controlling_parts)
            # now we need to do the checks
            if not len(responses_dict)>0:
                raise RuntimeError("OptimizationAlgorithm: could not associate control {} to any response !.".format(control))
            # remove responses that are not in the response list
            all_found_responses = list(responses_dict.keys())
            for response in all_found_responses:
                if response not in self.responses:
                    del responses_dict[response]

            all_found_objects = []
            for objects in responses_dict.values():
                all_found_objects.extend(objects)  
            all_found_objects = list(set(all_found_objects))
            for object in controls_controlling_parts:
                if not object in all_found_objects:
                    raise RuntimeError("OptimizationAlgorithm: could not associate control object {} of control {} to any response !".format(object,control))
       
            self.controls_responses_model_parts[control] = responses_dict
            # add data fields here
            for response,controlled_objects in responses_dict.items():
                control_controlling_root_model_parts = self.model_parts_controller.GetRootModelParts(controlled_objects)
                response_variable_name = self.responses_controller.GetResponseVariableName(response)
                response_control_gradient_field = "D_"+response_variable_name+"_D_"+control_variable_name
                for root_model in control_controlling_root_model_parts:
                    root_model.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(response_control_gradient_field))

        # check that we could associate all responses to controls
        all_found_responses = []
        for control_associated_responses in self.controls_responses_model_parts.values():
            all_found_responses.extend(list(control_associated_responses.keys()))

        if not set(all_found_responses) == set(self.responses):    
            raise RuntimeError("OptimizationAlgorithm: could not associate controls to any responses !")

        Logger.PrintInfo("::[OptimizationAlgorithm]:: ", "Variables ADDED")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop( self ):
        raise RuntimeError("Algorithm base class is called. Please check your implementation of the function >> InitializeOptimizationLoop << .")

    # --------------------------------------------------------------------------
    def RunOptimizationLoop( self ):
        raise RuntimeError("Algorithm base class is called. Please check your implementation of the function >> RunOptimizationLoop << .")

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop( self ):
        raise RuntimeError("Algorithm base class is called. Please check your implementation of the function >> FinalizeOptimizationLoop << .")

    def __CreateVTKIOs(self):
        self.controls_IOs={}
        for control,controlling_objects in self.controls_controlling_objects:
            vtk_parameters = Parameters()        

# ==============================================================================
