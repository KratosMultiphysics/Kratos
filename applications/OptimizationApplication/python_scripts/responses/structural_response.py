# importing the Kratos Library
from numpy import gradient
import KratosMultiphysics as KM
from KratosMultiphysics import Parameters, Logger
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

import time as timer
import numpy as np

# ==============================================================================
class StrainEnergyResponseFunction(ResponseFunctionInterface):
    """Linear strain energy response function. It triggers the primal analysis and
    uses the primal analysis results to calculate response value and gradient.

    Attributes
    ----------
    primal_model_part : Model part of the primal analysis object
    primal_analysis : Primal analysis object of the response function
    response_function_utility: Cpp utilities object doing the actual computation of response value and gradient.
    """

    def __init__(self,response_name, response_settings,response_analysis,model):

        self.response_settings = response_settings
        default_gradient_settings = KM.Parameters("""
        {
            "gradient_mode" : "semi_analytic",
            "step_size" : 1e-6
        }""")
        
        self.response_settings["gradient_settings"].ValidateAndAssignDefaults(default_gradient_settings)        

        self.supported_design_types = ["shape"]
        self.name = response_name
        self.primal_analysis = response_analysis
        self.model = model
        self.primal_model_part = self.primal_analysis._GetSolver().GetComputingModelPart()

        design_model_parts = self.response_settings["design_model_parts"].GetStringArray()

        for design_model_part in design_model_parts:
            design_model_part_splitted = design_model_part.split(".")
            if not design_model_part_splitted[0] == self.primal_model_part.Name:
                raise RuntimeError("StrainEnergyResponseFunction:init: root design_model_part {} of response '{}' is not the analysis model!".format(design_model_part_splitted[0],self.name))
                
            self.primal_model_part.AddNodalSolutionStepVariable(KM.SHAPE_SENSITIVITY)

        self.response_function_utility = StructuralMechanicsApplication.StrainEnergyResponseFunctionUtility(self.primal_model_part, self.response_settings["gradient_settings"])

    def Initialize(self):

        self.evaluate_model_parts = self.response_settings["evaluate_model_parts"].GetStringArray()
        self.design_model_parts = self.response_settings["design_model_parts"].GetStringArray()
        self.design_types = self.response_settings["design_types"].GetStringArray()

        if not len(self.design_model_parts)>0 :
            raise RuntimeError("StrainEnergyResponseFunction:Initialize: 'design_model_parts' of response '{}' can not be empty !".format(self.name))

        if not len(self.design_types)>0 :
            raise RuntimeError("StrainEnergyResponseFunction:Initialize: 'design_types' of response '{}' can not be empty !".format(self.name))            

        if not len(self.evaluate_model_parts)>0:
            raise RuntimeError("StrainEnergyResponseFunction:Initialize: 'evaluate_model_parts' of response '{}' can not be empty !".format(self.name))

        if len(self.design_types) != len(self.design_model_parts):
            raise RuntimeError("StrainEnergyResponseFunction:Initialize: 'design_types' and 'design_model_parts' of response '{}' should be of the same size !".format(self.name))

        for evaluate_model_part in self.evaluate_model_parts:
            evaluate_model_part_splitted = evaluate_model_part.split(".")
            if not evaluate_model_part_splitted[0] == self.primal_model_part.Name:
                raise RuntimeError("StrainEnergyResponseFunction:Initialize: root evaluate_model_part {} of response '{}' is not the analysis model!".format(evaluate_model_part_splitted[0],self.name))
            if not self.model.HasModelPart(evaluate_model_part): 
                raise RuntimeError("StrainEnergyResponseFunction:Initialize: evaluate_model_part {} of response '{}' does not exist!".format(evaluate_model_part,self.name))

        for design_type in self.design_types:
            if not design_type in self.supported_design_types:
                raise RuntimeError("StrainEnergyResponseFunction: design type {} of response '{}' is not supported !".format(design_type,self.name))

        self.response_function_utility.Initialize()

    def CalculateValue(self):
        Logger.PrintInfo("StrainEnergyResponse:CalculateValue: Starting value calculation for response ", self.name)
        startTime = timer.time()
        self.value = self.response_function_utility.CalculateValue()
        Logger.PrintInfo("StrainEnergyResponse:CalculateValue: Time needed for calculating value ",round(timer.time() - startTime,2),"s")        
        return self.value

    def GetValue(self):
        self.value = self.response_function_utility.CalculateValue()
        return self.value        

    def CalculateGradients(self):
        Logger.PrintInfo("StrainEnergyResponse", "Starting gradient calculation for response ", self.name)
        startTime = timer.time()
        self.response_function_utility.CalculateGradient()
        Logger.PrintInfo("StrainEnergyResponse", "Time needed for calculating gradients ",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self,design_type_model_part_dict):

        if type(design_type_model_part_dict) is not dict or not bool(design_type_model_part_dict):
            raise RuntimeError("StrainEnergyResponseFunction:CalculateGradient: the input entry should be a dict of a pair of design type and model part ")

        design_type = design_type_model_part_dict.keys()[0]
        design_model_part_name = design_type_model_part_dict.values()[0]

        if self.design_types[design_type] != design_model_part_name or not design_type in self.design_types or not design_model_part_name in self.design_model_parts :
            raise RuntimeError("StrainEnergyResponseFunction:CalculateGradient: requested gradient pair {} does not match with {} of response {}".format(design_type_model_part_dict,dict(design_type,self.design_types_model_part_dict[design_type])),self.name)
        
        Logger.PrintInfo("StrainEnergyResponse", "Starting gradient calculation for response ", self.name, " over model part ",design_model_part_name)
        startTime = timer.time()
        self.response_function_utility.CalculateGradient()
        Logger.PrintInfo("StrainEnergyResponse", "Time needed for calculating gradients ",round(timer.time() - startTime,2),"s")

    def GetGradient(self, design_type_model_part_dict):

        if type(design_type_model_part_dict) is not dict or not bool(design_type_model_part_dict):
            raise RuntimeError("StrainEnergyResponseFunction:GetGradient: the input entry should be a dict of a pair of design type and model part ")

        design_type = design_type_model_part_dict.keys()[0]
        design_model_part_name = design_type_model_part_dict.values()[0]

        if self.design_types[design_type] != design_model_part_name or not design_type in self.design_types or not design_model_part_name in self.design_model_parts :
            raise RuntimeError("StrainEnergyResponseFunction:GetGradient: requested gradient pair {} does not match with {} of response {}".format(design_type_model_part_dict,dict(design_type,self.design_types_model_part_dict[design_type])),self.name)

        
        model_part = self.model.GetModelPart(design_model_part_name)
        gradients = np.zeros(3*model_part.NumberOfNodes())
        index = 0
        for node in model_part.Nodes:
            gradients[index] = node.GetSolutionStepValue(KM.SHAPE_SENSITIVITY_X)
            gradients[index+1] = node.GetSolutionStepValue(KM.SHAPE_SENSITIVITY_Y)
            gradients[index+2] = node.GetSolutionStepValue(KM.SHAPE_SENSITIVITY_Z)
            index +=3
        
        return gradients

    def GetGradients(self):
      
        tot_num_nodes = 0
        for design_model_part_name in self.design_model_parts:
            model_part = self.model.GetModelPart(design_model_part_name)
            tot_num_nodes += model_part.NumberOfNodes()        
        all_gradients = np.zeros(3*tot_num_nodes)

        index = 0
        for design_model_part_name in self.design_model_parts:
            model_part = self.model.GetModelPart(design_model_part_name)           
            for node in model_part.Nodes:
                all_gradients[index] = node.GetSolutionStepValue(KM.SHAPE_SENSITIVITY_X)
                all_gradients[index+1] = node.GetSolutionStepValue(KM.SHAPE_SENSITIVITY_Y)
                all_gradients[index+2] = node.GetSolutionStepValue(KM.SHAPE_SENSITIVITY_Z)
                index +=3

        return all_gradients