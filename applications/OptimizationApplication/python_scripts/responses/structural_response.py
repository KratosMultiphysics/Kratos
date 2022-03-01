# importing the Kratos Library
from numpy import gradient
import KratosMultiphysics as KM
from KratosMultiphysics import Parameters, Logger
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
import KratosMultiphysics.OptimizationApplication as KOA
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

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

        if not self.response_settings.Has("gradient_settings"):
            self.gradient_settings = KM.Parameters()
            self.gradient_settings.AddString("gradient_mode","semi_analytic")
            self.gradient_settings.AddDouble("step_size",1e-6)
        else:
            self.gradient_settings = self.response_settings["gradient_settings"]     

        self.supported_control_types = ["shape"]
        self.model = model
        self.name = response_name
        self.primal_analysis = response_analysis
        self.primal_model_part = self.primal_analysis._GetSolver().GetComputingModelPart()

        self.evaluated_model_parts = response_settings["evaluated_objects"].GetStringArray()
        self.controlled_model_parts = response_settings["controlled_objects"].GetStringArray()
        self.control_types = response_settings["control_types"].GetStringArray()  

        if len(self.evaluated_model_parts) != 1:
            raise RuntimeError("StrainEnergyResponseFunction: 'evaluated_objects' of response '{}' must have only one entry !".format(self.name)) 

        for control_type in self.control_types:
            if not control_type in self.supported_control_types:
                raise RuntimeError("StrainEnergyResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types)) 
        
        # add vars
        for control_type in self.control_types:
            if control_type == "shape":
                self.primal_model_part.AddNodalSolutionStepVariable(KM.SHAPE_SENSITIVITY)
                self.primal_model_part.AddNodalSolutionStepVariable(KOA.D_STRAIN_ENERGY_D_X)
                self.primal_model_part.AddNodalSolutionStepVariable(KOA.D_STRAIN_ENERGY_D_CX)

        # create response
        self.response_function_utility = StructuralMechanicsApplication.StrainEnergyResponseFunctionUtility(self.primal_model_part, self.gradient_settings)

    def Initialize(self):

        for evaluated_model_part in self.evaluated_model_parts:
            evaluated_model_part_splitted = evaluated_model_part.split(".")
            if not evaluated_model_part_splitted[0] == self.primal_model_part.Name:
                raise RuntimeError("StrainEnergyResponseFunction:Initialize: root evaluated_model_part {} of response '{}' is not the analysis model!".format(evaluated_model_part_splitted[0],self.name))
            if not self.model.HasModelPart(evaluated_model_part): 
                raise RuntimeError("StrainEnergyResponseFunction:Initialize: evaluated_model_part {} of response '{}' does not exist!".format(evaluated_model_part,self.name))

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

    def CalculateGradientsForTypeAndObjects(self,control_type,controlled_objects,raise_error=True):

        if raise_error:
            if not control_type in self.control_types:
                raise RuntimeError("StrainEnergyResponseFunction:CalculateGradientsForTypeAndObjects: control type ",control_type," is not supported ")
            if not set(controlled_objects) <=set(self.controlled_model_parts):
                raise RuntimeError("StrainEnergyResponseFunction:CalculateGradientsForTypeAndObjects: controlled_objects ",controlled_objects," do not belong to response ",self.name)

        Logger.PrintInfo("StrainEnergyResponse", "Starting ",control_type," gradient calculation of response ", self.name," for ",controlled_objects)
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