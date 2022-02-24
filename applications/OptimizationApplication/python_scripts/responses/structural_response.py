# importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics import Parameters, Logger
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

import time as timer

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
        self.response_function_utility = StructuralMechanicsApplication.StrainEnergyResponseFunctionUtility(self.primal_model_part, self.response_settings["gradient_settings"])

    def Initialize(self):

        self.evaluate_model_parts = self.response_settings["evaluate_model_parts"].GetStringArray()
        self.design_model_parts = self.response_settings["design_model_parts"].GetStringArray()
        self.design_types = self.response_settings["design_types"].GetStringArray()

        if not len(self.evaluate_model_parts)>0:
            raise RuntimeError("StrainEnergyResponseFunction: 'evaluate_model_parts' of response '{}' can not be empty !".format(self.name))

        for evaluate_model_part in self.evaluate_model_parts:
            evaluate_model_part_splitted = evaluate_model_part.split(".")
            if not evaluate_model_part_splitted[0] == self.primal_model_part.Name:
                raise RuntimeError("StrainEnergyResponseFunction: root evaluate_model_part {} of response '{}' is not the analysis model!".format(evaluate_model_part_splitted[0],self.name))
            if not self.model.HasModelPart(evaluate_model_part): 
                raise RuntimeError("StrainEnergyResponseFunction: evaluate_model_part {} of response '{}' does not exist!".format(evaluate_model_part,self.name))

        if not len(self.design_model_parts)>0 :
            raise RuntimeError("StrainEnergyResponseFunction: 'design_model_parts' of response '{}' can not be empty !".format(self.name))

        for design_model_part in self.design_model_parts:
            design_model_part_splitted = design_model_part.split(".")
            if not design_model_part_splitted[0] == self.primal_model_part.Name:
                raise RuntimeError("StrainEnergyResponseFunction: root design_model_part {} of response '{}' is not the analysis model!".format(design_model_part_splitted[0],self.name))
            if not self.model.HasModelPart(design_model_part): 
                raise RuntimeError("StrainEnergyResponseFunction: evaluate_model_part {} of response '{}' does not exist!".format(design_model_part,self.name))

        for design_type in self.design_types:
            if not design_type in self.supported_design_types:
                raise RuntimeError("StrainEnergyResponseFunction: design type {} of response '{}' is not supported !".format(design_type,self.name))

        self.response_function_utility.Initialize()

    def CalculateValue(self):
        return self.response_function_utility.CalculateValue()

    def CalculateGradient(self):
        pass
        # Logger.PrintInfo("StrainEnergyResponse", "Starting gradient calculation for response", self.identifier)

        # startTime = timer.time()
        # self.response_function_utility.CalculateGradient()
        # Logger.PrintInfo("StrainEnergyResponse", "Time needed for calculating gradients",round(timer.time() - startTime,2),"s")

    def FinalizeSolutionStep(self):
        pass
        # self.primal_analysis.FinalizeSolutionStep()
        # self.primal_analysis.OutputSolutionStep()

    def Finalize(self):
        pass
        # self.primal_analysis.Finalize()

    def GetValue(self):
        pass
        # return self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE]

    def GetNodalGradient(self, variable):
        pass
        # if variable != KratosMultiphysics.SHAPE_SENSITIVITY:
        #     raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))
        # gradient = {}
        # for node in self.primal_model_part.Nodes:
        #     gradient[node.Id] = node.GetSolutionStepValue(variable)
        # return gradient
