# importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics import Parameters, Logger
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

import time as timer

def _GetModelPart(model, solver_settings):
    #TODO can be removed once model is fully available
    model_part_name = solver_settings["model_part_name"].GetString()
    if not model.HasModelPart(model_part_name):
        model_part = model.CreateModelPart(model_part_name, 2)
        domain_size = solver_settings["domain_size"].GetInt()
        if domain_size < 0:
            raise Exception('Please specify a "domain_size" >= 0!')
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)
    else:
        model_part = model.GetModelPart(model_part_name)

    return model_part

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

    def __init__(self, response_settings,response_analyzer,response_analyzer_model_part,model):

        self.response_settings = response_settings
        default_gradient_settings = KM.Parameters("""
        {
            "gradient_mode" : "semi_analytic",
            "step_size" : 1e-6
        }""")
        
        self.response_settings["gradient_settings"].ValidateAndAssignDefaults(default_gradient_settings)        

        self.primal_analysis = response_analyzer
        self.model = model
        self.primal_model_part = response_analyzer_model_part
        self.response_function_utility = StructuralMechanicsApplication.StrainEnergyResponseFunctionUtility(self.primal_model_part, self.response_settings["gradient_settings"])

    def Initialize(self):
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
