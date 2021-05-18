"""
This module contains an interface to the available response functions
"""
import time as timer

import KratosMultiphysics as KM
from KratosMultiphysics import Parameters, Logger
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis


class AdjointResponseFunction(ResponseFunctionInterface):
    """Linear adjoint response function.
    - runs the primal analysis
    - primal results are transferred to adjoint model part via python
    - uses primal results to calculate value
    - uses primal results to calculate gradient by running the adjoint analysis

    Attributes
    ----------
    primal_analysis : Primal analysis object of the response function
    adjoint_analysis : Adjoint analysis object of the response function
    """
    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings

        # Create the primal solver
        with open(self.response_settings["primal_settings"].GetString(),'r') as parameter_file:
            primal_parameters = Parameters( parameter_file.read() )

        self.primal_analysis = ConvectionDiffusionAnalysis(model, primal_parameters)
        self.primal_model_part = model.GetModelPart(primal_parameters["solver_settings"]["model_part_name"].GetString())

        # Create the adjoint solver
        adjoint_parameters = self._GetAdjointParameters()

        self.adjoint_analysis = ConvectionDiffusionAnalysis(model, adjoint_parameters)
        self.adjoint_model_part = model.GetModelPart(adjoint_parameters["solver_settings"]["model_part_name"].GetString())

        self.primal_state_variables = [
            KM.CONDUCTIVITY,
            KM.TEMPERATURE,
            KM.HEAT_FLUX,
            KM.FACE_HEAT_FLUX
        ]

        self.value = None

    def Initialize(self):
        self.primal_analysis.Initialize()
        self.adjoint_analysis.Initialize()

    def InitializeSolutionStep(self):
        self.value = None

        # Run the primal analysis.
        Logger.PrintInfo(self._GetLabel(), "Starting primal analysis for response:", self.identifier)
        startTime = timer.time()
        if not self.primal_analysis.time < self.primal_analysis.end_time:
            self.primal_analysis.end_time += 1
        self.primal_analysis.RunSolutionLoop()
        Logger.PrintInfo(self._GetLabel(), "Time needed for solving the primal analysis = ",round(timer.time() - startTime,2),"s")

    def CalculateValue(self):
        startTime = timer.time()
        self.value = self._GetResponseFunctionUtility().CalculateValue(self.primal_model_part)
        Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        # synchronize the modelparts
        self._SynchronizeAdjointFromPrimal()
        startTime = timer.time()
        Logger.PrintInfo(self._GetLabel(), "Starting adjoint analysis for response:", self.identifier)
        if not self.adjoint_analysis.time < self.adjoint_analysis.end_time:
            self.adjoint_analysis.end_time += 1
        self.adjoint_analysis.RunSolutionLoop()
        Logger.PrintInfo(self._GetLabel(), "Time needed for solving the adjoint analysis = ",round(timer.time() - startTime,2),"s")

    def GetValue(self):
        return self.value

    def GetNodalGradient(self, variable):
        if variable != KM.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))
        gradient = {}
        for node in self.adjoint_model_part.Nodes:
            gradient[node.Id] = node.GetSolutionStepValue(variable)
        return gradient

    def Finalize(self):
        self.primal_analysis.Finalize()
        self.adjoint_analysis.Finalize()

    def _GetResponseFunctionUtility(self):
        return self.adjoint_analysis._GetSolver().response_function

    def _SynchronizeAdjointFromPrimal(self):
        Logger.PrintInfo(self._GetLabel(), "Synchronize primal and adjoint modelpart for response:", self.identifier)

        if len(self.primal_model_part.Nodes) != len(self.adjoint_model_part.Nodes):
            raise RuntimeError("_SynchronizeAdjointFromPrimal: Model parts have a different number of nodes!")

        for primal_node, adjoint_node in zip(self.primal_model_part.Nodes, self.adjoint_model_part.Nodes):
            adjoint_node.X0 = primal_node.X0
            adjoint_node.Y0 = primal_node.Y0
            adjoint_node.Z0 = primal_node.Z0
            adjoint_node.X = primal_node.X
            adjoint_node.Y = primal_node.Y
            adjoint_node.Z = primal_node.Z

        # Put primal solution on adjoint model
        Logger.PrintInfo(self._GetLabel(), "Transfer primal state to adjoint model part.")
        variable_utils = KM.VariableUtils()
        for variable in self.primal_state_variables:
            variable_utils.CopyModelPartNodalVar(variable, self.primal_model_part, self.adjoint_model_part, 0)

    def _GetAdjointParameters(self):

        with open(self.response_settings["adjoint_settings"].GetString(),'r') as parameter_file:
            adjoint_parameters = Parameters( parameter_file.read() )

        if self.response_settings["transfer_response_parameters"].GetBool():

            # sensitivity settings
            if adjoint_parameters["solver_settings"].Has("sensitivity_settings"):
                Logger.PrintWarning(self._GetLabel(), "Adjoint settings already have 'sensitivity_settings' - Will be overwritten!")
                adjoint_parameters["solver_settings"].RemoveValue("sensitivity_settings")
            adjoint_parameters["solver_settings"].AddValue("sensitivity_settings", self.response_settings["sensitivity_settings"])

            # response settings
            if adjoint_parameters["solver_settings"].Has("response_function_settings"):
                Logger.PrintWarning(self._GetLabel(), "Adjoint settings already have 'response_function_settings' - Will be overwritten!")
                adjoint_parameters["solver_settings"].RemoveValue("response_function_settings")
            adjoint_parameters["solver_settings"].AddValue("response_function_settings", self.response_settings)
        
        return adjoint_parameters

    def _GetLabel(self):
        type_labels = {
            "point_temperature" : "LocalTemperatureAverageResponseFunction"
        }
        response_type = self.response_settings["response_type"].GetString()
        return "Adjoint" + type_labels[response_type] + "Response"
