"""This module contains the available structural response functions and their base class"""
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
import structural_mechanics_analysis

import time as timer

class ResponseFunctionBase(object):
    """The base class for structural response functions. Each response function
    is able to calculate its response value and gradient.
    All the necessary steps have to be implemented, like e.g. initializing,
    solving of primal (and adjoint) analysis ...
    """

    def Initialize(self):
        pass

    def CalculateValue(self):
        raise NotImplementedError("CalculateValue needs to be implemented by the base class")

    def CalculateGradient(self):
        raise NotImplementedError("CalculateGradient needs to be implemented by the base class")

    def GetShapeGradient(self):
        raise NotImplementedError("GetShapeGradient needs to be implemented by the base class")

    def Finalize(self):
        pass

class StrainEnergyResponseFunction(ResponseFunctionBase):
    """Linear strain energy response function. It triggers the primal analysis and
    uses the primal analysis results to evaluate response value and gradient.

    Attributes
    ----------
    primal_analysis : Primal analysis object of the response function
    response_function_utility: Cpp utilities object doing the actual computation of response value and gradient.
    """

    def __init__(self, identifier, response_settings, model_part = None):
        self.identifier = identifier

        self.response_function_utility = StructuralMechanicsApplication.StrainEnergyResponseFunctionUtility(model_part, response_settings)

        with open(response_settings["primal_settings"].GetString()) as parameters_file:
            ProjectParametersPrimal = Parameters(parameters_file.read())

        self.primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersPrimal, model_part)
        self.primal_analysis.GetModelPart().AddNodalSolutionStepVariable(SHAPE_SENSITIVITY)

    def Initialize(self):
        self.primal_analysis.Initialize()
        self.response_function_utility.Initialize()

    def CalculateValue(self):
        Logger.PrintInfo("\n> Starting primal analysis for response:", self.identifier)

        startTime = timer.time()
        self.primal_analysis.InitializeTimeStep()
        self.primal_analysis.SolveTimeStep()
        self.primal_analysis.FinalizeTimeStep()
        Logger.PrintInfo("> Time needed for solving the primal analysis = ",round(timer.time() - startTime,2),"s")

        startTime = timer.time()
        value = self.response_function_utility.CalculateValue()
        Logger.PrintInfo("> Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

        return value

    def CalculateGradient(self):
        self.response_function_utility.CalculateGradient()

    def GetShapeGradient(self):
        gradient = {}
        for node in self.primal_analysis.GetModelPart().Nodes:
            gradient[node.Id] = node.GetSolutionStepValue(SHAPE_SENSITIVITY)
        return gradient

    def Finalize(self):
        self.primal_analysis.Finalize()

class EigenFrequencyResponseFunction(StrainEnergyResponseFunction):
    """Eigenfrequency response function. The internal procedure is the same as
    for the StrainEnergyResponseFunction. It triggers the primal analysis and
    uses the primal analysis results to evaluate response value and gradient.
    Only the response_function_utility is a different object.

    Attributes
    ----------
    primal_analysis : Primal analysis object of the response function
    response_function_utility: Cpp utilities object doing the actual computation of response value and gradient.
    """

    def __init__(self, identifier, response_settings, model_part = None):
        self.identifier = identifier

        self.response_function_utility = StructuralMechanicsApplication.EigenfrequencyResponseFunctionUtility(model_part, response_settings)

        with open(response_settings["primal_settings"].GetString()) as parameters_file:
            ProjectParametersPrimal = Parameters(parameters_file.read())

        max_required_eigenfrequency = int(max(response_settings["traced_eigenfrequencies"].GetVector()))
        if max_required_eigenfrequency is not ProjectParametersPrimal["solver_settings"]["eigensolver_settings"]["number_of_eigenvalues"].GetInt():
            print("\n> WARNING: Specified number of eigenvalues in the primal analysis and the max required eigenvalue according the response settings do not match!!!")
            print("  Primal parameters were adjusted accordingly!\n")
            ProjectParametersPrimal["solver_settings"]["eigensolver_settings"]["number_of_eigenvalues"].SetInt(max_required_eigenfrequency)

        self.primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersPrimal, model_part)
        self.primal_analysis.GetModelPart().AddNodalSolutionStepVariable(SHAPE_SENSITIVITY)

class MassResponseFunction(ResponseFunctionBase):
    """Mass response function. It reads the materials for the model part and
    calculates response value and gradient.

    Attributes
    ----------
    model_part : Model part object of the response function
    response_function_utility: Cpp utilities object doing the actual computation of response value and gradient.
    """

    def __init__(self, identifier, response_settings, model_part):
        self.identifier = identifier
        self.response_settings = response_settings

        self.response_function_utility = StructuralMechanicsApplication.MassResponseFunctionUtility(model_part, response_settings)

        self.model_part = model_part
        self.model_part.AddNodalSolutionStepVariable(SHAPE_SENSITIVITY)

    def Initialize(self):
        import read_materials_process
        # Create a dictionary of model parts.
        model = Model()
        model.AddModelPart(self.model_part)
        # Add constitutive laws and material properties from json file to model parts.
        read_materials_process.ReadMaterialsProcess(model, self.response_settings["material_import_settings"])
        self.response_function_utility.Initialize()

    def CalculateValue(self):
        value = self.response_function_utility.CalculateValue()
        return value

    def CalculateGradient(self):
        self.response_function_utility.CalculateGradient()

    def GetShapeGradient(self):
        gradient = {}
        for node in self.model_part.Nodes:
            gradient[node.Id] = node.GetSolutionStepValue(SHAPE_SENSITIVITY)
        return gradient




