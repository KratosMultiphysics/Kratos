"""This module contains the available structural response functions and their base class"""
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
import structural_mechanics_analysis

import time as timer

# ==============================================================================
class ResponseFunctionBase(object):
    """The base class for structural response functions. Each response function
    is able to calculate its response value and gradient.
    All the necessary steps have to be implemented, like e.g. initializing,
    solving of primal (and adjoint) analysis ...
    """

    def RunCalculation(self, calculate_gradient):
        self.Initialize()
        self.InitializeSolutionStep()
        self.CalculateValue()
        if calculate_gradient:
            self.CalculateGradient()
        self.FinalizeSolutionStep()
        self.Finalize()

    def Initialize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def CalculateValue(self):
        raise NotImplementedError("CalculateValue needs to be implemented by the derived class")

    def CalculateGradient(self):
        raise NotImplementedError("CalculateGradient needs to be implemented by the derived class")

    def FinalizeSolutionStep(self):
        pass

    def Finalize(self):
        pass

    def GetValue(self):
        raise NotImplementedError("GetValue needs to be implemented by the derived class")

    def GetShapeGradient(self):
        raise NotImplementedError("GetShapeGradient needs to be implemented by the derived class")

# ==============================================================================
class StrainEnergyResponseFunction(ResponseFunctionBase):
    """Linear strain energy response function. It triggers the primal analysis and
    uses the primal analysis results to calculate response value and gradient.

    Attributes
    ----------
    primal_model_part : Model part of the primal analysis object
    primal_analysis : Primal analysis object of the response function
    response_function_utility: Cpp utilities object doing the actual computation of response value and gradient.
    """

    def __init__(self, identifier, response_settings, model_part = None):
        self.identifier = identifier

        self.response_function_utility = StructuralMechanicsApplication.StrainEnergyResponseFunctionUtility(model_part, response_settings)

        with open(response_settings["primal_settings"].GetString()) as parameters_file:
            ProjectParametersPrimal = Parameters(parameters_file.read())


        self.primal_model_part = model_part
        model = Model()
        model.AddModelPart(self.primal_model_part)
        self.primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model, ProjectParametersPrimal)
        self.primal_model_part.AddNodalSolutionStepVariable(SHAPE_SENSITIVITY)

    def Initialize(self):
        self.primal_analysis.Initialize()
        self.response_function_utility.Initialize()

    def InitializeSolutionStep(self):
        self.primal_analysis.time = self.primal_analysis.solver.AdvanceInTime(self.primal_analysis.time)
        self.primal_analysis.InitializeSolutionStep()

    def CalculateValue(self):
        Logger.PrintInfo("\n> Starting primal analysis for response", self.identifier)

        startTime = timer.time()
        self.primal_analysis.solver.Predict()
        self.primal_analysis.solver.SolveSolutionStep()
        Logger.PrintInfo("> Time needed for solving the primal analysis",round(timer.time() - startTime,2),"s")

        startTime = timer.time()
        value = self.response_function_utility.CalculateValue()
        self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE] = value
        Logger.PrintInfo("> Time needed for calculating the response value",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        Logger.PrintInfo("\n> Starting gradient calculation for response", self.identifier)

        startTime = timer.time()
        self.response_function_utility.CalculateGradient()
        Logger.PrintInfo("> Time needed for calculating gradients",round(timer.time() - startTime,2),"s")

    def FinalizeSolutionStep(self):
        self.primal_analysis.FinalizeSolutionStep()
        self.primal_analysis.OutputSolutionStep()

    def Finalize(self):
        self.primal_analysis.Finalize()

    def GetValue(self):
        return self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE]

    def GetShapeGradient(self):
        gradient = {}
        for node in self.primal_model_part.Nodes:
            gradient[node.Id] = node.GetSolutionStepValue(SHAPE_SENSITIVITY)
        return gradient

# ==============================================================================
class EigenFrequencyResponseFunction(StrainEnergyResponseFunction):
    """Eigenfrequency response function. The internal procedure is the same as
    for the StrainEnergyResponseFunction. It triggers the primal analysis and
    uses the primal analysis results to calculate response value and gradient.
    Only the response_function_utility is a different object.

    Attributes
    ----------
    primal_model_part : Model part of the primal analysis object
    primal_analysis : Primal analysis object of the response function
    response_function_utility: Cpp utilities object doing the actual computation of response value and gradient.
    """

    def __init__(self, identifier, response_settings, model_part):
        self.identifier = identifier

        self.primal_model_part = model_part
        self.response_function_utility = StructuralMechanicsApplication.EigenfrequencyResponseFunctionUtility(model_part, response_settings)

        with open(response_settings["primal_settings"].GetString()) as parameters_file:
            ProjectParametersPrimal = Parameters(parameters_file.read())

        eigen_solver_settings = ProjectParametersPrimal["solver_settings"]["eigensolver_settings"]

        max_required_eigenfrequency = int(max(response_settings["traced_eigenfrequencies"].GetVector()))
        if max_required_eigenfrequency is not eigen_solver_settings["number_of_eigenvalues"].GetInt():
            print("\n> WARNING: Specified number of eigenvalues in the primal analysis and the max required eigenvalue according the response settings do not match!!!")
            print("  Primal parameters were adjusted accordingly!\n")
            eigen_solver_settings["number_of_eigenvalues"].SetInt(max_required_eigenfrequency)

        if not eigen_solver_settings.Has("normalize_eigenvectors"):
            eigen_solver_settings.AddEmptyValue("normalize_eigenvectors")
            eigen_solver_settings["normalize_eigenvectors"].SetBool(True)
            print("\n> WARNING: Eigenfrequency response function requires mass normalization of eigenvectors!")
            print("  Primal parameters were adjusted accordingly!\n")

        if not eigen_solver_settings["normalize_eigenvectors"].GetBool():
            eigen_solver_settings["normalize_eigenvectors"].SetBool(True)
            print("\n> WARNING: Eigenfrequency response function requires mass normalization of eigenvectors!")
            print("  Primal parameters were adjusted accordingly!\n")

        model = Model()
        model.AddModelPart(self.primal_model_part)
        self.primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model, ProjectParametersPrimal)
        self.primal_model_part.AddNodalSolutionStepVariable(SHAPE_SENSITIVITY)

# ==============================================================================
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
        Logger.PrintInfo("\n> Starting primal analysis for response", self.identifier)

        startTime = timer.time()
        value = self.response_function_utility.CalculateValue()
        self.model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE] = value
        Logger.PrintInfo("> Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

    def CalculateGradient(self):
        Logger.PrintInfo("\n> Starting gradient calculation for response", self.identifier)

        startTime = timer.time()
        self.response_function_utility.CalculateGradient()
        Logger.PrintInfo("> Time needed for calculating gradients",round(timer.time() - startTime,2),"s")

    def GetValue(self):
        return self.model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE]

    def GetShapeGradient(self):
        gradient = {}
        for node in self.model_part.Nodes:
            gradient[node.Id] = node.GetSolutionStepValue(SHAPE_SENSITIVITY)
        return gradient

# ==============================================================================




