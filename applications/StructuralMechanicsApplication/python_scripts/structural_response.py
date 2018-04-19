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

    def __init__(self, identifier, project_parameters):
        self.identifier = identifier
        self.project_parameters = project_parameters

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
    primal_analysis : Object
        Primal analysis object of the response function
    response_function_utility: Object
        Cpp utilities object that does the actual computation of response value and
        gradient.
    """

    def __init__(self, identifier, project_parameters, response_function_utility, model_part = None):
        super(StrainEnergyResponseFunction, self).__init__(identifier, project_parameters)

        with open(project_parameters["primal_settings"].GetString()) as parameters_file:
            ProjectParametersPrimal = Parameters( parameters_file.read() )
        self.primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersPrimal, model_part)
        self.primal_analysis.GetModelPart().AddNodalSolutionStepVariable(SHAPE_SENSITIVITY)
        self.response_function_utility = response_function_utility

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

    Attributes
    ----------
    primal_analysis : Object
        Primal analysis object of the response function
    response_function_utility: Object
        Cpp utilities object that does the actual computation of response value and
        gradient.
    """

    def __init__(self, identifier, project_parameters, response_function_utility, model_part = None):
        super(EigenFrequencyResponseFunction, self).__init__(identifier, project_parameters, response_function_utility, model_part)

class MassResponseFunction(ResponseFunctionBase):
    """Mass response function. It reads the materials for the model part and
    calculates response value and gradient.

    Attributes
    ----------
    model_part : Object
        Model part of the response function
    response_function_utility: Object
        Cpp utilities object that does the actual computation of response value and
        gradient.
    """

    def __init__(self, identifier, project_parameters, response_function_utility, model_part):
        super(MassResponseFunction, self).__init__(identifier, project_parameters)
        self.response_function_utility = response_function_utility
        self.model_part = model_part
        self.model_part.AddNodalSolutionStepVariable(SHAPE_SENSITIVITY)

    def Initialize(self):
        import read_materials_process
        # Create a dictionary of model parts.
        model = Model()
        model.AddModelPart(self.model_part)
        # Add constitutive laws and material properties from json file to model parts.
        read_materials_process.ReadMaterialsProcess(model, self.project_parameters["material_import_settings"])
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




