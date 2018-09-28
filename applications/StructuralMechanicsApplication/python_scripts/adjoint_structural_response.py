"""This module contains the available structural response functions and their base class"""
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
import structural_mechanics_analysis

import time as timer

def _GetModelPart(model, solver_settings):
    #TODO can be removed once model is fully available
    model_part_name = solver_settings["model_part_name"].GetString()
    if not model.HasModelPart(model_part_name):
        model_part = ModelPart(model_part_name)
        domain_size = solver_settings["domain_size"].GetInt()
        if domain_size < 0:
            raise Exception('Please specify a "domain_size" >= 0!')
        model_part.ProcessInfo.SetValue(DOMAIN_SIZE, domain_size)
        model.AddModelPart(model_part)
    else:
        model_part = model.GetModelPart(model_part_name)

    return model_part

# ==============================================================================
class AdjointResponseFunctionBase(object):
    """The base class for structural response functions. Each response function
    is able to calculate its response value and gradient.
    All the necessary steps have to be implemented, like e.g. initializing,
    solving of primal (and adjoint) analysis ...
    """

    def RunCalculation(self, calculate_gradient):
        self.Initialize()
        self.InitializeSolutionStep()
        self.SolveAdjointProblem()
        if calculate_gradient:
            self.PerformSensitivityPostprocess()
        self.FinalizeSolutionStep()
        self.Finalize()

    def Initialize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def SolveAdjointProblem(self):
        raise NotImplementedError("SolveAdjointProblem needs to be implemented by the derived class")

    def PerformSensitivityPostprocess(self):
        raise NotImplementedError("PerformSensitivityPostprocess needs to be implemented by the derived class")

    def FinalizeSolutionStep(self):
        pass

    def Finalize(self):
        pass

# ==============================================================================
class AdjointResponseFunction(AdjointResponseFunctionBase):
    """Linear static adjointresponse function.
    - runs the primal analysis (writes the primal results to an .h5 file)
    - reads the primal results from the .h5 file into the adjoint model part
    - uses primal results to calculate value
    - uses primal results to calculate gradient by running the adjoint analysis

    Attributes
    ----------
    primal_analysis : Primal analysis object of the response function
    adjoint_analysis : Adjoint analysis object of the response function
    """
    def __init__(self, identifier, project_parameters, model, solve_primal = True):
        self.identifier = identifier

        self.solve_primal = solve_primal

        if solve_primal:
            # Create the primal solver
            with open(project_parameters["primal_settings"].GetString(),'r') as parameter_file:
                ProjectParametersPrimal = Parameters( parameter_file.read() )

            self.primal_model_part = _GetModelPart(model, ProjectParametersPrimal["solver_settings"])

            self.primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model, ProjectParametersPrimal)

        # Create the adjoint solver
        with open(project_parameters["adjoint_settings"].GetString(),'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read() )
        ProjectParametersAdjoint["solver_settings"].AddValue("response_function_settings", project_parameters)

        adjoint_model = Model()

        self.adjoint_model_part = _GetModelPart(adjoint_model, ProjectParametersAdjoint["solver_settings"])

        # TODO find out why it is not possible to use the same model_part
        self.adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(adjoint_model, ProjectParametersAdjoint)

        self.response_function_settings = project_parameters.Clone()

    def Initialize(self):
        if self.solve_primal:
            self.primal_analysis.Initialize()
        self.adjoint_analysis.Initialize()
        response_function = self._GetResponseFunctionUtility()
        self.adjoint_postprocess = StructuralMechanicsApplication.AdjointPostprocess(self.adjoint_model_part, response_function, self.response_function_settings)

    def InitializeSolutionStep(self):

        # Run the primal analysis if necessary.
        if self.solve_primal:
            Logger.PrintInfo("\n> Starting primal analysis for response:", self.identifier)
            startTime = timer.time()
            if not self.primal_analysis.time < self.primal_analysis.end_time:
                self.primal_analysis.end_time += 1
            self.primal_analysis.RunSolutionLoop()
            Logger.PrintInfo("> Time needed for solving the primal analysis = ",round(timer.time() - startTime,2),"s")

        # TODO the response value calculation for stresses currently only works on the adjoint modelpart
        # this needs to be improved, also the response value should be calculated on the PRIMAL modelpart!!
        self.adjoint_analysis.time = self.adjoint_analysis._GetSolver().AdvanceInTime(self.adjoint_analysis.time)
        self.adjoint_analysis.InitializeSolutionStep()

    def SolveAdjointProblem(self):
        Logger.PrintInfo("\n> Starting to solve adjoint problem for response:", self.identifier)
        startTime = timer.time()
        self.adjoint_analysis._GetSolver().Predict()
        self.adjoint_analysis._GetSolver().SolveSolutionStep()
        self.adjoint_postprocess.UpdateSensitivities()
        Logger.PrintInfo("> Time needed for solving the adjoint problem = ",round(timer.time() - startTime,2),"s")

    def PerformSensitivityPostprocess(self):
        Logger.PrintInfo("\n> Starting to compute sensivities in postprocessing step:", self.identifier)
        startTime = timer.time()
        self.adjoint_postprocess.UpdateSensitivities()
        Logger.PrintInfo("> Time needed for computing adjoint sensivities = ",round(timer.time() - startTime,2),"s")

    def FinalizeSolutionStep(self):
        self.adjoint_analysis.FinalizeSolutionStep()
        self.adjoint_analysis.OutputSolutionStep()

    def Finalize(self):
        if self.solve_primal:
            self.primal_analysis.Finalize()
        self.adjoint_analysis.Finalize()

    def _GetResponseFunctionUtility(self):
        return self.adjoint_analysis._GetSolver().response_function



