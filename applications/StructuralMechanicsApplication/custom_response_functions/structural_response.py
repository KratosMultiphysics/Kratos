from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
import structural_mechanics_analysis

import time as timer

class ResponseFunctionBase(object):
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

class AdjointResponseFunction(ResponseFunctionBase):
    def __init__(self, identifier, project_parameters, model_part = None):
        # TODO use **kwargs for optional model_part, primal_analysis, adjoint_analysis
        # TODO remember the name of the hdf5 file written
        super(AdjointResponseFunction, self).__init__(identifier, project_parameters)

        # Create the primal solver
        ProjectParametersPrimal = Parameters( open(project_parameters["primal_settings"].GetString(),'r').read() )
        self.primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersPrimal, model_part)

        # Create the adjoint solver
        ProjectParametersAdjoint = Parameters( open(project_parameters["adjoint_settings"].GetString(),'r').read() )
        ProjectParametersAdjoint["solver_settings"].AddValue("response_function_settings", project_parameters)
        self.adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersAdjoint)
        # TODO find out why it is not possible to use the same model_part
    def Initialize(self):
        self.primal_analysis.Initialize()
        self.adjoint_analysis.Initialize()

        # TODO should be created here, not in solver!
        self.response_function_utility = self.adjoint_analysis.GetSolver().response_function
    def CalculateValue(self):
        for node in self.adjoint_analysis.GetModelPart().Nodes:
            ref_node = self.primal_analysis.GetModelPart().Nodes[node.Id]
            node.X0 = ref_node.X0
            node.Y0 = ref_node.Y0
            node.Z0 = ref_node.Z0
            node.X = ref_node.X
            node.Y = ref_node.Y
            node.Z = ref_node.Z

        print("\n> Starting primal analysis for response:", self.identifier)
        startTime = timer.time()
        self.primal_analysis.InitializeTimeStep()
        self.primal_analysis.SolveTimeStep()
        self.primal_analysis.FinalizeTimeStep()
        print("> Time needed for solving the primal analysis = ",round(timer.time() - startTime,2),"s")
        startTime = timer.time()
        value = self.response_function_utility.CalculateValue(self.primal_analysis.GetModelPart())
        print("> Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")
        return value
    def CalculateGradient(self):
        print("\n> Starting adjoint analysis for response:", self.identifier)
        startTime = timer.time()
        self.adjoint_analysis.InitializeTimeStep()
        self.adjoint_analysis.SolveTimeStep()
        self.adjoint_analysis.FinalizeTimeStep()
        print("> Time needed for solving the adjoint analysis = ",round(timer.time() - startTime,2),"s")
    def GetShapeGradient(self):
        gradient = {}
        for node in self.adjoint_analysis.GetModelPart().Nodes:
            gradient[node.Id] = node.GetSolutionStepValue(SHAPE_SENSITIVITY)
        return gradient
        # TODO reset DISPLACEMENT, ROTATION ADJOINT_DISPLACEMENT and ADJOINT_ROTATION
    def Finalize(self):
        self.primal_analysis.Finalize()
        self.adjoint_analysis.Finalize()

class SimpleResponseFunctionWrapper(ResponseFunctionBase):
    def __init__(self, identifier, project_parameters, response_function_utility, model_part = None):
        super(SimpleResponseFunctionWrapper, self).__init__(identifier, project_parameters)

        # Create the primal solver
        ProjectParametersPrimal = Parameters( open(project_parameters["primal_settings"].GetString(),'r').read() )
        self.primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersPrimal, model_part)

        self.primal_analysis.GetModelPart().AddNodalSolutionStepVariable(SHAPE_SENSITIVITY)

        self.response_function_utility = response_function_utility
    def Initialize(self):
        self.primal_analysis.Initialize()
        self.response_function_utility.Initialize()
    def CalculateValue(self):
        print("\n> Starting primal analysis for response:", self.identifier)
        startTime = timer.time()
        self.primal_analysis.InitializeTimeStep()
        self.primal_analysis.SolveTimeStep()
        self.primal_analysis.FinalizeTimeStep()
        print("> Time needed for solving the primal analysis = ",round(timer.time() - startTime,2),"s")
        startTime = timer.time()
        value = self.response_function_utility.CalculateValue()
        print("> Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")
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

class MassResponseFunctionWrapper(ResponseFunctionBase):
    def __init__(self, identifier, project_parameters, response_function_utility, model_part):
        super(MassResponseFunctionWrapper, self).__init__(identifier, project_parameters)
        self.response_function_utility = response_function_utility
        self.model_part = model_part
        self.model_part.AddNodalSolutionStepVariable(SHAPE_SENSITIVITY)
    def Initialize(self):
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
    def Finalize(self):
        pass
