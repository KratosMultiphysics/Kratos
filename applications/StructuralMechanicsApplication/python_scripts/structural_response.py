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

    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier

        with open(response_settings["primal_settings"].GetString()) as parameters_file:
            ProjectParametersPrimal = Parameters(parameters_file.read())

        self.primal_model_part = _GetModelPart(model, ProjectParametersPrimal["solver_settings"])

        self.primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model, ProjectParametersPrimal)
        self.primal_model_part.AddNodalSolutionStepVariable(SHAPE_SENSITIVITY)

        self.response_function_utility = StructuralMechanicsApplication.StrainEnergyResponseFunctionUtility(self.primal_model_part, response_settings)

    def Initialize(self):
        self.primal_analysis.Initialize()
        self.response_function_utility.Initialize()

    def InitializeSolutionStep(self):
        self.primal_analysis.time = self.primal_analysis._GetSolver().AdvanceInTime(self.primal_analysis.time)
        self.primal_analysis.InitializeSolutionStep()

    def CalculateValue(self):
        Logger.PrintInfo("\n> Starting primal analysis for response", self.identifier)

        startTime = timer.time()
        self.primal_analysis._GetSolver().Predict()
        self.primal_analysis._GetSolver().SolveSolutionStep()
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

    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier

        with open(response_settings["primal_settings"].GetString()) as parameters_file:
            ProjectParametersPrimal = Parameters(parameters_file.read())

        eigen_solver_settings = ProjectParametersPrimal["solver_settings"]["eigensolver_settings"]

        max_required_eigenfrequency = int(max(response_settings["traced_eigenfrequencies"].GetVector()))
        if max_required_eigenfrequency is not eigen_solver_settings["number_of_eigenvalues"].GetInt():
            Logger.PrintWarning("\n> WARNING: Specified number of eigenvalues in the primal analysis and the max required eigenvalue according the response settings do not match!!!")
            Logger.PrintWarning("  Primal parameters were adjusted accordingly!\n")
            eigen_solver_settings["number_of_eigenvalues"].SetInt(max_required_eigenfrequency)

        if not eigen_solver_settings.Has("normalize_eigenvectors"):
            eigen_solver_settings.AddEmptyValue("normalize_eigenvectors")
            eigen_solver_settings["normalize_eigenvectors"].SetBool(True)
            Logger.PrintWarning("\n> WARNING: Eigenfrequency response function requires mass normalization of eigenvectors!")
            Logger.PrintWarning("  Primal parameters were adjusted accordingly!\n")

        if not eigen_solver_settings["normalize_eigenvectors"].GetBool():
            eigen_solver_settings["normalize_eigenvectors"].SetBool(True)
            Logger.PrintWarning("\n> WARNING: Eigenfrequency response function requires mass normalization of eigenvectors!")
            Logger.PrintWarning("  Primal parameters were adjusted accordingly!\n")

        self.primal_model_part = _GetModelPart(model, ProjectParametersPrimal["solver_settings"])

        self.primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model, ProjectParametersPrimal)
        self.primal_model_part.AddNodalSolutionStepVariable(SHAPE_SENSITIVITY)

        self.response_function_utility = StructuralMechanicsApplication.EigenfrequencyResponseFunctionUtility(self.primal_model_part, response_settings)

# ==============================================================================
class MassResponseFunction(ResponseFunctionBase):
    """Mass response function. It reads the materials for the model part and
    calculates response value and gradient.

    Attributes
    ----------
    model_part : Model part object of the response function
    response_function_utility: Cpp utilities object doing the actual computation of response value and gradient.
    """

    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings
        self.model = model
        self.model_part_needs_to_be_imported = False

        input_type = response_settings["model_import_settings"]["input_type"].GetString()
        model_part_name = response_settings["model_import_settings"]["input_filename"].GetString()
        if input_type == "mdpa":
            self.model_part = ModelPart(model_part_name)
            self.model.AddModelPart(self.model_part)
            self.model_part_needs_to_be_imported = True
        elif input_type == "use_input_model_part":
            self.model_part = self.model.GetModelPart(model_part_name)
        else:
            raise Exception("Other model part input options are not yet implemented.")

        self.response_function_utility = StructuralMechanicsApplication.MassResponseFunctionUtility(self.model_part, response_settings)

        self.model_part.AddNodalSolutionStepVariable(SHAPE_SENSITIVITY)

    def Initialize(self):
        import read_materials_process

        if self.model_part_needs_to_be_imported:
            # import model part
            model_part_io = ModelPartIO(self.model_part.Name)
            self.model_part.ProcessInfo.SetValue(DOMAIN_SIZE, 3)
            model_part_io.ReadModelPart(self.model_part)

        # Add constitutive laws and material properties from json file to model parts.
        read_materials_process.ReadMaterialsProcess(self.model, self.response_settings["material_import_settings"])
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
class AdjointResponseFunction(ResponseFunctionBase):
    """Linear static adjoint strain energy response function.
    - runs the primal analysis (writes the primal results to an .h5 file)
    - reads the primal results from the .h5 file into the adjoint model part
    - uses primal results to calculate value
    - uses primal results to calculate gradient by running the adjoint analysis

    Attributes
    ----------
    primal_analysis : Primal analysis object of the response function
    adjoint_analysis : Adjoint analysis object of the response function
    """
    def __init__(self, identifier, project_parameters, model):
        self.identifier = identifier

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

    def Initialize(self):
        self.primal_analysis.Initialize()
        self.adjoint_analysis.Initialize()

    def InitializeSolutionStep(self):
        # synchronize the modelparts # TODO this should happen automatically
        Logger.PrintInfo("\n> Synchronize primal and adjoint modelpart for response:", self.identifier)

        self._SynchronizeAdjointFromPrimal()

        # Run the primal analysis.
        # TODO if primal_analysis.status==solved: return
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


    def CalculateValue(self):
        startTime = timer.time()
        value = self._GetResponseFunctionUtility().CalculateValue(self.adjoint_model_part)
        Logger.PrintInfo("> Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

        self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE] = value


    def CalculateGradient(self):
        Logger.PrintInfo("\n> Starting adjoint analysis for response:", self.identifier)
        startTime = timer.time()
        self.adjoint_analysis._GetSolver().Predict()
        self.adjoint_analysis._GetSolver().SolveSolutionStep()
        Logger.PrintInfo("> Time needed for solving the adjoint analysis = ",round(timer.time() - startTime,2),"s")


    def GetValue(self):
        return self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE]


    def GetShapeGradient(self):
        gradient = {}
        for node in self.adjoint_model_part.Nodes:
            gradient[node.Id] = node.GetSolutionStepValue(SHAPE_SENSITIVITY)
        return gradient


    def FinalizeSolutionStep(self):
        self.adjoint_analysis.FinalizeSolutionStep()
        self.adjoint_analysis.OutputSolutionStep()


    def Finalize(self):
        self.primal_analysis.Finalize()
        self.adjoint_analysis.Finalize()


    def _GetResponseFunctionUtility(self):
        return self.adjoint_analysis._GetSolver().response_function


    def _SynchronizeAdjointFromPrimal(self):
        if len(self.primal_model_part.Nodes) != len(self.adjoint_model_part.Nodes):
            raise RuntimeError("_SynchronizeAdjointFromPrimal: Model parts have a different number of nodes!")

        for primal_node, adjoint_node in zip(self.primal_model_part.Nodes, self.adjoint_model_part.Nodes):
            adjoint_node.X0 = primal_node.X0
            adjoint_node.Y0 = primal_node.Y0
            adjoint_node.Z0 = primal_node.Z0
            adjoint_node.X = primal_node.X
            adjoint_node.Y = primal_node.Y
            adjoint_node.Z = primal_node.Z

# ==============================================================================
class AdjointLinearStrainEnergyResponse(AdjointResponseFunction):
    def __init__(self, identifier, project_parameters, model):
        super(AdjointLinearStrainEnergyResponse, self).__init__(identifier, project_parameters, model)

    def CalculateValue(self):
        startTime = timer.time()
        #The linear strain energy response needs the primal model part to calculate the response value!
        value = self._GetResponseFunctionUtility().CalculateValue(self.primal_model_part)
        Logger.PrintInfo("> Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

        self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE] = value
