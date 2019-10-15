"""This module contains the available structural response functions and their base class"""
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics import Parameters, Logger
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

        self.primal_analysis = StructuralMechanicsAnalysis(model, ProjectParametersPrimal)
        self.primal_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)

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
            gradient[node.Id] = node.GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY)
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

        self.primal_analysis = StructuralMechanicsAnalysis(model, ProjectParametersPrimal)
        self.primal_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)

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

        model_part_name = response_settings["model_part_name"].GetString()
        input_type = response_settings["model_import_settings"]["input_type"].GetString()
        if input_type == "mdpa":
            self.model_part = self.model.CreateModelPart(model_part_name, 2)
            domain_size = response_settings["domain_size"].GetInt()
            if domain_size not in [2, 3]:
                raise Exception("MassResponseFunction: Invalid 'domain_size': {}".format(domain_size))
            self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)
            self.model_part_needs_to_be_imported = True
        elif input_type == "use_input_model_part":
            self.model_part = self.model.GetModelPart(model_part_name)
        else:
            raise Exception("Other model part input options are not yet implemented.")

        self.response_function_utility = StructuralMechanicsApplication.MassResponseFunctionUtility(self.model_part, response_settings)

        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)

    def Initialize(self):
        if self.model_part_needs_to_be_imported:
            # import model part
            model_part_io = KratosMultiphysics.ModelPartIO(self.response_settings["model_import_settings"]["input_filename"].GetString())
            model_part_io.ReadModelPart(self.model_part)

        # Add constitutive laws and material properties from json file to model parts.
        material_settings = KratosMultiphysics.Parameters("""{"Parameters": {} }""")
        materials_file_name = self.response_settings["material_import_settings"]["materials_filename"]
        material_settings["Parameters"].AddValue("materials_filename", materials_file_name)
        KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)
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
            gradient[node.Id] = node.GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY)
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
    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings

        # Create the primal solver
        with open(self.response_settings["primal_settings"].GetString(),'r') as parameter_file:
            primal_parameters = Parameters( parameter_file.read() )

        self.primal_model_part = _GetModelPart(model, primal_parameters["solver_settings"])

        self.primal_analysis = StructuralMechanicsAnalysis(model, primal_parameters)

        # Create the adjoint solver
        adjoint_parameters = self._GetAdjointParameters()
        adjoint_model = KratosMultiphysics.Model()
        self.adjoint_model_part = _GetModelPart(adjoint_model, adjoint_parameters["solver_settings"])

        # TODO find out why it is not possible to use the same model_part
        self.adjoint_analysis = StructuralMechanicsAnalysis(adjoint_model, adjoint_parameters)

        self.primal_state_variables = [KratosMultiphysics.DISPLACEMENT]
        if primal_parameters["solver_settings"].Has("rotation_dofs"):
            if primal_parameters["solver_settings"]["rotation_dofs"].GetBool():
                self.primal_state_variables.append(KratosMultiphysics.ROTATION)

    def Initialize(self):
        self.primal_analysis.Initialize()
        self.adjoint_analysis.Initialize()

    def InitializeSolutionStep(self):
        # Run the primal analysis.
        # TODO if primal_analysis.status==solved: return
        Logger.PrintInfo("\n> Starting primal analysis for response:", self.identifier)
        startTime = timer.time()
        if not self.primal_analysis.time < self.primal_analysis.end_time:
            self.primal_analysis.end_time += 1
        self.primal_analysis.RunSolutionLoop()
        Logger.PrintInfo("> Time needed for solving the primal analysis = ",round(timer.time() - startTime,2),"s")

    def CalculateValue(self):
        startTime = timer.time()
        value = self._GetResponseFunctionUtility().CalculateValue(self.primal_model_part)
        Logger.PrintInfo("> Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

        self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE] = value

    def CalculateGradient(self):
        # synchronize the modelparts
        self._SynchronizeAdjointFromPrimal()
        startTime = timer.time()
        Logger.PrintInfo("\n> Starting adjoint analysis for response:", self.identifier)
        if not self.adjoint_analysis.time < self.adjoint_analysis.end_time:
            self.adjoint_analysis.end_time += 1
        self.adjoint_analysis.RunSolutionLoop()
        Logger.PrintInfo("> Time needed for solving the adjoint analysis = ",round(timer.time() - startTime,2),"s")

    def GetValue(self):
        return self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE]

    def GetShapeGradient(self):
        gradient = {}
        for node in self.adjoint_model_part.Nodes:
            gradient[node.Id] = node.GetSolutionStepValue(KratosMultiphysics.SHAPE_SENSITIVITY)
        return gradient

    def Finalize(self):
        self.primal_analysis.Finalize()
        self.adjoint_analysis.Finalize()

    def _GetResponseFunctionUtility(self):
        return self.adjoint_analysis._GetSolver().response_function

    def _SynchronizeAdjointFromPrimal(self):
        Logger.PrintInfo("\n> Synchronize primal and adjoint modelpart for response:", self.identifier)

        if len(self.primal_model_part.Nodes) != len(self.adjoint_model_part.Nodes):
            raise RuntimeError("_SynchronizeAdjointFromPrimal: Model parts have a different number of nodes!")

        # TODO this should happen automatically
        for primal_node, adjoint_node in zip(self.primal_model_part.Nodes, self.adjoint_model_part.Nodes):
            adjoint_node.X0 = primal_node.X0
            adjoint_node.Y0 = primal_node.Y0
            adjoint_node.Z0 = primal_node.Z0
            adjoint_node.X = primal_node.X
            adjoint_node.Y = primal_node.Y
            adjoint_node.Z = primal_node.Z

        # Put primal solution on adjoint model - for "auto" setting, else it has to be done by the user e.g. using hdf5 process
        if self.response_settings["adjoint_settings"].GetString() == "auto":
            Logger.PrintInfo("> Transfer primal state to adjoint model part.")
            variable_utils = KratosMultiphysics.VariableUtils()
            for variable in self.primal_state_variables:
                variable_utils.CopyModelPartNodalVar(variable, self.primal_model_part, self.adjoint_model_part, 0)


    def _GetAdjointParameters(self):

        adjoint_settings = self.response_settings["adjoint_settings"].GetString()

        if adjoint_settings == "auto":
            Logger.PrintInfo("\n> Automatic set up adjoint parameters for response:", self.identifier)

            with open(self.response_settings["primal_settings"].GetString(),'r') as parameter_file:
                primal_parameters = Parameters( parameter_file.read() )

            # check that HDF5 process is not there
            if primal_parameters["processes"].Has("list_other_processes"):
                for i in range(0,primal_parameters["processes"]["list_other_processes"].size()):
                    process = primal_parameters["processes"]["list_other_processes"][i]
                    raise Exception("Auto setup of adjoint parameters does not support {} in list_other_processes".format(process["python_module"].GetString()))

            # clone primal settings as base for adjoint
            adjoint_parameters = primal_parameters.Clone()

            # analysis settings
            solver_settings = adjoint_parameters["solver_settings"]
            primal_solver_type = solver_settings["solver_type"].GetString()
            if primal_solver_type != "static":
                raise Exception("Auto setup of adjoint parameters does not support {} solver_type. Only available for 'static'".format(primal_solver_type))
            solver_settings["solver_type"].SetString("adjoint_"+primal_solver_type)

            if not solver_settings.Has("compute_reactions"):
                solver_settings.AddEmptyValue("compute_reactions")
            solver_settings["compute_reactions"].SetBool(False)

            if not solver_settings.Has("move_mesh_flag"):
                solver_settings.AddEmptyValue("move_mesh_flag")
            solver_settings["move_mesh_flag"].SetBool(False)

            if solver_settings.Has("scheme_settings"):
                depr_msg = '\nDEPRECATION-WARNING: "scheme_settings" is deprecated, please remove it from your json parameters.\n'
                KratosMultiphysics.Logger.PrintWarning(__name__, depr_msg)
                solver_settings.RemoveValue("scheme_settings")

            if solver_settings["model_import_settings"]["input_type"].GetString() == "use_input_model_part":
                solver_settings["model_import_settings"]["input_type"].SetString("mdpa")
                solver_settings["model_import_settings"].AddEmptyValue("input_filename")
                model_part_name = solver_settings["model_part_name"].GetString()
                solver_settings["model_import_settings"]["input_filename"].SetString(model_part_name)

            # Dirichlet conditions: change variables
            for i in range(0,primal_parameters["processes"]["constraints_process_list"].size()):
                process = adjoint_parameters["processes"]["constraints_process_list"][i]
                variable_name = process["Parameters"]["variable_name"].GetString()
                process["Parameters"]["variable_name"].SetString("ADJOINT_"+variable_name)

            # Neumann conditions - do not modify to read the same load values as in primal:

            # Output process:
            # TODO how to add the output process? How find out about the variables?
            if adjoint_parameters.Has("output_processes"):
                Logger.PrintInfo("> Output process is removed for adjoint analysis. To enable it define adjoint_parameters yourself.")
                adjoint_parameters.RemoveValue("output_processes")

            # sensitivity settings
            adjoint_parameters["solver_settings"].AddValue("sensitivity_settings", self.response_settings["sensitivity_settings"])

            # response settings
            adjoint_parameters["solver_settings"].AddValue("response_function_settings", self.response_settings)

        else: # adjoint parameters file is explicitely given - do not change it.
            with open(self.response_settings["adjoint_settings"].GetString(),'r') as parameter_file:
                adjoint_parameters = Parameters( parameter_file.read() )

        return adjoint_parameters
