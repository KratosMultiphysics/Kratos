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

            # Dirichlet conditions: change variables
            for i in range(0,primal_parameters["processes"]["constraints_process_list"].size()):
                process = adjoint_parameters["processes"]["constraints_process_list"][i]
                variable_name = process["Parameters"]["variable_name"].GetString()
                process["Parameters"]["variable_name"].SetString("ADJOINT_"+variable_name)

            # Neumann conditions - do not modify to read the same load values as in primal:

            # Output process:
            # TODO how to add the output process? How find out about the variables?
            if adjoint_parameters.Has("output_configuration"):
                Logger.PrintInfo("> Output process is removed for adjoint analysis. To enable it define adjoint_parameters yourself.")
                adjoint_parameters.RemoveValue("output_configuration")

            # sensitivity settings
            adjoint_parameters["solver_settings"].AddValue("sensitivity_settings", self.response_settings["sensitivity_settings"])

            # response settings
            adjoint_parameters["solver_settings"].AddValue("response_function_settings", self.response_settings)

        else: # adjoint parameters file is explicitely given - do not change it.
            with open(self.response_settings["adjoint_settings"].GetString(),'r') as parameter_file:
                adjoint_parameters = Parameters( parameter_file.read() )

        return adjoint_parameters


# ==============================================================================
class AdjointBeamNormalStressResponseFunction(ResponseFunctionBase):
    """TODO MFusseder add description!
    """
    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings

        # Create the primal solver
        with open(self.response_settings["primal_settings"].GetString(),'r') as parameter_file:
            primal_parameters = Parameters( parameter_file.read() )

        self.primal_model_part = _GetModelPart(model, primal_parameters["solver_settings"])

        self.primal_analysis = StructuralMechanicsAnalysis(model, primal_parameters)

        # automatic creation is not provided for this response
        with open(self.response_settings["adjoint_settings"].GetString(),'r') as parameter_file:
                adjoint_parameters = Parameters( parameter_file.read() )

        self.__SetUpAdjointSubAnalysis()

        # adjustment of influence function. Needed to create adjoint fields.
        self.add_particular_solution = False
        if adjoint_parameters["solver_settings"]["response_function_settings"].Has("add_particular_solution"):
            self.add_particular_solution = adjoint_parameters["solver_settings"]["response_function_settings"]["add_particular_solution"].GetBool()
            # Set now the the bool to False in order to prevent wrong behaviour in dummy response
            adjoint_parameters["solver_settings"]["response_function_settings"]["add_particular_solution"].SetBool(False)

        # use here c++ local stress response as meaningful dummy
        adjoint_parameters["solver_settings"]["response_function_settings"]["response_type"].SetString("adjoint_local_stress")
        adjoint_parameters["solver_settings"]["response_function_settings"]["stress_type"].SetString("FX")
        adjoint_model = KratosMultiphysics.Model()
        self.adjoint_model_part = _GetModelPart(adjoint_model, adjoint_parameters["solver_settings"])
        self.adjoint_analysis = StructuralMechanicsAnalysis(adjoint_model, adjoint_parameters)
        self.sensitivity_settings = adjoint_parameters["solver_settings"]["sensitivity_settings"].Clone()
        #self.sensitivity_settings = self.response_settings["sensitivity_settings"]
        if adjoint_parameters["solver_settings"]["sensitivity_settings"].Has("element_cross_section_sensitivity_variables"):
            adjoint_parameters["solver_settings"]["sensitivity_settings"].RemoveValue("element_cross_section_sensitivity_variables")

        self.primal_state_variables = [KratosMultiphysics.DISPLACEMENT]
        self.adjoint_state_variables = [StructuralMechanicsApplication.ADJOINT_DISPLACEMENT]
        if primal_parameters["solver_settings"].Has("rotation_dofs"):
            if primal_parameters["solver_settings"]["rotation_dofs"].GetBool():
                self.primal_state_variables.append(KratosMultiphysics.ROTATION)
                self.adjoint_state_variables.append(StructuralMechanicsApplication.ADJOINT_ROTATION)

        self.traced_element_id = adjoint_parameters["solver_settings"]["response_function_settings"]["traced_element_id"].GetInt()
        # response_cross_section is the cross section of the traced element
        self.response_cross_section = None
        # this are cross sections of some sub model parts. Can be set by the function 'SetCrossSections'
        self.cross_sections = []

        # here the position within the cross section is defined where the the stress is traced
        self.stress_position_y = adjoint_parameters["solver_settings"]["response_function_settings"]["yp"].GetDouble()
        self.stress_position_z = adjoint_parameters["solver_settings"]["response_function_settings"]["zp"].GetDouble()
        y_stress_position = adjoint_parameters["solver_settings"]["response_function_settings"]["y_stress_position"].GetString()
        z_stress_position = adjoint_parameters["solver_settings"]["response_function_settings"]["z_stress_position"].GetString()

        if y_stress_position == "fixed_coordinate":
            self.adaptive_y_coord = False
        elif y_stress_position == "adaptive_coordinate": # here we assume a symmetric cross-section in local y-direction
            self.adaptive_y_coord = True
            if -0.5 <= self.stress_position_y <= 0.5: # coordinate has to be between -0.5 and 0.5 times the cross sectional dimension in y-direction
                pass
            else:
                raise RuntimeError('stress_position_y has to be between -0.5 and 0.5')
        else:
            raise RuntimeError('Unknown setting for y_stress_position!')

        if z_stress_position == "fixed_coordinate":
            self.adaptive_z_coord = False
        elif z_stress_position == "adaptive_coordinate": # here we assume a symmetric cross-section in local z-direction
            self.adaptive_z_coord = True
            if -0.5 <= self.stress_position_z <= 0.5: # coordinate has to be between -0.5 and 0.5 times the cross sectional dimension in z-direction
                pass
            else:
                raise RuntimeError('stress_position_z has to be between -0.5 and 0.5')
        else:
            raise RuntimeError('Unknown setting for z_stress_position!')

    # --------------------------------------------------------------------------
    # MFusseder TODO find a better solution to make cross sections available
    def SetCrossSections(self, cross_sections):
        self.cross_sections = cross_sections

    # --------------------------------------------------------------------------
    def SetResponseCrossSection(self, response_cross_section):
        self.response_cross_section = response_cross_section

    # --------------------------------------------------------------------------
    def Initialize(self):
        self._Check()
        self.primal_analysis.Initialize()
        self.adjoint_analysis.Initialize() # Here the kratos response is created. what to here in adjoint solver?
        self.adjoint_analysis_fx.Initialize()
        self.adjoint_analysis_my.Initialize()
        self.adjoint_analysis_mz.Initialize()

    # --------------------------------------------------------------------------
    def InitializeSolutionStep(self):
        # Run the primal analysis.
        # TODO if primal_analysis.status==solved: return
        Logger.PrintInfo("\n> Starting primal analysis for response:", self.identifier)
        startTime = timer.time()
        self._RunSolutionLoopPrimalAnalysis()
        Logger.PrintInfo("> Time needed for solving the primal analysis = ",round(timer.time() - startTime,2),"s")

    # --------------------------------------------------------------------------
    def CalculateValue(self):
        startTime = timer.time()
        FX = self.adjoint_analysis_fx._GetSolver().response_function.CalculateValue(self.primal_model_part)
        MY = self.adjoint_analysis_my._GetSolver().response_function.CalculateValue(self.primal_model_part)
        MZ = self.adjoint_analysis_mz._GetSolver().response_function.CalculateValue(self.primal_model_part)
        yp, zp = self._GetStressPositionWithinCrossSection()
        A = self.primal_model_part.GetElement(self.traced_element_id).Properties[StructuralMechanicsApplication.CROSS_AREA]
        I22 = self.primal_model_part.GetElement(self.traced_element_id).Properties[StructuralMechanicsApplication.I22]
        I33 = self.primal_model_part.GetElement(self.traced_element_id).Properties[StructuralMechanicsApplication.I33]

        value = FX / A + MY / I22 * zp + MZ / I33 * yp # TODO evaluate signs of stress resultants
        Logger.PrintInfo("> Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")

        self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE] = value

    # --------------------------------------------------------------------------
    def CalculateGradient(self):
        # synchronize the modelparts
        #self._SynchronizeAdjointFromPrimal() TODO is this needed?
        startTime = timer.time()
        Logger.PrintInfo("\n> Starting adjoint analysis for response:", self.identifier)

        # compute adjoint subproblems (solve adjoint analysis for fx, my and mz as responses)
        self._RunSolutionLoopAdjointSubProblems()

        # compute sensitivities w.r.t. specific variables which are not known by kratos elements (e.g. SECTION_HEIGTH_SENSITIVITY)
        if self.sensitivity_settings.Has("element_cross_section_sensitivity_variables"):
            variables, _ = GenerateVariableListFromInput(self.sensitivity_settings["element_cross_section_sensitivity_variables"], False)
            ComputeSpecificCrossSectionSensitivities(variables, self.cross_sections, self.adjoint_model_part_fx)
            ComputeSpecificCrossSectionSensitivities(variables, self.cross_sections, self.adjoint_model_part_my)
            ComputeSpecificCrossSectionSensitivities(variables, self.cross_sections, self.adjoint_model_part_mz)

        # here the results of the adjoint analysisn are combinded
        self._RunSolutionLoopMainProblem()
        self._UpdateSensitivities()

        Logger.PrintInfo("> Time needed for solving the adjoint analysis = ", round(timer.time() - startTime,2),"s")

    # --------------------------------------------------------------------------
    def GetValue(self):
        return self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE]

    # --------------------------------------------------------------------------
    def FinalizeSolutionStep(self):
        # add particular solution in order to compute and visualize the adjoint fields within 'OutputSolutionStep'.
        if self.add_particular_solution:
            self._AddParticularSolutionToStressInfluenceFunction()
        self.adjoint_analysis.OutputSolutionStep()

    # --------------------------------------------------------------------------
    def Finalize(self):
        self.primal_analysis.Finalize()
        self.adjoint_analysis.Finalize()
        self.adjoint_analysis_fx.Finalize()
        self.adjoint_analysis_my.Finalize()
        self.adjoint_analysis_mz.Finalize()

    """
    **************************************************************************************************************
    PROTECTED MEMBER FUNCTIONS
    **************************************************************************************************************
    """
    # --------------------------------------------------------------------------
    def _Check(self):
        if (self.response_cross_section is None) and ((self.adaptive_y_coord == True) or (self.adaptive_z_coord == True)):
            raise RuntimeError('adaptive_coordinate is only possible with a given response_cross_section!')

    # --------------------------------------------------------------------------
    def _RunSolutionLoopPrimalAnalysis(self):
        if not self.primal_analysis.time < self.primal_analysis.end_time:
            self.primal_analysis.end_time += 1
        self.primal_analysis.RunSolutionLoop()

    # --------------------------------------------------------------------------
    def _RunSolutionLoopAdjointSubProblems(self):
        # compute here all contributions to sensitivities and influence functions
        self.adjoint_analysis_fx.RunSolutionLoop()
        self.adjoint_analysis_my.RunSolutionLoop()
        self.adjoint_analysis_mz.RunSolutionLoop()

    # --------------------------------------------------------------------------
    def _RunSolutionLoopMainProblem(self):
        # TODO are all this steps of analysis stage necessary?
        self.adjoint_analysis.time = self.adjoint_analysis._GetSolver().AdvanceInTime(self.adjoint_analysis.time)
        self.adjoint_analysis.InitializeSolutionStep()
        self.adjoint_analysis._GetSolver().Predict()

        A = self.primal_model_part.GetElement(self.traced_element_id).Properties[StructuralMechanicsApplication.CROSS_AREA]
        I22 = self.primal_model_part.GetElement(self.traced_element_id).Properties[StructuralMechanicsApplication.I22]
        I33 = self.primal_model_part.GetElement(self.traced_element_id).Properties[StructuralMechanicsApplication.I33]
        yp, zp = self._GetStressPositionWithinCrossSection()

        # solve here adjoint problem by superposition TODO improve this
        for variable in self.adjoint_state_variables:
            for node, node_fx, node_my, node_mz in zip(self.adjoint_model_part.Nodes, self.adjoint_model_part_fx.Nodes, self.adjoint_model_part_my.Nodes, self.adjoint_model_part_mz.Nodes):
                adj_disp_fx = node_fx.GetSolutionStepValue(variable)
                adj_disp_my = node_my.GetSolutionStepValue(variable)
                adj_disp_mz = node_mz.GetSolutionStepValue(variable)
                # TODO evaluate signs of stress resultants
                adjoint_quantity = adj_disp_fx / A + adj_disp_my / I22 * zp + adj_disp_mz / I33 * yp
                node.SetSolutionStepValue(variable, adjoint_quantity)

        # it is not possible to use the sensitivity builder of the adjoint analysis even the adjoint variables are correct now.
        # The c++ dummy response functions delivers not the needed partial derivatives w.r.t. design variable.
        # This is the reason why the function '_UpdateSensitivities' is needed.
        self.adjoint_analysis.FinalizeSolutionStep()

    # --------------------------------------------------------------------------
    def _UpdateSensitivities(self):
        nodal_variables, _ = GenerateVariableListFromInput(self.sensitivity_settings["nodal_solution_step_sensitivity_variables"], False)
        element_variables, element_sensitivity_variables = GenerateVariableListFromInput(self.sensitivity_settings["element_data_value_sensitivity_variables"], True)
        condition_variables, condition_sensitivity_variables = GenerateVariableListFromInput(self.sensitivity_settings["condition_data_value_sensitivity_variables"], True)
        sen_mp = GetSensitivityModelPart(self.sensitivity_settings["sensitivity_model_part_name"].GetString(), self.adjoint_model_part)
        sen_mp_fx = GetSensitivityModelPart(self.sensitivity_settings["sensitivity_model_part_name"].GetString(), self.adjoint_model_part_fx)
        sen_mp_my = GetSensitivityModelPart(self.sensitivity_settings["sensitivity_model_part_name"].GetString(), self.adjoint_model_part_my)
        sen_mp_mz = GetSensitivityModelPart(self.sensitivity_settings["sensitivity_model_part_name"].GetString(), self.adjoint_model_part_mz)
        A = self.primal_model_part.GetElement(self.traced_element_id).Properties[StructuralMechanicsApplication.CROSS_AREA]
        I22 = self.primal_model_part.GetElement(self.traced_element_id).Properties[StructuralMechanicsApplication.I22]
        I33 = self.primal_model_part.GetElement(self.traced_element_id).Properties[StructuralMechanicsApplication.I33]
        yp, zp = self._GetStressPositionWithinCrossSection()
        # response superposition for nodal design variables
        for var in nodal_variables:
            if var.Name() == 'SHAPE_SENSITIVITY':
                for node, node_fx, node_my, node_mz in zip(sen_mp.Nodes, sen_mp_fx.Nodes, sen_mp_my.Nodes, sen_mp_mz.Nodes):
                    sen_fx = node_fx.GetSolutionStepValue(var)
                    sen_my = node_my.GetSolutionStepValue(var)
                    sen_mz = node_mz.GetSolutionStepValue(var)
                    sensitivity =  sen_fx / A + sen_my / I22 * zp + sen_mz / I33 * yp
                    node.SetSolutionStepValue(var, sensitivity)
            else:
                raise RuntimeError(var.Name(), ' not available!')

        # response superposition for elemental design variables
        for var, sen_var in zip(element_variables, element_sensitivity_variables):
            for elem, elem_fx, elem_my, elem_mz in zip(sen_mp.Elements, sen_mp_fx.Elements, sen_mp_my.Elements, sen_mp_mz.Elements):
                sen_fx = elem_fx.GetValue(sen_var)
                sen_my = elem_my.GetValue(sen_var)
                sen_mz = elem_mz.GetValue(sen_var)
                sensitivity =  sen_fx / A + sen_my / I22 * zp + sen_mz / I33 * yp

                # add additional contribution from partial derivative w.r.t. design variable
                # In this cases it is not considered that 'zp' or 'yp' can be dependent on the the
                # following cross sectional properties
                if (elem.Id == self.traced_element_id) and (var.Name() == 'CROSS_AREA'):
                    FX = self.adjoint_analysis_fx._GetSolver().response_function.CalculateValue(self.primal_model_part)
                    sensitivity -= FX / A**2
                if (elem.Id == self.traced_element_id) and (var.Name() == 'I22'):
                    MY = self.adjoint_analysis_my._GetSolver().response_function.CalculateValue(self.primal_model_part)
                    sensitivity -= MY / I22**2 * zp
                if (elem.Id == self.traced_element_id) and (var.Name() == 'I33'):
                    MZ = self.adjoint_analysis_mz._GetSolver().response_function.CalculateValue(self.primal_model_part)
                    sensitivity -= MZ / I33**2 * yp

                elem.SetValue(sen_var, sensitivity)

        # response superposition for conditional design variables
        for var, sen_var in zip(condition_variables, condition_sensitivity_variables):
            if var.Name() == 'POINT_LOAD':
                for cond, cond_fx, cond_my, cond_mz in zip(sen_mp.Conditions, sen_mp_fx.Conditions, sen_mp_my.Conditions, sen_mp_mz.Conditions):
                    sen_fx = cond_fx.GetValue(sen_var)
                    sen_my = cond_my.GetValue(sen_var)
                    sen_mz = cond_mz.GetValue(sen_var)
                    sensitivity =  sen_fx / A + sen_my / I22 * zp + sen_mz / I33 * yp
                    cond.SetValue(sen_var, sensitivity)
            else:
                raise RuntimeError(sen_var.Name(), ' not available!')

        if self.sensitivity_settings.Has("element_cross_section_sensitivity_variables"):
            elem_cs_variables, _ = GenerateVariableListFromInput(self.sensitivity_settings["element_cross_section_sensitivity_variables"], False)
            for var in elem_cs_variables:
                for elem, elem_fx, elem_my, elem_mz in zip(sen_mp.Elements, sen_mp_fx.Elements, sen_mp_my.Elements, sen_mp_mz.Elements):
                    sen_fx = elem_fx.GetValue(var)
                    sen_my = elem_my.GetValue(var)
                    sen_mz = elem_mz.GetValue(var)
                    sensitivity =  sen_fx / A + sen_my / I22 * zp + sen_mz / I33 * yp

                    if (elem.Id == self.traced_element_id) and ((var.Name() == 'SECTION_HEIGTH_SENSITIVITY') or (var.Name() == 'SECTION_WIDTH_SENSITIVITY')):
                        if self.response_cross_section is None:
                            raise RuntimeError('It is not possible to compute ' + var.Name() + ' without a defined cross section!')
                        FX = self.adjoint_analysis_fx._GetSolver().response_function.CalculateValue(self.primal_model_part)
                        MY = self.adjoint_analysis_my._GetSolver().response_function.CalculateValue(self.primal_model_part)
                        MZ = self.adjoint_analysis_mz._GetSolver().response_function.CalculateValue(self.primal_model_part)
                        dA_dX = self.response_cross_section.ComputeCrossAreaDerivative(var)
                        dI22_dX = self.response_cross_section.ComputeI22Derivative(var)
                        dI33_dX = self.response_cross_section.ComputeI33Derivative(var)
                        sensitivity -= (FX / A**2 * dA_dX +
                                        MY / I22**2 * zp * dI22_dX +
                                        MZ / I33**2 * yp * dI33_dX)
                        if self.adaptive_y_coord and (var.Name() == 'SECTION_WIDTH_SENSITIVITY'):
                            sensitivity += MZ / I33 * self.stress_position_y
                        if self.adaptive_z_coord and (var.Name() == 'SECTION_HEIGTH_SENSITIVITY'):
                            sensitivity += MY / I22 * self.stress_position_z
                    elem.SetValue(var, sensitivity)

    # --------------------------------------------------------------------------
    def _GetStressPositionWithinCrossSection(self):
        yp = self.stress_position_y
        zp = self.stress_position_z
        if self.adaptive_y_coord:
            dim_y = self.response_cross_section.GetCharacteristicDimensionY()
            yp *= dim_y
        if self.adaptive_z_coord:
            dim_z = self.response_cross_section.GetCharacteristicDimensionZ()
            zp *= dim_z
        return yp, zp

    # --------------------------------------------------------------------------
    def _AddParticularSolutionToStressInfluenceFunction(self):
        # in this function the negative unit pre-deformations of the the stress resultant influence functions
        # are scaled according to the stress formula and summed up leading to the particular solution of the
        # normal stress influence function. Needed to compute the adjoint fields within the method of
        # generalized influence functions.
        yp, zp = self._GetStressPositionWithinCrossSection()
        A = self.primal_model_part.GetElement(self.traced_element_id).Properties[StructuralMechanicsApplication.CROSS_AREA]
        I22 = self.primal_model_part.GetElement(self.traced_element_id).Properties[StructuralMechanicsApplication.I22]
        I33 = self.primal_model_part.GetElement(self.traced_element_id).Properties[StructuralMechanicsApplication.I33]

        elem_fx = self.adjoint_model_part_fx.GetElement(self.traced_element_id)
        elem_my = self.adjoint_model_part_my.GetElement(self.traced_element_id)
        elem_mz = self.adjoint_model_part_mz.GetElement(self.traced_element_id)
        elem_stress = self.adjoint_model_part.GetElement(self.traced_element_id)

        if  (elem_fx.Has(StructuralMechanicsApplication.ADJOINT_PARTICULAR_DISPLACEMENT) == False or
                elem_my.Has(StructuralMechanicsApplication.ADJOINT_PARTICULAR_DISPLACEMENT) == False or
                elem_mz.Has(StructuralMechanicsApplication.ADJOINT_PARTICULAR_DISPLACEMENT) == False):
            raise RuntimeError('Some or all adjoint sub analysis have no particular solution!')

        part_sol_fx = elem_fx.GetValue(StructuralMechanicsApplication.ADJOINT_PARTICULAR_DISPLACEMENT)
        part_sol_my = elem_my.GetValue(StructuralMechanicsApplication.ADJOINT_PARTICULAR_DISPLACEMENT)
        part_sol_mz = elem_mz.GetValue(StructuralMechanicsApplication.ADJOINT_PARTICULAR_DISPLACEMENT)
        complete_part_sol = (part_sol_fx / A +
                             part_sol_my / I22 * zp +
                             part_sol_mz / I33 * yp )
        elem_stress.SetValue(StructuralMechanicsApplication.ADJOINT_PARTICULAR_DISPLACEMENT, complete_part_sol)


    """
    **************************************************************************************************************
    PRIVATE MEMBER FUNCTIONS
    **************************************************************************************************************
    """
    # --------------------------------------------------------------------------
    def __SetUpAdjointSubAnalysis(self):
        with open(self.response_settings["adjoint_settings"].GetString(),'r') as parameter_file:
            adjoint_parameters_fx = Parameters( parameter_file.read() )
        with open(self.response_settings["adjoint_settings"].GetString(),'r') as parameter_file:
            adjoint_parameters_my = Parameters( parameter_file.read() )
        with open(self.response_settings["adjoint_settings"].GetString(),'r') as parameter_file:
            adjoint_parameters_mz = Parameters( parameter_file.read() )

        # adjoint analysis for FX
        adjoint_parameters_fx["solver_settings"]["response_function_settings"]["response_type"].SetString("adjoint_local_stress")
        adjoint_parameters_fx["solver_settings"]["response_function_settings"]["stress_type"].SetString("FX")
        adjoint_parameters_fx.RemoveValue("output_processes")
        if adjoint_parameters_fx["solver_settings"]["sensitivity_settings"].Has("element_cross_section_sensitivity_variables"):
            adjoint_parameters_fx["solver_settings"]["sensitivity_settings"].RemoveValue("element_cross_section_sensitivity_variables")
        adjoint_model_fx = KratosMultiphysics.Model()
        self.adjoint_model_part_fx = _GetModelPart(adjoint_model_fx, adjoint_parameters_fx["solver_settings"])
        self.adjoint_analysis_fx = StructuralMechanicsAnalysis(adjoint_model_fx, adjoint_parameters_fx)
        # adjoint analysis for MY
        adjoint_parameters_my["solver_settings"]["response_function_settings"]["response_type"].SetString("adjoint_local_stress")
        adjoint_parameters_my["solver_settings"]["response_function_settings"]["stress_type"].SetString("MY")
        adjoint_parameters_my.RemoveValue("output_processes")
        if adjoint_parameters_my["solver_settings"]["sensitivity_settings"].Has("element_cross_section_sensitivity_variables"):
            adjoint_parameters_my["solver_settings"]["sensitivity_settings"].RemoveValue("element_cross_section_sensitivity_variables")
        adjoint_model_my = KratosMultiphysics.Model()
        self.adjoint_model_part_my = _GetModelPart(adjoint_model_my, adjoint_parameters_my["solver_settings"])
        self.adjoint_analysis_my = StructuralMechanicsAnalysis(adjoint_model_my, adjoint_parameters_my)
        # adjoint analysis for MZ
        adjoint_parameters_mz["solver_settings"]["response_function_settings"]["response_type"].SetString("adjoint_local_stress")
        adjoint_parameters_mz["solver_settings"]["response_function_settings"]["stress_type"].SetString("MZ")
        adjoint_parameters_mz.RemoveValue("output_processes")
        if adjoint_parameters_mz["solver_settings"]["sensitivity_settings"].Has("element_cross_section_sensitivity_variables"):
            adjoint_parameters_mz["solver_settings"]["sensitivity_settings"].RemoveValue("element_cross_section_sensitivity_variables")
        adjoint_model_mz = KratosMultiphysics.Model()
        self.adjoint_model_part_mz = _GetModelPart(adjoint_model_mz, adjoint_parameters_mz["solver_settings"])
        self.adjoint_analysis_mz = StructuralMechanicsAnalysis(adjoint_model_mz, adjoint_parameters_mz)

# *************************************************************************************************************
# ADDITIONAL FUNCTIONS
# *************************************************************************************************************
# --------------------------------------------------------------------------
def GenerateVariableListFromInput(parameter, create_sensitivity_variable_list):
    if not parameter.IsArray():
        raise Exception("Error: Variable list is unreadable")
    variable_list = []
    sensitivity_variable_list = []
    for i in range(0, parameter.size()):
        variable_list.append(KratosMultiphysics.KratosGlobals.GetVariable( parameter[i].GetString() ))
        if create_sensitivity_variable_list:
            sensitivity_variable_list.append(KratosMultiphysics.KratosGlobals.GetVariable( parameter[i].GetString() + "_SENSITIVITY"))
    return variable_list, sensitivity_variable_list

# --------------------------------------------------------------------------
def GetSensitivityModelPart(model_part_name, root_mp):
    if root_mp.HasSubModelPart(model_part_name):
        model_part = root_mp.GetSubModelPart(model_part_name)
        return model_part
    elif root_mp.Name == model_part_name: # TODO is this necessary?
        return root_mp
    else:
        raise RuntimeError('Given model part ' + model_part_name + ' is not available!')

# --------------------------------------------------------------------------
def ComputeSpecificCrossSectionSensitivities(variables_list, cross_sections, model_part):
    for var_i in variables_list:
        for sec_i in cross_sections:
            sec_i.ComputeSpecificCrossSectionSensitivities(var_i, model_part)

