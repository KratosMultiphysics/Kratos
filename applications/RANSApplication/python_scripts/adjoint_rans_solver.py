import numpy as np

from pathlib import Path

# Importing the Kratos Library
import KratosMultiphysics as Kratos

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.RANSApplication as KratosRANS

from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
from KratosMultiphysics.RANSApplication.formulations.utilities import CalculateNormalsOnConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import GetKratosObjectPrototype
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateBlockBuilderAndSolver
from KratosMultiphysics.RANSApplication.formulations.utilities import GetTimeDerivativeVariablesRecursively
from KratosMultiphysics.RANSApplication.formulations.utilities import ExecutionScope
from KratosMultiphysics.RANSApplication.formulations.utilities import SolveProblem

from KratosMultiphysics.FluidDynamicsApplication.adjoint_stabilization_utilities import ComputeStabilizationCoefficient

try:
    from KratosMultiphysics.HDF5Application.single_mesh_temporal_input_process import Factory as InputFactory
    from KratosMultiphysics.HDF5Application.single_mesh_temporal_output_process import Factory as OutputFactory
except ImportError:
    InputFactory = None
    OutputFactory = None

# case specific imports
if (IsDistributedRun() and CheckIfApplicationsAvailable("TrilinosApplication")):
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos     # MPI solvers
elif (IsDistributedRun()):
    raise Exception("Distributed run requires TrilinosApplication")

# Import base class
from KratosMultiphysics.RANSApplication.coupled_rans_solver import CoupledRANSSolver

def CreateSolver(main_model_part, custom_settings):
    return AdjointRANSSolver(main_model_part, custom_settings)

class AdjointSolutionController:
    def __init__(self, main_model_part, parameters):
        self.main_model_part = main_model_part
        self.parameters = parameters
        self.check_step = -1
        self.inteval_steps = 0
        self.diffusion_coefficient = 0.0

    def AddTimeSchemeParameters(self, parameters):
        pass

    def InitializeSolutionStep(self):
        if (self.check_step == -1):
            self.check_step = self.main_model_part.ProcessInfo[Kratos.STEP]

        self.inteval_steps = self.main_model_part.ProcessInfo[Kratos.STEP] - self.check_step

    def FinalizeSolutionStep(self):
        pass

class DiffusionSolutionController(AdjointSolutionController):
    def __init__(self, main_model_part, parameters):
        super().__init__(main_model_part, parameters)

        default_parameters = Kratos.Parameters("""{
            "solution_control_method"                : "diffusion",
            "stabilization_coefficient"              : 0.0,
            "stabilization_residual_scaling_vector"  : [1.0],
            "stabilization_derivative_scaling_vector": [1.0],
            "elemental_stabilization_coefficient_variable_name": "PLEASE_SPECIFY_SCALAR_VARIABLE_NAME",
            "stabilization_computation_settings"     :{}
        }""")

        if self.parameters.Has("stabilization_coefficient"):
            if self.parameters["stabilization_coefficient"].IsString():
                default_parameters["stabilization_coefficient"].SetString("")

        self.parameters.ValidateAndAssignDefaults(default_parameters)

    def AddTimeSchemeParameters(self, parameters):
        if self.parameters.Has("stabilization_coefficient"):
            if self.parameters["stabilization_coefficient"].IsString():
                if self.parameters["stabilization_coefficient"].GetString() == "computed":
                    from KratosMultiphysics.RANSApplication.adjoint_rans_analysis import AdjointRANSAnalysis
                    stabilization_coefficient = ComputeStabilizationCoefficient(
                                                    AdjointRANSAnalysis,
                                                    self.parameters["stabilization_computation_settings"],
                                                    DiffusionSolutionController.__ExecuteAnalysis)
                else:
                    raise Exception("Only \"computed\" is supported as a string literal for stabilization_coefficient, or use double value. [ stabilization_coefficient = {:s} ].".format(self.parameters["stabilization_coefficient"].GetString()))
            elif self.parameters["stabilization_coefficient"].IsDouble():
                stabilization_coefficient = self.parameters["stabilization_coefficient"].GetDouble()
            else:
                raise Exception("Only doubel and \"computed\" is supported as values for stabilization_coefficient.")

        # now add only necessary settings
        parameters.AddValue("stabilization_residual_scaling_vector", self.parameters["stabilization_residual_scaling_vector"])
        parameters.AddValue("stabilization_derivative_scaling_vector", self.parameters["stabilization_derivative_scaling_vector"])
        parameters.AddEmptyValue("stabilization_coefficient")
        parameters.AddValue("elemental_stabilization_coefficient_variable_name", self.parameters["elemental_stabilization_coefficient_variable_name"])
        parameters["stabilization_coefficient"].SetDouble(stabilization_coefficient)

    @staticmethod
    def __ExecuteAnalysis(analysis_class_type, adjoint_parameters, stabilization_coefficient, solve_id):
        Kratos.Logger.PrintInfo("Stabilization Analysis", "Solving adjoint problem with {:e} stabilization coefficient.".format(stabilization_coefficient))
        scheme_settings = adjoint_parameters["solver_settings"]["adjoint_solution_controller_settings"]
        if not scheme_settings.Has("stabilization_coefficient"):
            scheme_settings.AddEmptyValue("stabilization_coefficient")
        scheme_settings["stabilization_coefficient"].SetDouble(stabilization_coefficient)

        with ExecutionScope("adjoint_stabilization_analysis/{:d}".format(solve_id)):
            execute_analysis = not Path("adjoint_stabilization_data.dat").is_file()
            if not execute_analysis:
                if Path("stabilization_coefficient.json").is_file():
                    with open("stabilization_coefficient.json", "r") as file_input:
                        kratos_parameters = Kratos.Parameters(file_input.read())
                    execute_analysis = kratos_parameters["solver_settings"]["adjoint_solution_controller_settings"]["stabilization_coefficient"].GetDouble() != stabilization_coefficient
                else:
                    execute_analysis = True

            if execute_analysis:
                _, simulation = SolveProblem(analysis_class_type, adjoint_parameters, "stabilization_coefficient")
                return simulation.time_series_data
            else:
                Kratos.Logger.PrintInfo("Stabilization Analysis", "Found existing solution at {:s}".format(str(Path(".").absolute())))
                return np.loadtxt("adjoint_stabilization_data.dat")

class ResetToZeroAdjointSolutionController(AdjointSolutionController):
    def __init__(self, main_model_part, parameters):
        super().__init__(main_model_part, parameters)

        default_parameters = Kratos.Parameters("""{
            "solution_control_method"  : "reset_to_zero",
            "reset_interval_steps"     : 1
        }""")

        self.parameters.ValidateAndAssignDefaults(default_parameters)

        self.reset_interval_steps = self.parameters["reset_interval_steps"].GetInt()

    def FinalizeSolutionStep(self):
        if (self.inteval_steps >= self.reset_interval_steps):
            self.check_step = self.main_model_part.ProcessInfo[Kratos.STEP]

            # reset the adjoint variables
            var_utils = Kratos.VariableUtils()

            var_utils.SetHistoricalVariableToZero(KratosCFD.ADJOINT_FLUID_VECTOR_1, self.main_model_part.Nodes)
            var_utils.SetHistoricalVariableToZero(KratosCFD.ADJOINT_FLUID_VECTOR_2, self.main_model_part.Nodes)
            var_utils.SetHistoricalVariableToZero(KratosCFD.ADJOINT_FLUID_VECTOR_3, self.main_model_part.Nodes)
            var_utils.SetHistoricalVariableToZero(KratosCFD.AUX_ADJOINT_FLUID_VECTOR_1, self.main_model_part.Nodes)
            var_utils.SetHistoricalVariableToZero(KratosCFD.ADJOINT_FLUID_SCALAR_1, self.main_model_part.Nodes)

            var_utils.SetHistoricalVariableToZero(KratosRANS.RANS_SCALAR_1_ADJOINT_1, self.main_model_part.Nodes)
            var_utils.SetHistoricalVariableToZero(KratosRANS.RANS_SCALAR_1_ADJOINT_2, self.main_model_part.Nodes)
            var_utils.SetHistoricalVariableToZero(KratosRANS.RANS_SCALAR_1_ADJOINT_3, self.main_model_part.Nodes)
            var_utils.SetHistoricalVariableToZero(KratosRANS.RANS_AUX_ADJOINT_SCALAR_1, self.main_model_part.Nodes)

            var_utils.SetHistoricalVariableToZero(KratosRANS.RANS_SCALAR_2_ADJOINT_1, self.main_model_part.Nodes)
            var_utils.SetHistoricalVariableToZero(KratosRANS.RANS_SCALAR_2_ADJOINT_2, self.main_model_part.Nodes)
            var_utils.SetHistoricalVariableToZero(KratosRANS.RANS_SCALAR_2_ADJOINT_3, self.main_model_part.Nodes)
            var_utils.SetHistoricalVariableToZero(KratosRANS.RANS_AUX_ADJOINT_SCALAR_2, self.main_model_part.Nodes)

            Kratos.Logger.PrintInfo(self.__class__.__name__, "Resetted adjoint solution to zero.")

class AdjointRANSSolver(CoupledRANSSolver):
    def __init__(self, model, custom_settings):
        # adjoint settings is validated here without going through
        # GetDefaultParameters, because GetDefaultParameters are
        # used to validate base class primal problem settings.
        default_settings = Kratos.Parameters("""
        {
            "solver_type" : "AdjointRANSSolver",
            "primal_problem_project_parameters_file_name": "PLEASE_SPECIFY_PRIMAL_PROJECT_PARAMETERS_FILE_NAME",
            "response_function_settings" : {
                "response_type" : "drag"
            },
            "sensitivity_settings" : {},
            "linear_solver_settings" : {
                "solver_type" : "amgcl"
            },
            "adjoint_solution_controller_settings": {
                "solution_control_method"  : "none"
            },
            "echo_level": 0,
            "compute_response_function_interpolation_error" : false,
            "response_function_interpolation_error_settings": {}
        }""")

        self.adjoint_settings = custom_settings
        self.adjoint_settings.ValidateAndAssignDefaults(default_settings)

        # open primal problem project parameters
        with open(self.adjoint_settings["primal_problem_project_parameters_file_name"].GetString(), "r") as file_input:
            self.primal_problem_project_parameters = Kratos.Parameters(file_input.read())

        self.primal_problem_solver_settings = self.primal_problem_project_parameters["solver_settings"]
        self._validate_settings_in_baseclass=True # To be removed eventually
        super().__init__(model, self.primal_problem_solver_settings)

        self.adjoint_element_map = {
            ("RansCircularConvectionRFC",) : "RansCircularConvectionRFCAdjoint",
            ("RansDiffusionRFC",) : "RansDiffusionRFCAdjoint",
            ("QSVMS", ) : "QSVMSAdjoint",
            ("QSVMS", "RansKEpsilonKRFC", "RansKEpsilonEpsilonRFC") : "RansKEpsilonQSVMSRFCAdjoint",
            ("QSVMS", "RansKOmegaKRFC", "RansKOmegaOmegaRFC") : "RansKOmegaQSVMSRFCAdjoint",
            ("QSVMS", "RansKOmegaSSTKRFC", "RansKOmegaSSTOmegaRFC") : "RansKOmegaSSTQSVMSRFCAdjoint"
        }

        self.adjoint_condition_map = {
            ('',) : "RansScalarEquationAdjoint",
            # pure cfd conditions
            ("RansVMSMonolithicKBasedWall",) : "AdjointMonolithicWallCondition",
            # k-epsilon conditions
            ("RansVMSMonolithicKBasedWall", "", "RansKEpsilonEpsilonKBasedWall") : "RansKEpsilonVMSKBasedEpsilonKBasedWallAdjoint",
            ("RansVMSMonolithicKBasedWall", "", "RansKEpsilonEpsilonUBasedWall") : "RansKEpsilonVMSKBasedEpsilonUBasedWallAdjoint",
            # k-omega conditions
            ("RansVMSMonolithicKBasedWall", "", "RansKOmegaOmegaKBasedWall") : "RansKOmegaVMSKBasedOmegaKBasedWallAdjoint",
            ("RansVMSMonolithicKBasedWall", "", "RansKOmegaOmegaUBasedWall") : "RansKOmegaVMSKBasedOmegaUBasedWallAdjoint",
            # k-omega-sst conditions
            ("RansVMSMonolithicKBasedWall", "", "RansKOmegaSSTOmegaKBasedWall") : "RansKOmegaSSTVMSKBasedOmegaKBasedWallAdjoint",
            ("RansVMSMonolithicKBasedWall", "", "RansKOmegaSSTOmegaUBasedWall") : "RansKOmegaSSTVMSKBasedOmegaUBasedWallAdjoint"
        }

        self.min_buffer_size = 2
        self.main_model_part.ProcessInfo[KratosRANS.RANS_IS_STEADY] = self.is_steady

        adjoint_solution_controller_settings = self.adjoint_settings["adjoint_solution_controller_settings"]
        if (adjoint_solution_controller_settings.Has("solution_control_method")):
            adjoint_solution_controller_type = adjoint_solution_controller_settings["solution_control_method"].GetString()
            solution_control_methods = [
                "none",
                "reset_to_zero",
                "diffusion"
            ]
            if (not adjoint_solution_controller_type in solution_control_methods):
                msg = "Unsupported adjoint solution control method provided. Requested \"adjoint_solution_controller_settings\" = \"{:s}\".".format(adjoint_solution_controller_type)
                msg += "Supported controller types are:\n\t" + "\n\t".join(solution_control_methods)
                raise RuntimeError(msg)

            if (adjoint_solution_controller_type == "none"):
                self.adjoint_solution_controller = AdjointSolutionController(self.main_model_part, adjoint_solution_controller_settings)
            elif (adjoint_solution_controller_type == "reset_to_zero"):
                self.adjoint_solution_controller = ResetToZeroAdjointSolutionController(self.main_model_part, adjoint_solution_controller_settings)
            elif (adjoint_solution_controller_type == "diffusion"):
                self.adjoint_solution_controller = DiffusionSolutionController(self.main_model_part, adjoint_solution_controller_settings)
        else:
            self.adjoint_solution_controller = AdjointSolutionController(self.main_model_part, adjoint_solution_controller_settings)

        self.adjoint_list_of_time_steps = None
        self.compute_transient_response_function_interpolation_error = False
        if (self.adjoint_settings["compute_response_function_interpolation_error"].GetBool()):
            time_scheme_settings = self.settings["time_scheme_settings"]
            time_scheme_type = time_scheme_settings["scheme_type"].GetString()
            if (time_scheme_type == "bossak"):
                self.adjoint_list_of_time_steps = []
                self.compute_transient_response_function_interpolation_error = True
                if (InputFactory is None):
                    raise Exception("Error importing single_mesh_temporal_input_process from HDF5Application which is required for response function interpolation error computation with bossak time integration scheme.")
                if (OutputFactory is None):
                    raise Exception("Error importing single_mesh_temporal_output_process from HDF5Application which is required for response function interpolation error computation with bossak time integration scheme.")
            elif (time_scheme_type != "steady"):
                raise Exception("Unsupported time scheme type provided for response function interpolation error computation [ time_scheme_type = {:s} ]. Supported time scheme types are: bossak.".format(time_scheme_type))

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Construction of AdjointRANSSolver finished.")

    def AddVariables(self):
        # add primal variables
        super().AddVariables()

        # add adjoint specific variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_1)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_2)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_3)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.AUX_ADJOINT_FLUID_VECTOR_1)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_SCALAR_1)

        self.main_model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_SCALAR_1_ADJOINT_1)
        self.main_model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_SCALAR_1_ADJOINT_2)
        self.main_model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_SCALAR_1_ADJOINT_3)
        self.main_model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_AUX_ADJOINT_SCALAR_1)

        self.main_model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_SCALAR_2_ADJOINT_1)
        self.main_model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_SCALAR_2_ADJOINT_2)
        self.main_model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_SCALAR_2_ADJOINT_3)
        self.main_model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_AUX_ADJOINT_SCALAR_2)

        # add sensitivity variables
        self.main_model_part.AddNodalSolutionStepVariable(Kratos.SHAPE_SENSITIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(Kratos.NORMAL_SENSITIVITY)

        if (self.compute_transient_response_function_interpolation_error):
            self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUX_ADJOINT_FLUID_VECTOR_1)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Adjoint fluid solver variables added correctly.")

    def AddDofs(self):
        Kratos.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_VECTOR_1_X, self.main_model_part)
        Kratos.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_VECTOR_1_Y, self.main_model_part)
        Kratos.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_VECTOR_1_Z, self.main_model_part)
        Kratos.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_SCALAR_1, self.main_model_part)

        Kratos.VariableUtils().AddDof(KratosRANS.RANS_SCALAR_1_ADJOINT_1, self.main_model_part)
        Kratos.VariableUtils().AddDof(KratosRANS.RANS_SCALAR_2_ADJOINT_1, self.main_model_part)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Adjoint fluid solver DOFs added correctly.")

    def PrepareModelPart(self):
        if not self.is_restarted():
            ## Set fluid properties from materials json file
            materials_imported = self._SetPhysicalProperties()
            if not materials_imported:
                Kratos.Logger.PrintWarning(self.__class__.__name__, "Material properties have not been imported. Check \'material_import_settings\' in your ProjectParameters.json.")
            ## Replace default elements and conditions
            self._ReplaceElementsAndConditions()
            ## remove fluid_computational_model_part if it exists. It will be there if adaptive mesh refinement is used.
            if (self.main_model_part.HasSubModelPart("fluid_computational_model_part")):
                self.main_model_part.RemoveSubModelPart("fluid_computational_model_part")
            ## Executes the check and prepare model process
            self._ExecuteCheckAndPrepare()
            ## Set buffer size
            self.main_model_part.SetBufferSize(self.min_buffer_size)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Model reading finished.")

    def Initialize(self):
        if (IsDistributedRun()):
            self.communicator = KratosTrilinos.CreateCommunicator()
        else:
            self.communicator = None

        # check whether periodic conditions are present, if so replace them with dummy LineCondition2D2N
        for condition in self.main_model_part.Conditions:
            if (condition.Is(Kratos.PERIODIC)):
                id = condition.Id
                node_ids = [node.Id for node in condition.GetGeometry()]
                properties = condition.Properties
                self.main_model_part.RemoveCondition(condition)
                self.main_model_part.CreateNewCondition("LineCondition2D2N", id, node_ids, properties)

        domain_size = self.main_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        CalculateNormalsOnConditions(self.main_model_part)
        Kratos.NormalCalculationUtils().CalculateNormalShapeDerivativesOnSimplex(self.main_model_part.Conditions, domain_size)

        # Construct and set the solution strategy
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Initialize the strategy and adjoint utilities
        solution_strategy.Initialize()
        self.GetResponseFunction().Initialize()
        self.GetSensitivityBuilder().Initialize()

        model_part = self.GetComputingModelPart()

        # initialize variables for response function interpolation error calculation
        solving_variables_list = self.formulation.GetSolvingVariables()
        number_of_dofs_per_node = len(solving_variables_list)
        zero_vector = Kratos.Vector(number_of_dofs_per_node, 0.0)

        variable_utils = Kratos.VariableUtils()
        variable_utils.SetNonHistoricalVariable(Kratos.RESPONSE_FUNCTION_INTERPOLATION_ERROR, zero_vector, model_part.Nodes)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def InitializeSolutionStep(self):
        self._GetSolutionStrategy().InitializeSolutionStep()
        self.GetResponseFunction().InitializeSolutionStep()
        self.GetSensitivityBuilder().InitializeSolutionStep()

        if (not self.is_steady):
            self.adjoint_solution_controller.InitializeSolutionStep()

    def Predict(self):
        self._GetSolutionStrategy().Predict()

    def SolveSolutionStep(self):
        return self._GetSolutionStrategy().SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self._GetSolutionStrategy().FinalizeSolutionStep()
        self.GetResponseFunction().FinalizeSolutionStep()

        original_fractional_step = self.main_model_part.ProcessInfo[Kratos.FRACTIONAL_STEP]
        self.main_model_part.ProcessInfo[Kratos.FRACTIONAL_STEP] = 200
        self.GetSensitivityBuilder().UpdateSensitivities()
        self.main_model_part.ProcessInfo[Kratos.FRACTIONAL_STEP] = original_fractional_step

        self.GetSensitivityBuilder().FinalizeSolutionStep()

        if (not self.is_steady):
            self.adjoint_solution_controller.FinalizeSolutionStep()

        if (self.compute_transient_response_function_interpolation_error):
            # save response function interpolation data
            if (not hasattr(self, "response_function_interpolation_data_output_process")):
                # hdf5 file settings
                hdf5_settings = Kratos.Parameters("""{
                    "Parameters": {
                        "model_part_name": "",
                        "file_settings": {
                            "file_name": "response_function_interpolation_data/adjoint_<model_part_name>-<time>.h5",
                            "time_format": "0.6f",
                            "max_files_to_keep": "unlimited",
                            "file_access_mode": "truncate",
                            "echo_level": 0
                        },
                        "output_time_settings": {
                            "step_frequency": 1,
                            "time_frequency": 1.0
                        },
                        "element_data_value_settings" : {
                            "list_of_variables": [
                            "RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1",
                            "RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_2"
                            ]
                        }
                    }
                }""")
                hdf5_settings["Parameters"]["model_part_name"].SetString(self.GetComputingModelPart().FullName())
                hdf5_settings["Parameters"]["output_time_settings"]["time_frequency"].SetDouble(self._ComputeDeltaTime())
                self.response_function_interpolation_data_output_process = OutputFactory(hdf5_settings, self.GetComputingModelPart().GetModel())
                self.response_function_interpolation_data_output_process.ExecuteInitialize()

                Kratos.Logger.PrintInfo(self.__class__.__name__, "Created HDF5 temporal output process for response function interpolation error computation.")

            self.adjoint_list_of_time_steps.append(self.GetComputingModelPart().ProcessInfo[Kratos.TIME])
            self.response_function_interpolation_data_output_process.ExecuteFinalizeSolutionStep()

    def Finalize(self):
        super().Finalize()

        if (self.compute_transient_response_function_interpolation_error):
            model_part = self.GetComputingModelPart()

            # initialize variables for response function interpolation error calculation
            solving_variables_list = self.formulation.GetSolvingVariables()
            number_of_dofs_per_node = len(solving_variables_list)
            zero_vector = Kratos.Vector(number_of_dofs_per_node, 0.0)
            max_vector = Kratos.Vector(number_of_dofs_per_node, 1e+50)

            variable_utils = Kratos.VariableUtils()
            # set initial error to zero assuming, we dont have any errors before start of the error calculation
            # this corresponds to 0th step errors
            variable_utils.SetNonHistoricalVariable(Kratos.RESPONSE_FUNCTION_INTERPOLATION_ERROR, max_vector, model_part.Nodes)
            variable_utils.SetNonHistoricalVariable(Kratos.RESPONSE_FUNCTION_INTERPOLATION_ERROR, zero_vector, model_part.Elements)
            variable_utils.SetNonHistoricalVariable(KratosRANS.RANS_RESPONSE_FUNCTION_DOFS_INTERPOLATION_ERROR_RATE, zero_vector, model_part.Elements)

            hdf5_settings = Kratos.Parameters("""
            {
                "Parameters": {
                    "model_part_name": "",
                    "file_settings": {
                        "file_name": "response_function_interpolation_data/adjoint_<model_part_name>-<time>.h5",
                        "time_format": "0.6f",
                        "echo_level": 1
                    },
                    "element_data_value_settings" : {
                        "list_of_variables": [
                            "RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1",
                            "RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_2"
                        ]
                    }
                }
            }""")
            hdf5_settings["Parameters"]["model_part_name"].SetString(self.GetComputingModelPart().FullName())
            response_function_interpolation_data_input_process = InputFactory(hdf5_settings, self.GetComputingModelPart().GetModel())
            response_function_interpolation_data_input_process.ExecuteInitialize()
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Created HDF5 temporal input process for response function interpolation error computation.")

            # calculate constants
            delta_time = self._ComputeDeltaTime()
            gamma = Kratos.TimeDiscretization.Bossak(model_part.ProcessInfo[Kratos.BOSSAK_ALPHA]).GetGamma()

            vtk_settings = Kratos.Parameters("""{
                "model_part_name"                             : "PLEASE_SPECIFY_MODEL_PART_NAME",
                "file_format"                                 : "binary",
                "output_precision"                            : 7,
                "output_control_type"                         : "step",
                "output_interval"                             : 1.0,
                "output_sub_model_parts"                      : false,
                "output_path"                                 : "goal_error_vtk_output",
                "nodal_data_value_variables"                  : [
                    "RESPONSE_FUNCTION_INTERPOLATION_ERROR"
                ],
                "element_data_value_variables"                : [
                    "RESPONSE_FUNCTION_INTERPOLATION_ERROR",
                    "RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1",
                    "RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_2",
                    "RANS_RESPONSE_FUNCTION_DOFS_INTERPOLATION_ERROR_RATE"
                ]
            }""")

            vtk_settings["model_part_name"].SetString(self.main_model_part.FullName())
            vtk_output = Kratos.VtkOutput(self.main_model_part, vtk_settings)

            # calculate interpolation errors in forward time
            model_part.ProcessInfo[Kratos.TIME] = 0.0

            self.adjoint_list_of_time_steps = sorted(self.adjoint_list_of_time_steps)
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Running forward interpolation error calculation with {:f} relaxation.".format(self.interpolation_error_relaxation))
            for index, current_time_step in enumerate(self.adjoint_list_of_time_steps):
                model_part.ProcessInfo[Kratos.TIME] = current_time_step
                model_part.ProcessInfo[Kratos.STEP] = index + 1

                # read RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1 and RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_2
                response_function_interpolation_data_input_process.ExecuteInitializeSolutionStep()

                # compute interpolation error in forward time
                KratosRANS.RansAdjointUtilities.CalculateTransientReponseFunctionInterpolationError(
                    model_part,
                    gamma,
                    delta_time)

                vtk_output.PrintOutput()

                if (self.echo_level > 0):
                    Kratos.Logger.PrintInfo(self.__class__.__name__, "Computed response function interpolation error at {:f}s".format(current_time_step))

    def Check(self):
        self._GetSolutionStrategy().Check()

    def Clear(self):
        self._GetSolutionStrategy().Clear()

    def _ReplaceElementsAndConditions(self):
        domain_size = self.main_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        self.new_elem_name = "{:s}{:d}D{:d}N".format(
            self.adjoint_element_map[tuple(self.formulation.GetElementNames())],
            domain_size,
            domain_size + 1
        )

        self.new_cond_name = "{:s}{:d}D{:d}N".format(
            self.adjoint_condition_map[tuple(self.formulation.GetConditionNames())],
            domain_size,
            domain_size
        )

        ## Set the element and condition names in the Json parameters
        self.settings.AddValue("element_replace_settings", Kratos.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(self.new_elem_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(self.new_cond_name)

        ## Call the replace elements and conditions process
        Kratos.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

    def _ComputeDeltaTime(self):
        if self.settings["time_stepping"]["automatic_time_step"].GetBool():
            raise Exception("Automatic time stepping is not supported by adjoint RANS solver.")

        delta_time = self.settings["time_stepping"]["time_step"].GetDouble()
        return -1.0 * delta_time

    # TODO: I THINK THIS SHOULD BE MOVED TO THE BASE PYTHON SOLVER
    def is_restarted(self):
        # this function avoids the long call to ProcessInfo and is also safer
        # in case the detection of a restart is changed later
        return self.main_model_part.ProcessInfo[Kratos.IS_RESTARTED]

    def _GetScheme(self):
        if not hasattr(self, '_scheme'):
            self._scheme = self._CreateScheme()
        return self._scheme

    def _GetLinearSolver(self):
        if not hasattr(self, '_linear_solver'):
            self._linear_solver = self._CreateLinearSolver()
        return self._linear_solver

    def _GetBuilderAndSolver(self):
        if not hasattr(self, '_builder_and_solver'):
            self._builder_and_solver = self._CreateBuilderAndSolver()
        return self._builder_and_solver

    def _GetSolutionStrategy(self):
        if not hasattr(self, '_solution_strategy'):
            self._solution_strategy = self._CreateSolutionStrategy()
        return self._solution_strategy

    def GetResponseFunction(self):
        if not hasattr(self, '_response_function'):
            self._response_function = self.__CreateResponseFunction()
        return self._response_function

    def GetSensitivityBuilder(self):
        if not hasattr(self, '_sensitivity_builder'):
            self._sensitivity_builder = self.__CreateSensitivityBuilder()
        return self._sensitivity_builder

    def __CreateResponseFunction(self):
        domain_size = self.main_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        response_type = self.adjoint_settings["response_function_settings"]["response_type"].GetString()

        if domain_size == 2:
            drag_response_function_type = KratosCFD.DragResponseFunction2D
            drag_frequency_response_function_type = KratosCFD.DragFrequencyResponseFunction2D
            residual_response_function_type = KratosCFD.ResidualResponseFunction2D
            domain_integrated_3d_vector_magnitude_square_p_mean_response_function_type = KratosCFD.DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction2D
        elif domain_size == 3:
            drag_response_function_type = KratosCFD.DragResponseFunction3D
            drag_frequency_response_function_type = KratosCFD.DragFrequencyResponseFunction3D
            residual_response_function_type = KratosCFD.ResidualResponseFunction3D
            domain_integrated_3d_vector_magnitude_square_p_mean_response_function_type = KratosCFD.DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction3D
        else:
            raise RuntimeError("Invalid DOMAIN_SIZE: " + str(domain_size))

        if response_type == "drag":
            response_function = drag_response_function_type(
                self.adjoint_settings["response_function_settings"]["custom_settings"],
                self.main_model_part)
        elif response_type == "norm_square":
            response_function = KratosCFD.VelocityPressureNormSquareResponseFunction(
                self.adjoint_settings["response_function_settings"]["custom_settings"],
                self.main_model_part.GetModel())
        elif response_type == "drag_frequency":
            response_function = drag_frequency_response_function_type(
                self.adjoint_settings["response_function_settings"]["custom_settings"],
                self.main_model_part)
        elif response_type == "residual":
            response_function = residual_response_function_type(
                self.adjoint_settings["response_function_settings"]["custom_settings"],
                self.main_model_part)
        elif response_type == "domain_integrated":
            response_function = KratosCFD.DomainIntegratedResponseFunction(
                self.adjoint_settings["response_function_settings"]["custom_settings"],
                self.main_model_part)
        elif response_type == "domain_integrated_3d_vector_magnitude_square_power_mean":
            response_function = domain_integrated_3d_vector_magnitude_square_p_mean_response_function_type(
                self.adjoint_settings["response_function_settings"]["custom_settings"],
                self.main_model_part)
        else:
            raise Exception("Invalid response_type: " + response_type + ". Available response functions: \'drag\'.")
        return response_function

    def __CreateSensitivityBuilder(self):
        response_function = self.GetResponseFunction()

        domain_size = self.main_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        block_size = len(self.formulation.GetSolvingVariables())

        time_scheme_settings = self.settings["time_scheme_settings"]
        time_scheme_type = time_scheme_settings["scheme_type"].GetString()

        if (time_scheme_type == "steady"):
            self.sensitivity_builder_scheme = KratosCFD.SimpleSteadySensitivityBuilderScheme(domain_size, block_size)
            self.is_steady = True
        elif (time_scheme_type == "bossak"):
            self.sensitivity_builder_scheme = KratosCFD.VelocityBossakSensitivityBuilderScheme(time_scheme_settings["alpha_bossak"].GetDouble(), domain_size, block_size)
            self.is_steady = False

        sensitivity_builder = Kratos.SensitivityBuilder(
            self.adjoint_settings["sensitivity_settings"],
            self.main_model_part,
            response_function,
            self.sensitivity_builder_scheme
            )
        return sensitivity_builder

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        time_scheme = self._GetScheme()
        builder_and_solver = self._GetBuilderAndSolver()
        calculate_reaction_flag = False
        reform_dof_set_at_each_step = False
        calculate_norm_dx_flag = False
        move_mesh_flag = False
        return GetKratosObjectPrototype("ResidualBasedLinearStrategy")(
            computing_model_part,
            time_scheme,
            builder_and_solver,
            calculate_reaction_flag,
            reform_dof_set_at_each_step,
            calculate_norm_dx_flag,
            move_mesh_flag)

    def _CreateLinearSolver(self):
        linear_solver_factory = GetKratosObjectPrototype("LinearSolverFactory")
        return linear_solver_factory(self.adjoint_settings["linear_solver_settings"])

    def _CreateScheme(self):
        response_function = self.GetResponseFunction()

        domain_size = self.main_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        block_size = len(self.formulation.GetSolvingVariables())

        self.element_refinement_process = None
        if (self.adjoint_settings["compute_response_function_interpolation_error"].GetBool()):
            refinement_parameters = self.adjoint_settings["response_function_interpolation_error_settings"]
            if (refinement_parameters.Has("relaxation_factor")):
                self.interpolation_error_relaxation = refinement_parameters["relaxation_factor"].GetDouble()
                refinement_parameters.RemoveValue("relaxation_factor")
            else:
                self.interpolation_error_relaxation = 1.0
            self._CheckAndAddValue(refinement_parameters, "model_part_name", self.GetComputingModelPart().FullName())
            self._CheckAndAddValue(refinement_parameters, "refined_element_name", self.new_elem_name)
            self._CheckAndAddValue(refinement_parameters, "refined_condition_name", self.new_cond_name)

            if (not refinement_parameters.Has("nodal_interpolation_settings")):
                refinement_parameters.AddEmptyValue("nodal_interpolation_settings")

            nodal_interpolation_settings = refinement_parameters["nodal_interpolation_settings"]
            list_of_variable_names = []
            # add primal variables
            for variable in self.formulation.GetSolvingVariables():
                list_of_variable_names.extend([var.Name() for var in GetTimeDerivativeVariablesRecursively(variable)])

            # add primal auxiliary variables
            self._CheckAndAddVariableIfExisting(list_of_variable_names, Kratos.NORMAL)
            self._CheckAndAddVariableIfExisting(list_of_variable_names, Kratos.DISTANCE)
            self._CheckAndAddVariableIfExisting(list_of_variable_names, Kratos.BODY_FORCE)

            # add adjoint variables
            list_of_variable_names.extend([
                "ADJOINT_FLUID_VECTOR_1",
                "ADJOINT_FLUID_SCALAR_1",
                "RANS_SCALAR_1_ADJOINT_1",
                "RANS_SCALAR_2_ADJOINT_1"
            ])
            self._CheckAndAddValueToList(nodal_interpolation_settings, "historical_variables_list", list_of_variable_names)
            self._CheckAndAddValueToList(nodal_interpolation_settings, "non_hitsorical_variables_list", ["RELAXED_ACCELERATION"])
            self._CheckAndAddValueToList(nodal_interpolation_settings, "flags_list", ["SLIP", "INLET", "OUTLET", "STRUCTURE"])

            self.element_refinement_process = KratosCFD.ElementRefinementProcess(self.model, refinement_parameters)

        scheme_type = self.settings["time_scheme_settings"]["scheme_type"].GetString()
        if scheme_type == "bossak":
            self.adjoint_solution_controller.AddTimeSchemeParameters(self.settings["time_scheme_settings"])
            scheme = GetKratosObjectPrototype("VelocityBossakAdjointScheme")(
                self.settings["time_scheme_settings"],
                response_function,
                domain_size,
                block_size,
                self.element_refinement_process,
                self.adjoint_settings["echo_level"].GetInt())
        elif scheme_type == "steady":
            scheme = GetKratosObjectPrototype("SimpleSteadyAdjointScheme")(
                response_function,
                domain_size,
                block_size,
                self.element_refinement_process)
        else:
            raise Exception("Invalid scheme_type: " + scheme_type)

        return scheme

    def _CreateBuilderAndSolver(self):
        linear_solver = self._GetLinearSolver()
        return CreateBlockBuilderAndSolver(
            linear_solver,
            self.settings["consider_periodic_conditions"].GetBool(),
            self.communicator)

    def _CheckAndAddValue(self, parameters, check_string, required_value):
        if (not parameters.Has(check_string)):
            parameters.AddEmptyValue(check_string)
        parameters[check_string].SetString(required_value)

    def _CheckAndAddValueToList(self, parameters, check_string, required_values):
        if (parameters.Has(check_string)):
            values = parameters[check_string].GetStringArray()
            for required_value in required_values:
                if (required_value not in values):
                    values.append(required_value)

            parameters[check_string].SetStringArray(values)
        else:
            parameters.AddEmptyValue(check_string)
            parameters[check_string].SetStringArray(required_values)

    def _CheckAndAddVariableIfExisting(self, list_of_variable_names, variable):
        if (self.GetComputingModelPart().HasNodalSolutionStepVariable(variable)):
            list_of_variable_names.append(variable.Name())