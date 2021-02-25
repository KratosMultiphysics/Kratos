# Importing the Kratos Library
import KratosMultiphysics as Kratos

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.RANSApplication as KratosRANS

import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

from KratosMultiphysics.RANSApplication.formulations.utilities import CalculateNormalsOnConditions
from KratosMultiphysics.RANSApplication.coupled_rans_solver import CoupledRANSSolver

def CreateSolver(main_model_part, custom_settings):
    return AdjointRANSSolver(main_model_part, custom_settings)

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
            }
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
            ("QSVMS", ) : "QSVMSAdjoint",
            ("QSVMS", "RansKEpsilonKRFC", "RansKEpsilonEpsilonRFC") : "RansKEpsilonQSVMSRFCAdjoint"
        }

        self.adjoint_condition_map = {
            ("RansVMSMonolithicKBasedWall",) : "AdjointMonolithicWallCondition",
            ("RansVMSMonolithicKBasedWall", "", "RansKEpsilonEpsilonKBasedWall") : "RansKEpsilonVMSKBasedEpsilonKBasedWallAdjoint"
        }

        self.min_buffer_size = 2
        self.main_model_part.ProcessInfo[KratosRANS.RANS_IS_STEADY] = self.is_steady

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
            ## Executes the check and prepare model process
            self._ExecuteCheckAndPrepare()
            ## Set buffer size
            self.main_model_part.SetBufferSize(self.min_buffer_size)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Model reading finished.")

    def Initialize(self):
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

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def InitializeSolutionStep(self):
        self._GetSolutionStrategy().InitializeSolutionStep()
        self.GetResponseFunction().InitializeSolutionStep()
        self.GetSensitivityBuilder().InitializeSolutionStep()

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

    def Check(self):
        self._GetSolutionStrategy().Check()

    def Clear(self):
        self._GetSolutionStrategy().Clear()

    def _ReplaceElementsAndConditions(self):
        domain_size = self.main_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        new_elem_name = "{:s}{:d}D{:d}N".format(
            self.adjoint_element_map[tuple(self.formulation.GetElementNames())],
            domain_size,
            domain_size + 1
        )

        new_cond_name = "{:s}{:d}D{:d}N".format(
            self.adjoint_condition_map[tuple(self.formulation.GetConditionNames())],
            domain_size,
            domain_size
        )

        ## Set the element and condition names in the Json parameters
        self.settings.AddValue("element_replace_settings", Kratos.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_cond_name)

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
        if response_type == "drag":
            if domain_size == 2:
                response_function = KratosCFD.DragResponseFunction2D(
                    self.adjoint_settings["response_function_settings"]["custom_settings"],
                    self.main_model_part)
            elif domain_size == 3:
                response_function = KratosCFD.DragResponseFunction3D(
                    self.adjoint_settings["response_function_settings"]["custom_settings"],
                    self.main_model_part)
            else:
                raise Exception("Invalid DOMAIN_SIZE: " + str(domain_size))
        elif response_type == "norm_square":
            response_function = KratosCFD.VelocityPressureNormSquareResponseFunction(
                self.adjoint_settings["response_function_settings"]["custom_settings"],
                self.main_model_part)
        else:
            raise Exception("Invalid response_type: " + response_type + ". Available response functions: \'drag\'.")
        return response_function

    def __CreateSensitivityBuilder(self):
        response_function = self.GetResponseFunction()

        element_names = tuple(self.formulation.GetElementNames())
        domain_size = self.main_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        if (domain_size == 2):
            if (element_names == (("QSVMS"), )):
                bossak_scheme_type = KratosCFD.VelocityBossakSensitivityBuilderScheme2D
                steady_scheme_type = KratosCFD.SimpleSteadySensitivityBuilderScheme2D
            else:
                bossak_scheme_type = KratosRANS.RansVelocityBossakSensitivityBuilderScheme2D
                steady_scheme_type = KratosRANS.RansSimpleSteadySensitivityBuilderScheme2D
        elif (domain_size == 3):
            if (element_names == (("QSVMS"), )):
                bossak_scheme_type = KratosCFD.VelocityBossakSensitivityBuilderScheme3D
                steady_scheme_type = KratosCFD.SimpleSteadySensitivityBuilderScheme3D
            else:
                bossak_scheme_type = KratosRANS.RansVelocityBossakSensitivityBuilderScheme3D
                steady_scheme_type = KratosRANS.RansSimpleSteadySensitivityBuilderScheme3D
        else:
            raise Exception("Invalid DOMAIN_SIZE: " + str(domain_size))

        time_scheme_settings = self.settings["time_scheme_settings"]
        time_scheme_type = time_scheme_settings["scheme_type"].GetString()

        if (time_scheme_type == "steady"):
            self.sensitivity_builder_scheme = steady_scheme_type()
        elif (time_scheme_type == "bossak"):
            self.sensitivity_builder_scheme = bossak_scheme_type(time_scheme_settings["alpha_bossak"].GetDouble())

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
        return Kratos.ResidualBasedLinearStrategy(
            computing_model_part,
            time_scheme,
            builder_and_solver,
            calculate_reaction_flag,
            reform_dof_set_at_each_step,
            calculate_norm_dx_flag,
            move_mesh_flag)

    def _CreateLinearSolver(self):
        linear_solver_configuration = self.adjoint_settings["linear_solver_settings"]
        return linear_solver_factory.ConstructSolver(linear_solver_configuration)

    def _CreateScheme(self):
        response_function = self.GetResponseFunction()

        element_names = tuple(self.formulation.GetElementNames())
        domain_size = self.main_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        if (domain_size == 2):
            if (element_names == (("QSVMS"), )):
                bossak_scheme_type = KratosCFD.VelocityBossakAdjointScheme2D
                steady_scheme_type = KratosCFD.SimpleSteadyAdjointScheme2D
            else:
                bossak_scheme_type = KratosRANS.RansVelocityBossakAdjointScheme2D
                steady_scheme_type = KratosRANS.RansSimpleSteadyAdjointScheme2D
        elif (domain_size == 3):
            if (element_names == (("QSVMS"), )):
                bossak_scheme_type = KratosCFD.VelocityBossakAdjointScheme3D
                steady_scheme_type = KratosCFD.SimpleSteadyAdjointScheme3D
            else:
                bossak_scheme_type = KratosRANS.RansVelocityBossakAdjointScheme3D
                steady_scheme_type = KratosRANS.RansSimpleSteadyAdjointScheme3D
        else:
            raise Exception("Invalid DOMAIN_SIZE: " + str(domain_size))

        scheme_type = self.settings["time_scheme_settings"]["scheme_type"].GetString()
        domain_size = self.main_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        if scheme_type == "bossak":
            scheme = bossak_scheme_type(self.settings["time_scheme_settings"], response_function)
        elif scheme_type == "steady":
            scheme = steady_scheme_type(response_function)
        else:
            raise Exception("Invalid scheme_type: " + scheme_type)

        return scheme

    def _CreateBuilderAndSolver(self):
        linear_solver = self._GetLinearSolver()
        if self.settings["consider_periodic_conditions"].GetBool():
            builder_and_solver = KratosCFD.ResidualBasedBlockBuilderAndSolverPeriodic(
                linear_solver,
                KratosCFD.PATCH_INDEX)
        else:
            builder_and_solver = Kratos.ResidualBasedBlockBuilderAndSolver(linear_solver)
        return builder_and_solver
