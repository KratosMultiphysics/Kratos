# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.adjoint_fluid_solver import AdjointFluidSolver

def CreateSolver(main_model_part, custom_settings):
    return AdjointVMSMonolithicSolver(main_model_part, custom_settings)

class AdjointVMSMonolithicSolver(AdjointFluidSolver):

    @classmethod
    def GetDefaultParameters(cls):
        # default settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type" : "adjoint_vmsmonolithic_solver",
            "model_part_name": "",
            "domain_size": -1,
            "scheme_settings" : {
                "scheme_type" : "bossak"
            },
            "response_function_settings" : {
                "response_type" : "drag"
            },
            "sensitivity_settings" : {},
            "model_import_settings" : {
                "input_type"     : "mdpa",
                "input_filename" : "unknown_name"
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "linear_solver_settings" : {
                "solver_type" : "amgcl"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts"  : [""],
            "no_skin_parts"  : [""],
            "dynamic_tau" : 0.0,
            "oss_switch"  : 0,
            "echo_level"  : 0,
            "time_stepping"               : {
                "automatic_time_step" : false,
                "time_step"           : -0.1
            },
            "consider_periodic_conditions": false,
            "assign_neighbour_elements_to_conditions": true
        }""")

        default_settings.AddMissingParameters(super(AdjointVMSMonolithicSolver, cls).GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):
        super().__init__(model,custom_settings)
        self.element_has_nodal_properties = True

        self.element_name = "VMSAdjointElement"
        self.condition_name = "AdjointMonolithicWallCondition"

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH, self.settings["oss_switch"].GetInt())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of AdjointVMSMonolithicSolver finished.")

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_1)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_2)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_3)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.AUX_ADJOINT_FLUID_VECTOR_1)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_SCALAR_1)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_SENSITIVITY)

        if self.settings["consider_periodic_conditions"].GetBool() == True:
            self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.PATCH_INDEX)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Adjoint fluid solver variables added correctly.")

    def Initialize(self):
        # Construct and set the solution strategy
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Initialize the strategy and adjoint utilities
        solution_strategy.Initialize()
        self.GetResponseFunction().Initialize()
        self.GetSensitivityBuilder().Initialize()

        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.main_model_part.Conditions, domain_size)
        KratosMultiphysics.NormalCalculationUtils().CalculateNormalShapeDerivativesOnSimplex(self.main_model_part.Conditions, domain_size)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def _CreateScheme(self):
        response_function = self.GetResponseFunction()
        scheme_type = self.settings["scheme_settings"]["scheme_type"].GetString()
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        # the schemes used in fluid supports SLIP conditions which rotates element/condition
        # matrices based on nodal NORMAL. Hence, the consistent adjoints also need to
        # rotate adjoint element/condition matrices accordingly. Therefore, following
        # schemes are used.
        if scheme_type == "bossak":
            scheme = KratosCFD.VelocityBossakAdjointScheme(self.settings["scheme_settings"], response_function, domain_size, domain_size + 1)
        elif scheme_type == "steady":
            scheme = KratosCFD.SimpleSteadyAdjointScheme(response_function, domain_size, domain_size + 1)
        else:
            raise Exception("Invalid scheme_type: " + scheme_type)
        return scheme

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        time_scheme = self._GetScheme()
        builder_and_solver = self._GetBuilderAndSolver()
        calculate_reaction_flag = False
        reform_dof_set_at_each_step = False
        calculate_norm_dx_flag = False
        move_mesh_flag = False
        return KratosMultiphysics.ResidualBasedLinearStrategy(
            computing_model_part,
            time_scheme,
            builder_and_solver,
            calculate_reaction_flag,
            reform_dof_set_at_each_step,
            calculate_norm_dx_flag,
            move_mesh_flag)
