from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.adjoint_fluid_solver import AdjointFluidSolver

from KratosMultiphysics.FluidDynamicsApplication.adjoint_turbulence_model_solver import CreateAdjointTurbulenceModel

def CreateSolver(main_model_part, custom_settings):
    return AdjointVMSMonolithicSolver(main_model_part, custom_settings)

class AdjointVMSMonolithicSolver(AdjointFluidSolver):

    @classmethod
    def GetDefaultSettings(cls):
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
            "adjoint_turbulence_model_solver_settings": {},
            "consider_periodic_conditions": false
        }""")

        default_settings.AddMissingParameters(super(AdjointVMSMonolithicSolver, cls).GetDefaultSettings())
        return default_settings

    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass=True # To be removed eventually
        super(AdjointVMSMonolithicSolver,self).__init__(model,custom_settings)
        self.element_has_nodal_properties = True

        # construct the linear solver
        import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        if not self.settings["adjoint_turbulence_model_solver_settings"].IsEquivalentTo(KratosMultiphysics.Parameters("{}")):
            # if not empty
            self._adjoint_turbulence_model_solver = CreateAdjointTurbulenceModel(model, self.settings["adjoint_turbulence_model_solver_settings"])
            self.element_name = self._adjoint_turbulence_model_solver.GetAdjointElementName()
            self.condition_name = self._adjoint_turbulence_model_solver.GetAdjointConditionName()
        else:
            self.element_name = "VMSAdjointElement"
            if self.settings["domain_size"].GetInt() == 2:
                self.condition_name = "LineCondition"
            elif self.settings["domain_size"].GetInt() == 3:
                self.condition_name = "SurfaceCondition"


        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of AdjointVMSMonolithicSolver finished.")

    def AddVariables(self):

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_1)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_2)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_3)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.AUX_ADJOINT_FLUID_VECTOR_1)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_SCALAR_1)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_SENSITIVITY)

        # Adding variables required for the turbulence modelling
        if hasattr(self, "_adjoint_turbulence_model_solver"):
            self._adjoint_turbulence_model_solver.fluid_model_part = self.main_model_part
            self._adjoint_turbulence_model_solver.AddVariables()

        if self.settings["consider_periodic_conditions"].GetBool() == True:
            self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.PATCH_INDEX)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Adjoint fluid solver variables added correctly.")

    def AddDofs(self):
        super(AdjointVMSMonolithicSolver, self).AddDofs()

        if hasattr(self, "_adjoint_turbulence_model_solver"):
            self._adjoint_turbulence_model_solver.AddDofs()

    def Initialize(self):

        self.computing_model_part = self.GetComputingModelPart()

        self.response_function = self.GetResponseFunction()

        self.sensitivity_builder = KratosMultiphysics.SensitivityBuilder(self.settings["sensitivity_settings"], self.main_model_part, self.response_function)

        if self.settings["scheme_settings"]["scheme_type"].GetString() == "bossak":
            self.time_scheme = KratosMultiphysics.ResidualBasedAdjointBossakScheme(self.settings["scheme_settings"], self.response_function)
        elif self.settings["scheme_settings"]["scheme_type"].GetString() == "steady":
            self.time_scheme = KratosMultiphysics.ResidualBasedAdjointSteadyScheme(self.response_function)
        else:
            raise Exception("invalid scheme_type: " + self.settings["scheme_settings"]["scheme_type"].GetString())

        if self.settings["consider_periodic_conditions"].GetBool() == True:
            builder_and_solver = KratosCFD.ResidualBasedBlockBuilderAndSolverPeriodic(
                    self.linear_solver, KratosCFD.PATCH_INDEX)
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = KratosMultiphysics.ResidualBasedLinearStrategy(self.main_model_part,
                                                                     self.time_scheme,
                                                                     self.linear_solver,
                                                                     builder_and_solver,
                                                                     False,
                                                                     False,
                                                                     False,
                                                                     False)

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())


        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH, self.settings["oss_switch"].GetInt())

        if hasattr(self, "_adjoint_turbulence_model_solver"):
            self._adjoint_turbulence_model_solver.Initialize()

        (self.solver).Initialize()
        (self.response_function).Initialize()
        (self.sensitivity_builder).Initialize()
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def GetResponseFunction(self):
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        if self.settings["response_function_settings"]["response_type"].GetString() == "drag":
            if (domain_size == 2):
                return KratosCFD.DragResponseFunction2D(self.settings["response_function_settings"]["custom_settings"], self.main_model_part)
            elif (domain_size == 3):
                return KratosCFD.DragResponseFunction3D(self.settings["response_function_settings"]["custom_settings"], self.main_model_part)
            else:
                raise Exception("Invalid DOMAIN_SIZE: " + str(domain_size))
        else:
            raise Exception("invalid response_type: " + self.settings["response_function_settings"]["response_type"].GetString())