from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

## Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver

from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
from KratosMultiphysics.FluidDynamicsApplication import check_and_prepare_model_process_fluid

def CreateSolver(model, custom_settings):
    return NavierStokesCompressibleExplicitSolver(model, custom_settings)

class NavierStokesCompressibleExplicitSolver(FluidSolver):
    def __init__(self, model, custom_settings):
        # Call base fluid solver constructor
        self._validate_settings_in_baseclass=True # To be removed eventually
        super(NavierStokesCompressibleExplicitSolver,self).__init__(model,custom_settings)

        # Define the formulation settings
        self.element_name = "CompressibleNavierStokesExplicit"
        if custom_settings["domain_size"].GetInt() == 2:
            self.condition_name = "LineCondition" # TODO: We need to create a Compressible NS condition (now using the base ones)
        elif custom_settings["domain_size"].GetInt() == 3:
            self.condition_name = "SurfaceCondition" # TODO: We need to create a Compressible NS condition (now using the base ones)
        else:
            err_msg = "Wrong domain size "
            raise Exception(err_msg)
        self.min_buffer_size = 2
        self.element_has_nodal_properties = False # Note that DENSITY is nodally stored but considered as a DOF

        KratosMultiphysics.Logger.PrintInfo("::[NavierStokesCompressibleExplicitSolver]:: ","Construction of NavierStokesCompressibleExplicitSolver finished.")

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "compressible_solver_from_defaults",
            "model_part_name": "FluidModelPart",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "",
                "reorder": false
            },
            "material_import_settings": {
                "materials_filename": "FluidMaterials.json"
            },
            "echo_level": 1,
            "time_order": 2,
            "move_mesh_flag": false,
            "shock_capturing": true,
            "compute_reactions": false,
            "reform_dofs_at_each_step" : false,
            "assign_neighbour_elements_to_conditions": true,
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : true,
                "CFL_number"          : 1.0,
                "minimum_delta_time"  : 1.0e-8,
                "maximum_delta_time"  : 1.0e-2
            },
            "use_oss" : true
        }""")

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def AddVariables(self):
        # Add DOF variables (formulation written in conservative form) and reactions
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY) # Density DOF
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENTUM) # Momentum DOF
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TOTAL_ENERGY) # Total energy DOF
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.REACTION_DENSITY) # Density DOF reaction
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION) # Momentum DOF reaction
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.REACTION_ENERGY) # Total energy DOF reaction

        # Required variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.MASS_SOURCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.HEAT_SOURCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        # Post-process variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.MACH)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)

        KratosMultiphysics.Logger.PrintInfo("::[NavierStokesCompressibleExplicitSolver]:: ","Explicit compressible fluid solver variables added correctly")

    def AddDofs(self):
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DENSITY, KratosFluid.REACTION_DENSITY, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MOMENTUM_X, KratosMultiphysics.REACTION_X, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MOMENTUM_Y, KratosMultiphysics.REACTION_Y, self.main_model_part)
        if domain_size == 3:
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MOMENTUM_Z, KratosMultiphysics.REACTION_Z, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.TOTAL_ENERGY, KratosFluid.REACTION_ENERGY, self.main_model_part)

    def Initialize(self):
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.OSS_SWITCH] = int(self.settings["use_oss"].GetBool())
        self.GetComputingModelPart().ProcessInfo[KratosFluid.SHOCK_CAPTURING_SWITCH] = int(self.settings["shock_capturing"].GetBool())

        self.solver = self._get_solution_strategy()
        self.solver.SetEchoLevel(self.settings["echo_level"].GetInt())
        self.solver.Initialize()

        KratosMultiphysics.Logger.PrintInfo("::[NavierStokesCompressibleExplicitSolver]:: ","Explicit compressible fluid solver initialization finished.")

    def _get_solution_strategy(self):
        if not hasattr(self, '_solution_strategy'):
            self._solution_strategy = self._create_solution_strategy()
        return self._solution_strategy

    def _create_solution_strategy(self):
        self.computing_model_part = self.GetComputingModelPart()
        strategy_settings = KratosMultiphysics.Parameters('''{}''')
        strategy_settings.AddEmptyValue("rebuild_level").SetInt(0 if self.settings["reform_dofs_at_each_step"].GetBool() else 1)
        strategy_settings.AddEmptyValue("move_mesh_flag").SetBool(self.settings["move_mesh_flag"].GetBool())
        strategy_settings.AddEmptyValue("shock_capturing").SetBool(self.settings["shock_capturing"].GetBool())

        strategy = KratosFluid.CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4(
            self.computing_model_part,
            strategy_settings)

        return strategy

    def _CreateEstimateDtUtility(self):
        """This method overloads FluidSolver in order to enforce:
        ```
        self.settings["time_stepping"]["consider_compressibility_in_CFL"] == True
        ```
        """
        if self.settings["time_stepping"].Has("consider_compressibility_in_CFL"):
            KratosMultiphysics.Logger.PrintWarning("", "User-specifed consider_compressibility_in_CFL will be overriden with TRUE")
        else:
            self.settings["time_stepping"].AddEmptyValue("consider_compressibility_in_CFL")

        self.settings["time_stepping"]["consider_compressibility_in_CFL"].SetBool(True)

        estimate_dt_utility = KratosFluid.EstimateDtUtility(
                self.GetComputingModelPart(),
                self.settings["time_stepping"])

        return estimate_dt_utility
