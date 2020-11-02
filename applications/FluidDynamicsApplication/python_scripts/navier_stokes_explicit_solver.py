# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication

## Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver

from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
from KratosMultiphysics.FluidDynamicsApplication import check_and_prepare_model_process_fluid

def CreateSolver(model, custom_settings):
    return NavierStokesExplicitSolver(model, custom_settings)

class NavierStokesExplicitSolver(FluidSolver):
    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass = True # To be removed eventually
        super(NavierStokesExplicitSolver,self).__init__(model,custom_settings)

        if custom_settings["formulation"]["element_type"].GetString() != "QSNavierStokesExplicit":
            raise Exception("NavierStokesExplicitSolver only accepts QSNavierStokesExplicit as the \"element_type\" in \"formulation\"")

        self.element_name = custom_settings["formulation"]["element_type"].GetString()
        self.condition_name = custom_settings["formulation"]["condition_type"].GetString()

        self.min_buffer_size = 2
        self.element_has_nodal_properties = False # TODO: check
        self.compute_reactions = self.settings["compute_reactions"].GetBool() # TODO: check

        if self.settings["use_oss"].GetBool():
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH,1)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesExplicitSolver finished.")

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "FractionalStep",
            "time_integration_method": "explicit",
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
            "move_mesh_flag": false,
            "compute_reactions": false,
            "reform_dofs_at_each_step" : false,
            "assign_neighbour_elements_to_conditions": false,
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step": false,
                "CFL_number": 1.0,
                "minimum_delta_time": 1.0e-8,
                "maximum_delta_time": 1.0e-2
            },
            "use_oss": true,
            "dynamic_tau": 1.0,
            "formulation": {
                "element_type": "QSNavierStokesExplicit",
                "condition_type": "WallCondition"
            }
        }""")

        default_settings.AddMissingParameters(super(NavierStokesExplicitSolver, cls).GetDefaultParameters())
        return default_settings

    def AddVariables(self):
        # Add DOF variables (formulation written in conservative form) and reactions
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Explicit fluid solver variables added correctly")

    def Initialize(self):
        self.solver = self._get_solution_strategy()
        self.solver.SetEchoLevel(self.settings["echo_level"].GetInt())
        self.solver.Initialize()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Explicit fluid solver initialization finished.")

    def _get_solution_strategy(self):
        if not hasattr(self, '_solution_strategy'):
            self._solution_strategy = self._create_solution_strategy()
        return self._solution_strategy

    def _create_solution_strategy(self):
        self.computing_model_part = self.GetComputingModelPart()
        strategy_settings = KratosMultiphysics.Parameters('''{}''')
        strategy_settings.AddEmptyValue("rebuild_level").SetInt(0 if self.settings["reform_dofs_at_each_step"].GetBool() else 1)
        strategy_settings.AddEmptyValue("move_mesh_flag").SetBool(self.settings["move_mesh_flag"].GetBool())

        strategy = FluidDynamicsApplication.NavierStokesExplicitSolvingStrategyRungeKutta4(
            self.computing_model_part,
            strategy_settings)

        return strategy
