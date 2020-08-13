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
        # Call base fluid solver constructor
        self._validate_settings_in_baseclass=True # To be removed eventually
        super(NavierStokesExplicitSolver,self).__init__(model,custom_settings)

        # Define the formulation settings
        self.element_name = "SymbolicExplicitQSNavierStokes"
        if custom_settings["domain_size"].GetInt() == 2:
            self.condition_name = "LineCondition"
        elif custom_settings["domain_size"].GetInt() == 3:
            self.condition_name = "SurfaceCondition"
        else:
            err_mgs = "Wrong domain size "
            raise Exception(err_msg)
        self.min_buffer_size = 2
        self.element_has_nodal_properties = False # TODO: check

        KratosMultiphysics.Logger.PrintInfo("::[NavierStokesExplicitSolver]:: ","Construction of NavierStokesExplicitSolver finished.")

    @classmethod
    def GetDefaultSettings(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "Explicit",
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
            "compute_reactions": false,
            "reform_dofs_at_each_step" : false,
            "assign_neighbour_elements_to_conditions": false,
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : false,
                "CFL_number"          : 1.0,
                "minimum_delta_time"  : 1.0e-8,
                "maximum_delta_time"  : 1.0e-2
            },
            "use_oss" : true,
            "transient_parameters" : {
                "dynamic_tau": 1.0
            }
        }""")

        default_settings.AddMissingParameters(super(NavierStokesExplicitSolver, cls).GetDefaultSettings())
        return default_settings

    def AddVariables(self):
        # Add DOF variables (formulation written in conservative form) and reactions
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)

        KratosMultiphysics.Logger.PrintInfo("::[NavierStokesExplicitSolver]:: ","Explicit fluid solver variables added correctly")

    def Initialize(self):
        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            raise Exception("ERROR: _GetAutomaticTimeSteppingUtility out of date")

        # Set the time discretization utility to compute the BDF coefficients
        time_order = self.settings["time_order"].GetInt()
        if time_order == 2:
            self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)
        else:
            raise Exception("Only \"time_order\" equal to 2 is supported. Provided \"time_order\": " + str(time_order))

        if self.settings["use_oss"].GetBool():
            self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.OSS_SWITCH] = 1

        self.solver = self._get_solution_strategy()
        self.solver.SetEchoLevel(self.settings["echo_level"].GetInt())
        self.solver.Initialize()

        KratosMultiphysics.Logger.PrintInfo("::[NavierStokesExplicitSolver]:: ","Explicit fluid solver initialization finished.")

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
