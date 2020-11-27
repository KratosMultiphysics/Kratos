# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication

from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_solver_fractionalstep import NavierStokesSolverFractionalStep

from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
from KratosMultiphysics.FluidDynamicsApplication import check_and_prepare_model_process_fluid

def CreateSolver(model, custom_settings):
    return NavierStokesExplicitSolver(model, custom_settings)

class NavierStokesExplicitSolver(NavierStokesSolverFractionalStep):
    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass = True
        super().__init__(model,custom_settings)

        if custom_settings["formulation"]["element_type"].GetString() != "QSNavierStokesExplicit":
            raise Exception("NavierStokesExplicitSolver only accepts QSNavierStokesExplicit as the \"element_type\" in \"formulation\"")

        self.min_buffer_size = 3

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesExplicitSolver finished.")

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "FractionalStep",
            "model_part_name": "",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "",
                "reorder": false
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "predictor_corrector": false,
            "maximum_velocity_iterations": 3,
            "maximum_pressure_iterations": 3,
            "velocity_tolerance": 1e-3,
            "pressure_tolerance": 1e-2,
            "dynamic_tau": 0.01,
            "oss_switch": 0,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "compute_reactions": false,
            "reform_dofs_at_each_step" : false,
            "pressure_linear_solver_settings":  {
                "solver_type"                    : "amgcl",
                "max_iteration"                  : 200,
                "tolerance"                      : 1e-6,
                "provide_coordinates"            : false,
                "smoother_type"                  : "ilu0",
                "krylov_type"                    : "cg",
                "gmres_krylov_space_dimension"   : 100,
                "use_block_matrices_if_possible" : false,
                "coarsening_type"                : "aggregation",
                "scaling"                        : true,
                "verbosity"                      : 0
            },
            "velocity_linear_solver_settings": {
                "solver_type"                    : "amgcl",
                "max_iteration"                  : 200,
                "tolerance"                      : 1e-6,
                "provide_coordinates"            : false,
                "smoother_type"                  : "ilu0",
                "krylov_type"                    : "lgmres",
                "gmres_krylov_space_dimension"   : 100,
                "use_block_matrices_if_possible" : false,
                "coarsening_type"                : "aggregation",
                "scaling"                        : true,
                "verbosity"                      : 0
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "assign_neighbour_elements_to_conditions": false,
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step": false,
                "CFL_number": 1.0,
                "minimum_delta_time": 1.0e-8,
                "maximum_delta_time": 1.0e-2
            },
            "move_mesh_flag": false,
            "use_slip_conditions": true,
            "formulation": {
                "element_type": "QSNavierStokesExplicit",
                "condition_type": "WallCondition"
            },
            "time_integration_method": "explicit"
        }""")

        default_settings.AddMissingParameters(super(NavierStokesExplicitSolver, cls).GetDefaultParameters())
        return default_settings

    def _CreateLinearSolver(self):
        # Create the pressure linear solver
        pressure_linear_solver_configuration = self.settings["pressure_linear_solver_settings"]
        pressure_linear_solver = linear_solver_factory.ConstructSolver(pressure_linear_solver_configuration)
        # Create the velocity explicit solver
        velocity_linear_solver_configuration = self.settings["velocity_linear_solver_settings"]
        velocity_linear_solver = linear_solver_factory.ConstructSolver(velocity_linear_solver_configuration)
        # Return a tuple containing both linear solvers
        return (pressure_linear_solver, velocity_linear_solver)

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        domain_size = computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        # Create the pressure and velocity solvers
        # Note that linear_solvers is a tuple. The first item is the pressure
        # linear solver. The second item is the velocity explicit solver.
        linear_solvers = self._GetLinearSolver()

        # Create the explicit fractional step settings instance
        # TODO: next part would be much cleaner if we passed directly the parameters to the c++
        if self.settings["consider_periodic_conditions"].GetBool():
            pass
            # fractional_step_settings = FluidDynamicsApplication.FractionalStepSettingsPeriodic(
            #     computing_model_part,
            #     domain_size,
            #     self.settings["time_order"].GetInt(),
            #     self.settings["use_slip_conditions"].GetBool(),
            #     self.settings["move_mesh_flag"].GetBool(),
            #     self.settings["reform_dofs_at_each_step"].GetBool(),
            #     FluidDynamicsApplication.PATCH_INDEX)
        else:
            time_order=4 # Runge-Kutta order 4
            rebuild_level = 0
            fractional_step_settings = FluidDynamicsApplication.FractionalStepSettings(
                computing_model_part,
                domain_size,
                time_order,
                self.settings["use_slip_conditions"].GetBool(),
                self.settings["move_mesh_flag"].GetBool(),
                self.settings["reform_dofs_at_each_step"].GetBool(),
                rebuild_level)

        # Set the strategy echo level
        fractional_step_settings.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Set the velocity and pressure fractional step strategy settings
        fractional_step_settings.SetStrategy(FluidDynamicsApplication.StrategyLabel.Pressure,
            linear_solvers[0],
            self.settings["pressure_tolerance"].GetDouble(),
            self.settings["maximum_pressure_iterations"].GetInt())

        fractional_step_settings.SetStrategy(FluidDynamicsApplication.StrategyLabel.Velocity,
            linear_solvers[1],
            self.settings["velocity_tolerance"].GetDouble(),
            self.settings["maximum_velocity_iterations"].GetInt())

        # Create the fractional step strategy
        if self.settings["consider_periodic_conditions"].GetBool() == True:
            pass
            # solution_strategy = FluidDynamicsApplication.FractionalStepStrategy(
            #     computing_model_part,
            #     fractional_step_settings,
            #     self.settings["predictor_corrector"].GetBool(),
            #     self.settings["compute_reactions"].GetBool(),
            #     FluidDynamicsApplication.PATCH_INDEX)
        else:
            solution_strategy = FluidDynamicsApplication.ExplicitFractionalStepStrategy(
                computing_model_part,
                fractional_step_settings,
                self.settings["predictor_corrector"].GetBool(),
                self.settings["compute_reactions"].GetBool())

        return solution_strategy