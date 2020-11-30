# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication

from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_solver_fractionalstep import NavierStokesSolverFractionalStep

from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
from KratosMultiphysics.FluidDynamicsApplication import check_and_prepare_model_process_fluid

def CreateSolver(model, custom_settings):
    return NavierStokesSolverFractionalStepSemiExplicit(model, custom_settings)

class NavierStokesSolverFractionalStepSemiExplicit(NavierStokesSolverFractionalStep):
    def __init__(self, model, custom_settings):
        super().__init__(model,custom_settings)

        if custom_settings["formulation"]["element_type"].GetString() != "QSNavierStokesSemiExplicit":
            raise Exception("NavierStokesSolverFractionalStepSemiExplicit only accepts QSNavierStokesSemiExplicit as the \"element_type\" in \"formulation\"")
        if custom_settings["domain_size"].GetInt() == 2:
            self.condition_name = "LineCondition"
        elif custom_settings["domain_size"].GetInt() == 3:
            self.condition_name = "SurfaceCondition"
        else:
            err_mgs = "Wrong domain size "
            raise Exception(err_msg)

        self.min_buffer_size = 2

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesSolverFractionalStepSemiExplicit finished.")

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
                "element_type": "QSNavierStokesSemiExplicit",
                "condition_type": ""
            },
            "time_integration_method": "semiexplicit"
        }""")

        default_settings.AddMissingParameters(super(NavierStokesSolverFractionalStepSemiExplicit, cls).GetDefaultParameters())
        return default_settings

    def AddDofs(self):
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_X, KratosMultiphysics.REACTION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.REACTION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE,self.main_model_part)
        if domain_size == 3:
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Z, KratosMultiphysics.REACTION_Z,self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver DOFs added correctly.")

    def _CreateLinearSolver(self):
        pressure_linear_solver_configuration = self.settings["pressure_linear_solver_settings"]
        pressure_linear_solver = linear_solver_factory.ConstructSolver(pressure_linear_solver_configuration)
        return pressure_linear_solver

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        domain_size = computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        # Create the pressure solver
        linear_solvers = self._GetLinearSolver()

        time_order= 4 # Runge-Kutta order 4
        fractional_step_settings = FluidDynamicsApplication.FractionalStepSettings(
            computing_model_part,
            domain_size,
            time_order,
            self.settings["use_slip_conditions"].GetBool(),
            self.settings["move_mesh_flag"].GetBool(),
            self.settings["reform_dofs_at_each_step"].GetBool())

        # Set the strategy echo level
        fractional_step_settings.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Set the velocity and pressure fractional step strategy settings
        fractional_step_settings.SetStrategy(FluidDynamicsApplication.StrategyLabel.Pressure,
            linear_solvers,
            self.settings["pressure_tolerance"].GetDouble(),
            self.settings["maximum_pressure_iterations"].GetInt())

        fractional_step_settings.SetExplicitStrategy(FluidDynamicsApplication.StrategyLabel.Velocity)

        solution_strategy = FluidDynamicsApplication.SemiExplicitFractionalStepStrategy(
            computing_model_part,
            fractional_step_settings,
            self.settings["predictor_corrector"].GetBool(),
            self.settings["compute_reactions"].GetBool())

        return solution_strategy