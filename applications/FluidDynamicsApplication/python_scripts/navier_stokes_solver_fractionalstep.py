# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver

def CreateSolver(model, custom_settings):
    return NavierStokesSolverFractionalStep(model, custom_settings)

class NavierStokesSolverFractionalStep(FluidSolver):

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "fractional_step",
            "model_part_name": "",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
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
            "time_order": 2,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "assign_neighbour_elements_to_conditions": true,
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
            "skin_parts":[""],
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : false,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01
            },
            "move_mesh_flag": false,
            "use_slip_conditions": true,
            "formulation": {
                "element_type": "FractionalStep",
                "condition_type": "WallCondition"
            }
        }""")

        default_settings.AddMissingParameters(super(NavierStokesSolverFractionalStep, cls).GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):
        super(NavierStokesSolverFractionalStep,self).__init__(model,custom_settings)

        if custom_settings["formulation"]["element_type"].GetString() != "FractionalStep":
            raise Exception("NavierStokesFractionalStepSolver only accepts FractionalStep as the \"element_type\" in \"formulation\"")

        self.element_name = custom_settings["formulation"]["element_type"].GetString()
        self.condition_name = custom_settings["formulation"]["condition_type"].GetString()
        self.element_has_nodal_properties = True

        self.min_buffer_size = 3

        self.compute_reactions = self.settings["compute_reactions"].GetBool()

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH, self.settings["oss_switch"].GetInt())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesSolverFractionalStep solver finished.")

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.Y_WALL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FRACT_VEL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE_OLD_IT)
        # The following are used for the calculation of projections
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESS_PROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONV_PROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DIVPROJ)

        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.Q_VALUE)
        if self.settings["consider_periodic_conditions"].GetBool() == True:
            self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.PATCH_INDEX)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver variables added correctly.")

    def Initialize(self):
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        solution_strategy.Initialize()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def _CreateScheme(self):
        pass

    def _CreateLinearSolver(self):
        # Create the pressure linear solver
        pressure_linear_solver_configuration = self.settings["pressure_linear_solver_settings"]
        pressure_linear_solver = linear_solver_factory.ConstructSolver(pressure_linear_solver_configuration)
        # Create the velocity linear solver
        velocity_linear_solver_configuration = self.settings["velocity_linear_solver_settings"]
        velocity_linear_solver = linear_solver_factory.ConstructSolver(velocity_linear_solver_configuration)
        # Return a tuple containing both linear solvers
        return (pressure_linear_solver, velocity_linear_solver)

    def _CreateConvergenceCriterion(self):
        pass

    def _CreateBuilderAndSolver(self):
        pass

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        domain_size = computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        # Create the pressure and velocity linear solvers
        # Note that linear_solvers is a tuple. The first item is the pressure
        # linear solver. The second item is the velocity linear solver.
        linear_solvers = self._GetLinearSolver()

        # Create the fractional step settings instance
        # TODO: next part would be much cleaner if we passed directly the parameters to the c++
        if self.settings["consider_periodic_conditions"].GetBool():
            fractional_step_settings = KratosCFD.FractionalStepSettingsPeriodic(
                computing_model_part,
                domain_size,
                self.settings["time_order"].GetInt(),
                self.settings["use_slip_conditions"].GetBool(),
                self.settings["move_mesh_flag"].GetBool(),
                self.settings["reform_dofs_at_each_step"].GetBool(),
                KratosCFD.PATCH_INDEX)
        else:
            fractional_step_settings = KratosCFD.FractionalStepSettings(
                computing_model_part,
                domain_size,
                self.settings["time_order"].GetInt(),
                self.settings["use_slip_conditions"].GetBool(),
                self.settings["move_mesh_flag"].GetBool(),
                self.settings["reform_dofs_at_each_step"].GetBool())

        # Set the strategy echo level
        fractional_step_settings.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Set the velocity and pressure fractional step strategy settings
        fractional_step_settings.SetStrategy(KratosCFD.StrategyLabel.Pressure,
            linear_solvers[0],
            self.settings["pressure_tolerance"].GetDouble(),
            self.settings["maximum_pressure_iterations"].GetInt())

        fractional_step_settings.SetStrategy(KratosCFD.StrategyLabel.Velocity,
            linear_solvers[1],
            self.settings["velocity_tolerance"].GetDouble(),
            self.settings["maximum_velocity_iterations"].GetInt())

        # Create the fractional step strategy
        if self.settings["consider_periodic_conditions"].GetBool() == True:
            solution_strategy = KratosCFD.FractionalStepStrategy(
                computing_model_part,
                fractional_step_settings,
                self.settings["predictor_corrector"].GetBool(),
                self.settings["compute_reactions"].GetBool(),
                KratosCFD.PATCH_INDEX)
        else:
            solution_strategy = KratosCFD.FractionalStepStrategy(
                computing_model_part,
                fractional_step_settings,
                self.settings["predictor_corrector"].GetBool(),
                self.settings["compute_reactions"].GetBool())

        return solution_strategy

    def _SetNodalProperties(self):
        # Get density and dynamic viscostity from the properties of the first element
        for el in self.main_model_part.Elements:
            rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            if rho <= 0.0:
                raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho,el.Properties.Id))
            dyn_viscosity = el.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
            if dyn_viscosity <= 0.0:
                raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(dyn_viscosity,el.Properties.Id))
            kin_viscosity = dyn_viscosity / rho
            break
        else:
            raise Exception("No fluid elements found in the main model part.")
        # Transfer the obtained properties to the nodes
        KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.VISCOSITY, kin_viscosity, self.main_model_part.Nodes)
