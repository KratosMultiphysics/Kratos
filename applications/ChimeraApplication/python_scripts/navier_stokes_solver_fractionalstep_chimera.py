from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ChimeraApplication as KratosChimera

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from fluid_solver import FluidSolver

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesSolverFractionalStepForChimera(main_model_part, custom_settings)

class NavierStokesSolverFractionalStepForChimera(FluidSolver):

    def _ValidateSettings(self, settings):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "FractionalStep",
            "model_part_name" :"",
            "domain_size" :-1,
            "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
            },
            "implementation"    : "MPC",
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
            "reform_dofs_at_each_step": true,
            "pressure_linear_solver_settings":  {
                "solver_type"                    : "AMGCL",
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
                "solver_type"                    : "AMGCL",
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
                "automatic_time_step" : true,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01
            },
            "move_mesh_flag": false,
            "use_slip_conditions": true,
            "reorder": false
        }""")

        ## Overwrite the default settings with user-provided parameters
        settings.ValidateAndAssignDefaults(default_settings)
        return settings

    def __init__(self, main_model_part, custom_settings):
        super(NavierStokesSolverFractionalStepForChimera,self).__init__(main_model_part,custom_settings)
        self.element_name = "FractionalStep"
        self.condition_name = "WallCondition"
        self.min_buffer_size = 3

        self._is_printing_rank = True

        ## Construct the linear solvers
        import linear_solver_factory
        self.pressure_linear_solver = linear_solver_factory.ConstructSolver(self.settings["pressure_linear_solver_settings"])
        self.velocity_linear_solver = linear_solver_factory.ConstructSolver(self.settings["velocity_linear_solver_settings"])

        self.compute_reactions = self.settings["compute_reactions"].GetBool()

        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverFractionalStepForChimera", "Construction of NavierStokesSolverFractionalStepForChimera solver finished.")


        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        #self.main_model_part = main_model_part

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
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

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)

        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.PRESSURE_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.SHEAR_FORCE)


        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.Q_VALUE)
        if self.settings["consider_periodic_conditions"].GetBool() == True:
            self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.PATCH_INDEX)

        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverFractionalStepForChimera", "Fluid solver variables added correctly.")

    def PrepareModelPart(self):
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self._set_physical_properties()
        super(NavierStokesSolverFractionalStepForChimera, self).PrepareModelPart()

    def Initialize(self):
        #self.computing_model_part = self.GetComputingModelPart()
        self.computing_model_part =self.main_model_part
        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._get_automatic_time_stepping_utility()

        #TODO: next part would be much cleaner if we passed directly the parameters to the c++
        if self.settings["consider_periodic_conditions"] == True:
            self.solver_settings = KratosCFD.FractionalStepSettingsPeriodic(self.computing_model_part,
                                                                            self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                            self.settings.GetInt(),
                                                                            self.settings["use_slip_conditions"].GetBool(),
                                                                            self.settings["move_mesh_flag"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            KratosCFD.PATCH_INDEX)

        else:
            if (self.settings["implementation"].GetString() == "MPC"):
                print("fractional step chimera : MPC ")
                self.solver_settings = KratosChimera.FractionalStepSettingsForChimera(self.computing_model_part,
                                                                        self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                        self.settings["time_order"].GetInt(),
                                                                        self.settings["use_slip_conditions"].GetBool(),
                                                                        self.settings["move_mesh_flag"].GetBool(),
                                                                        self.settings["reform_dofs_at_each_step"].GetBool())
            else:
                self.solver_settings = KratosChimera.FractionalStepSettingsForChimera(self.computing_model_part,
                                                                        self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                        self.settings["time_order"].GetInt(),
                                                                        self.settings["use_slip_conditions"].GetBool(),
                                                                        self.settings["move_mesh_flag"].GetBool(),
                                                                        self.settings["reform_dofs_at_each_step"].GetBool())

        self.solver_settings.SetEchoLevel(self.settings["echo_level"].GetInt())
        
        builder_and_solver_velocity = KratosChimera.ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera(self.velocity_linear_solver)
        builder_and_solver_pressure = KratosChimera.ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera(self.pressure_linear_solver)
       
        self.solver_settings.SetStrategy(KratosChimera.ChimeraStrategyLabel.Velocity,
                                         self.velocity_linear_solver,
                                         self.settings["velocity_tolerance"].GetDouble(),
                                         self.settings["maximum_velocity_iterations"].GetInt(),
                                         builder_and_solver_velocity)

        self.solver_settings.SetStrategy(KratosChimera.ChimeraStrategyLabel.Pressure,
                                         self.pressure_linear_solver,
                                         self.settings["pressure_tolerance"].GetDouble(),
                                         self.settings["maximum_pressure_iterations"].GetInt(),
                                         builder_and_solver_pressure)


        if self.settings["consider_periodic_conditions"].GetBool() == True:
            self.solver = KratosCFD.FSStrategy(self.computing_model_part,
                                               self.solver_settings,
                                               self.settings["predictor_corrector"].GetBool(),
                                               KratosCFD.PATCH_INDEX)
        else:
            self.solver = KratosChimera.FSStrategyForChimera(self.computing_model_part,
                                               self.solver_settings,
                                               self.settings["predictor_corrector"].GetBool())

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH, self.settings["oss_switch"].GetInt())

        (self.solver).Initialize()

        self.solver.Check()

        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverFractionalStepForChimera", "Solver initialization finished.")


    def SolveSolutionStep(self):
        print("solve solution step")
        if self._TimeBufferIsInitialized():
            self.solver.SolveSolutionStep()
            if self.compute_reactions:
                self.solver.CalculateReactions()


    def _set_physical_properties(self):
        # Transfer density and (kinematic) viscostity to the nodes
        for el in self.main_model_part.Elements:
            rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            if rho <= 0.0:
                raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho,el.Properties.Id))
            dyn_viscosity = el.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
            if dyn_viscosity <= 0.0:
                raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(dyn_viscosity,el.Properties.Id))
            kin_viscosity = dyn_viscosity / rho
            break

        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.VISCOSITY, kin_viscosity, self.main_model_part.Nodes)
