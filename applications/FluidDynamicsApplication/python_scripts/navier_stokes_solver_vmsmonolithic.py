from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
import navier_stokes_base_solver

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesSolver_VMSMonolithic(main_model_part, custom_settings)

class NavierStokesSolver_VMSMonolithic(navier_stokes_base_solver.NavierStokesBaseSolver):

    def __init__(self, main_model_part, custom_settings):

        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "navier_stokes_solver_vmsmonolithic",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "maximum_iterations": 10,
            "dynamic_tau": 0.01,
            "oss_switch": 0,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "compute_reactions": false,
            "reform_dofs_at_each_step": true,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"        : {
                    "solver_type" : "AMGCL",
                    "smoother_type":"ilu0",
                    "krylov_type": "gmres",
                    "coarsening_type": "aggregation",
                    "max_iteration": 200,
                    "provide_coordinates": false,
                    "gmres_krylov_space_dimension": 100,
                    "verbosity" : 0,
                    "tolerance": 1e-7,
                    "scaling": false,
                    "block_size": 1,
                    "use_block_matrices_if_possible" : true,
                    "coarse_enough" : 5000
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : true,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01
            },
            "alpha":-0.3,
            "move_mesh_strategy": 0,
            "periodic": "periodic",
            "move_mesh_flag": false,
            "turbulence_model": "None",
            "reorder": false
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        ## Construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        ## Set the element replace settings
        self.settings.AddEmptyValue("element_replace_settings")
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
            self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                {
                "element_name":"VMS3D4N",
                "condition_name": "MonolithicWallCondition3D"
                }
                """)
        elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                {
                "element_name":"VMS2D3N",
                "condition_name": "MonolithicWallCondition2D"
                }
                """)
        else:
            raise Exception("domain size is not 2 or 3")

        print("Construction of NavierStokesSolver_VMSMonolithic finished")


    def AddVariables(self):
        ## Add base class variables
        super(NavierStokesSolver_VMSMonolithic, self).AddVariables()
        ## Add specific variables needed for the monolithic solver
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)         # Kinematic viscosity value
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)

        print("Monolithic fluid solver variables added correctly")


    def Initialize(self):

        self.computing_model_part = self.GetComputingModelPart()

        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        # Creating the solution strategy
        self.conv_criteria = KratosCFD.VelPrCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                     self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                     self.settings["relative_pressure_tolerance"].GetDouble(),
                                                     self.settings["absolute_pressure_tolerance"].GetDouble())

        (self.conv_criteria).SetEchoLevel(self.settings["echo_level"].GetInt())

        if (self.settings["turbulence_model"].GetString() == "None"):
            if self.settings["consider_periodic_conditions"].GetBool() == True:
                self.time_scheme = KratosCFD.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(self.settings["alpha"].GetDouble(),
                                                                                                          self.settings["move_mesh_strategy"].GetInt(),
                                                                                                          self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                                                          KratosCFD.PATCH_INDEX)
            else:
                self.time_scheme = KratosCFD.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(self.settings["alpha"].GetDouble(),
                                                                                                          self.settings["move_mesh_strategy"].GetInt(),
                                                                                                          self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
        else:
            raise Exception("Turbulence models are not added yet.")

        if self.settings["consider_periodic_conditions"].GetBool() == True:
            builder_and_solver = KratosCFD.ResidualBasedBlockBuilderAndSolverPeriodic(self.linear_solver,
                                                                                KratosCFD.PATCH_INDEX)
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)


        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.main_model_part,
                                                                            self.time_scheme,
                                                                            self.linear_solver,
                                                                            self.conv_criteria,
                                                                            builder_and_solver,
                                                                            self.settings["maximum_iterations"].GetInt(),
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            self.settings["move_mesh_flag"].GetBool())

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())
        (self.solver).Check()

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH, self.settings["oss_switch"].GetInt())

        (self.solver).Initialize()

        print ("Monolithic solver initialization finished.")


    def InitializeSolutionStep(self):
        (self.solver).InitializeSolutionStep()


    def Solve(self):
        # Note that the first time step in buffer is filled when reading the model part
        (self.solver).Solve()


    def _ExecuteAfterReading(self):
        super(NavierStokesSolver_VMSMonolithic, self)._ExecuteAfterReading()

        # Read the KINEMATIC VISCOSITY
        for el in self.main_model_part.Elements:
            # rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            kin_viscosity = el.Properties.GetValue(KratosMultiphysics.VISCOSITY)
            break

        # KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.VISCOSITY, kin_viscosity, self.main_model_part.Nodes)
