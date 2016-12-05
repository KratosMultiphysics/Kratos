from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Import base class file
import navier_stokes_base_solver

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesSolver_FractionalStep(main_model_part, custom_settings)

class NavierStokesSolver_FractionalStep(navier_stokes_base_solver.NavierStokesBaseSolver):

    def __init__(self, main_model_part, custom_settings):

        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "navier_stokes_solver_fractionalstep",
            "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
            },
            "predictor_corrector": false,
            "maximum_velocity_iterations": 3,
            "maximum_pressure_iterations": 3,
            "velocity_tolerance": 1e-3,
            "pressure_tolerance": 1e-2,
            "echo_level": 1,
            "consider_periodic_conditions": false,
            "time_order": 2,
            "dynamic_tau": 0.001,
            "compute_reactions": false,
            "divergence_clearance_steps": 0,
            "reform_dofs_at_each_step": false,
            "pressure_linear_solver_settings":  {
                "solver_type"                    : "AMGCL",
                "max_iteration"                  : 200,
                "tolerance"                      : 1e-6,
                "provide_coordinates"            : false,
                "smoother_type"                  : "damped_jacobi",
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
                "krylov_type"                    : "bicgstab",
                "gmres_krylov_space_dimension"   : 100,
                "use_block_matrices_if_possible" : false,
                "coarsening_type"                : "aggregation",
                "scaling"                        : true,
                "verbosity"                      : 0
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts":[""],
            "no_skin_parts":[""]
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        ## Construct the linear solvers
        import linear_solver_factory
        self.pressure_linear_solver = linear_solver_factory.ConstructSolver(self.settings["pressure_linear_solver_settings"])
        self.velocity_linear_solver = linear_solver_factory.ConstructSolver(self.settings["velocity_linear_solver_settings"])

        self.compute_reactions = self.settings["compute_reactions"].GetBool()

        ## Set the element replace settings
        self.settings.AddEmptyValue("element_replace_settings")
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
            self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                {
                "element_name":"FractionalStep3D4N",
                "condition_name": "WallCondition3D3N"
                }
                """)
        elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                {
                "element_name":"FractionalStep2D3N",
                "condition_name": "WallCondition2D2N"
                }
                """)
        else:
            raise Exception("Domain size is not 2 or 3.")

        print("Construction of NavierStokesSolver_FractionalStep finished.")


    def AddVariables(self):
        ## Add base class variables
        super(NavierStokesSolver_FractionalStep, self).AddVariables()
        ## Add specific variables needed for the fractional step solver
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FRACT_VEL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE_OLD_IT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESS_PROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONV_PROJ)

        print("Fractional step fluid solver variables added correctly")


    def Initialize(self):
        compute_model_part = self.GetComputingModelPart()

        MoveMeshFlag = False

        self.use_slip_conditions = True

        #TODO: next part would be much cleaner if we passed directly the parameters to the c++
        if self.settings["consider_periodic_conditions"] == True:
            self.solver_settings = KratosCFD.FractionalStepSettingsPeriodic(compute_model_part,
                                                                            compute_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                            self.settings.GetInt(),
                                                                            self.use_slip_conditions,
                                                                            MoveMeshFlag,
                                                                            self.settings["reform_dofs_at_each_step]"].GetBool(),
                                                                            KratosCFD.PATCH_INDEX)

        else:
            self.solver_settings = KratosCFD.FractionalStepSettings(compute_model_part,
                                                                    compute_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                    self.settings["time_order"].GetInt(),
                                                                    self.use_slip_conditions,
                                                                    MoveMeshFlag,
                                                                    self.settings["reform_dofs_at_each_step"].GetBool())

        self.solver_settings.SetEchoLevel(self.settings["echo_level"].GetInt())

        self.solver_settings.SetStrategy(KratosCFD.StrategyLabel.Velocity,
                                         self.velocity_linear_solver,
                                         self.settings["velocity_tolerance"].GetDouble(),
                                         self.settings["maximum_velocity_iterations"].GetInt())

        self.solver_settings.SetStrategy(KratosCFD.StrategyLabel.Pressure,
                                         self.pressure_linear_solver,
                                         self.settings["pressure_tolerance"].GetDouble(),
                                         self.settings["maximum_pressure_iterations"].GetInt())


        if self.settings["consider_periodic_conditions"].GetBool() == True:
            self.solver = KratosCFD.FSStrategy(compute_model_part,
                                               self.solver_settings,
                                               self.settings["predictor_corrector"].GetBool(),
                                               KratosCFD.PATCH_INDEX)
        else:
            self.solver = KratosCFD.FSStrategy(compute_model_part,
                                               self.solver_settings,
                                               self.settings["predictor_corrector"].GetBool())

        (self.solver).Check()

        print ("Initialization NavierStokesSolver_FractionalStep finished.")


    def SolverFinalizeSolutionStep(self):
        (self.solver).FinalizeSolutionStep()
        if(self.compute_reactions):
                    (self.solver).CalculateReactions()


    def Solve(self):
        (self.solver).Solve()
        if(self.compute_reactions):
            (self.solver).CalculateReactions()
