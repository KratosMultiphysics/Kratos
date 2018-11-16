from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ChimeraApplication as KratosChimera

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from navier_stokes_solver_vmsmonolithic import NavierStokesSolverMonolithic

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesSolverMonolithicChimera(main_model_part, custom_settings)

class NavierStokesSolverMonolithicChimera(NavierStokesSolverMonolithic):
    def __init__(self, model, custom_settings):
        if custom_settings.Has("implementation"):
            self.implementation = custom_settings["implementation"].GetString()
            custom_settings.RemoveValue("implementation")
        else:
            raise Exception('No "implementation" specified!')

        super(NavierStokesSolverMonolithicChimera,self).__init__(model,custom_settings)

        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithicChimera", "Construction of NavierStokesSolverMonolithic finished.")


    def AddVariables(self):
        super(NavierStokesSolverMonolithicChimera,self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.PRESSURE_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.SHEAR_FORCE)

        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithicChimera", "Fluid solver variables added correctly.")

    def Initialize(self):

        #self.computing_model_part = self.GetComputingModelPart()
        self.computing_model_part = self.main_model_part
        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._get_automatic_time_stepping_utility()

        # Creating the solution strategy
        self.conv_criteria = KratosCFD.VelPrCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                     self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                     self.settings["relative_pressure_tolerance"].GetDouble(),
                                                     self.settings["absolute_pressure_tolerance"].GetDouble())

        (self.conv_criteria).SetEchoLevel(self.settings["echo_level"].GetInt())

        if (self.settings["turbulence_model"].GetString() == "None"):
            if self.settings["time_scheme"].GetString() == "bossak":
                if self.settings["consider_periodic_conditions"].GetBool() == True:
                    self.time_scheme = KratosCFD.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
                                        self.settings["alpha"].GetDouble(),
                                        self.settings["move_mesh_strategy"].GetInt(),
                                        self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                        KratosCFD.PATCH_INDEX)
                else:
                    self.time_scheme = KratosCFD.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
                                        self.settings["alpha"].GetDouble(),
                                        self.settings["move_mesh_strategy"].GetInt(),
                                        self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
            elif self.settings["time_scheme"].GetString() == "bdf2":
                self.time_scheme = KratosCFD.GearScheme()
            elif self.settings["time_scheme"].GetString() == "steady":
                self.time_scheme = KratosCFD.ResidualBasedSimpleSteadyScheme(
                                        self.settings["velocity_relaxation"].GetDouble(),
                                        self.settings["pressure_relaxation"].GetDouble(),
                                        self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
        else:
            raise Exception("Turbulence models are not added yet.")

        if self.settings["consider_periodic_conditions"].GetBool() == True:
            builder_and_solver = KratosCFD.ResidualBasedBlockBuilderAndSolverPeriodic(self.linear_solver,
                                                                                KratosCFD.PATCH_INDEX)
        else:
            if (self.implementation == "MPC"):
                builder_and_solver = KratosChimera.ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera(self.linear_solver)
                #builder_and_solver = KratosChimera.ResidualBasedBlockBuilderAndSolverWithMpcChimera(self.linear_solver)
            else:
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.computing_model_part,
                                                                            self.time_scheme,
                                                                            self.linear_solver,
                                                                            self.conv_criteria,
                                                                            builder_and_solver,
                                                                            self.settings["maximum_iterations"].GetInt(),
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            self.settings["move_mesh_flag"].GetBool())

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        self.stabilization.SetProcessInfo(self.computing_model_part)

        (self.solver).Initialize()

        self.solver.Check()

        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithicChimera", "Solver initialization finished.")
