from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication","KratosPFEM2Application")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.KratosPFEM2Application as KratosPFEM2

# Import base class file
from navier_stokes_solver_vmsmonolithic import NavierStokesSolverMonolithic

def CreateSolver(model, custom_settings):
    return PFEM2NavierStokesMonolithicSolver(model, custom_settings)

class PFEM2NavierStokesMonolithicSolver(NavierStokesSolverMonolithic):

    def __init__(self, model, custom_settings):
        super(PFEM2NavierStokesMonolithicSolver,self).__init__(model,custom_settings)

    def AddVariables(self):
        # FluidDynamics Variables
        super(PFEM2NavierStokesMonolithicSolver,self).AddVariables()

        # PFEM2 Variables
        model_part.AddNodalSolutionStepVariable(YP)
        model_part.AddNodalSolutionStepVariable(PRESS_PROJ); # Tihs variable is required but unused by ModeParticleUtilityPFEM2
        model_part.AddNodalSolutionStepVariable(PROJECTED_VELOCITY)
        model_part.AddNodalSolutionStepVariable(DELTA_VELOCITY)
        model_part.AddNodalSolutionStepVariable(MEAN_SIZE) # Tihs variable is required but unused by ModeParticleUtilityPFEM2

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("PFEM2NavierStokesMonolithicSolver", "Fluid solver variables added correctly.")

    def AddDofs(self):
        super(PFEM2NavierStokesMonolithicSolver,self).AddDofs()

        # PFEM2 Dofs # mms: we don't need DISTANCE dof (it is the particles pseudo-levelset)
        # KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISTANCE, self.main_model_part)

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("PFEM2NavierStokesMonolithicSolver", "Fluid solver DOFs added correctly.")

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
            if self.settings["time_scheme"].GetString() == "bossak":
                if self.settings["consider_periodic_conditions"].GetBool() == True:
                    self.time_scheme = KratosPFEM2.ResidualBasedPredictorCorrectorVelocityBossakAleScheme(
                                        self.settings["alpha"].GetDouble(),
                                        self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                        KratosCFD.PATCH_INDEX)
                else:
                    self.time_scheme = KratosPFEM2.ResidualBasedPredictorCorrectorVelocityBossakAleScheme(
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

        self.formulation.SetProcessInfo(self.computing_model_part)

        (self.solver).Initialize()

        KratosMultiphysics.Logger.PrintInfo("PFEM2NavierStokesMonolithicSolver", "Solver initialization finished.")
