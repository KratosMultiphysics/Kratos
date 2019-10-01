from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ChimeraApplication as KratosChimera
# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_solver_vmsmonolithic import NavierStokesSolverMonolithic

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesSolverMonolithicChimera(main_model_part, custom_settings)

class NavierStokesSolverMonolithicChimera(NavierStokesSolverMonolithic):
    def __init__(self, model, custom_settings):
        super(NavierStokesSolverMonolithicChimera,self).__init__(model,custom_settings)

        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithicChimera", "Construction of NavierStokesSolverMonolithic finished.")

    def AddVariables(self):
        super(NavierStokesSolverMonolithicChimera,self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)

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

       # Creating the time integration scheme
        if (self.element_integrates_in_time):
            # "Fake" scheme for those cases in where the element manages the time integration
            # It is required to perform the nodal update once the current time step is solved
            self.time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(
                self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]+1)
            # In case the BDF2 scheme is used inside the element, the BDF time discretization utility is required to update the BDF coefficients
            if (self.settings["time_scheme"].GetString() == "bdf2"):
                time_order = 2
                self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)
            else:
                err_msg = "Requested elemental time scheme \"" + self.settings["time_scheme"].GetString()+ "\" is not available.\n"
                err_msg += "Available options are: \"bdf2\""
                raise Exception(err_msg)
        else:
            if not hasattr(self, "_turbulence_model_solver"):
                # Bossak time integration scheme
                if self.settings["time_scheme"].GetString() == "bossak":
                    if self.settings["consider_periodic_conditions"].GetBool() == True:
                        self.time_scheme = KratosCFD.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
                            self.settings["alpha"].GetDouble(),
                            self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                            KratosCFD.PATCH_INDEX)
                    else:
                        self.time_scheme = KratosCFD.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
                            self.settings["alpha"].GetDouble(),
                            self.settings["move_mesh_strategy"].GetInt(),
                            self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
                # BDF2 time integration scheme
                elif self.settings["time_scheme"].GetString() == "bdf2":
                    self.time_scheme = KratosCFD.GearScheme()
                # Time scheme for steady state fluid solver
                elif self.settings["time_scheme"].GetString() == "steady":
                    self.time_scheme = KratosCFD.ResidualBasedSimpleSteadyScheme(
                            self.settings["velocity_relaxation"].GetDouble(),
                            self.settings["pressure_relaxation"].GetDouble(),
                            self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
                else:
                    err_msg = "Requested time scheme " + self.settings["time_scheme"].GetString() + " is not available.\n"
                    err_msg += "Available options are: \"bossak\", \"bdf2\" and \"steady\""
                    raise Exception(err_msg)
            else:
                KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithicForChimera turbulent solver is not possible.")
                raise NotImplementedError

        if self.settings["consider_periodic_conditions"].GetBool() == True:
            KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithicForChimera Periodic conditions are not implemented in this case .")
            raise NotImplementedError
        else:
            builder_and_solver = KratosChimera.ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera(self.linear_solver)

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

        self.solver.Check()

        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithicChimera", "Solver initialization finished.")
