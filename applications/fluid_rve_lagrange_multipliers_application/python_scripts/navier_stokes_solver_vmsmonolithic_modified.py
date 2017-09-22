from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.FluidRveLagrangeMultipliersApplication as KratosRVEapp

## Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Import base class file
import navier_stokes_solver_vmsmonolithic

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesSolver_VMSMonolithic_Modified(main_model_part, custom_settings)


BaseAlgorithm = navier_stokes_solver_vmsmonolithic.NavierStokesSolver_VMSMonolithic


class NavierStokesSolver_VMSMonolithic_Modified(BaseAlgorithm):

    def __init__(self, main_model_part, custom_settings):
        super(NavierStokesSolver_VMSMonolithic_Modified,self).__init__(main_model_part, custom_settings)


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

        builder_and_solver = KratosRVEapp.ResidualBasedEliminationBuilderAndSolverPeriodicNormalOnly(self.linear_solver,
                                                                                KratosCFD.PATCH_INDEX)



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

