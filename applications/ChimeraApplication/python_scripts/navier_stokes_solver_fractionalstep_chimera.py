from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ChimeraApplication as KratosChimera

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_solver_fractionalstep import NavierStokesSolverFractionalStep


def CreateSolver(model, custom_settings):
    return NavierStokesSolverFractionalStepForChimera(model, custom_settings)

class NavierStokesSolverFractionalStepForChimera(NavierStokesSolverFractionalStep):


    def __init__(self, model, custom_settings):
        super(NavierStokesSolverFractionalStepForChimera,self).__init__(model,custom_settings)
        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverFractionalStepForChimera", "Construction of NavierStokesSolverFractionalStepForChimera finished.")

    def AddVariables(self):
        super(NavierStokesSolverFractionalStepForChimera,self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)

        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverFractionalStepForChimera", "Fluid solver variables added correctly.")


    def Initialize(self):
        #self.computing_model_part = self.GetComputingModelPart()
        self.computing_model_part =self.main_model_part
        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._get_automatic_time_stepping_utility()

        #TODO: next part would be much cleaner if we passed directly the parameters to the c++
        if self.settings["consider_periodic_conditions"] == True:
            KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverFractionalStepForChimera Periodic conditions are not implemented in this case .")
            raise NotImplementedError
        else:
            self.solver_settings = KratosChimera.FractionalStepSettings(self.computing_model_part,
                                                                    self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                    self.settings["time_order"].GetInt(),
                                                                    self.settings["use_slip_conditions"].GetBool(),
                                                                    self.settings["move_mesh_flag"].GetBool(),
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
