from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

## Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Import base class file
import navier_stokes_solver_vmsmonolithic

def CreateSolver(model, custom_settings):
    return SteadyNavierStokesSolver_VMSMonolithic(model, custom_settings)

class SteadyNavierStokesSolver_VMSMonolithic(navier_stokes_solver_vmsmonolithic.NavierStokesSolverMonolithic):

    def __init__(self, model, custom_settings):

        # parse and strip parameters that do not exist in base class. we need to remove
        # extra parameters so base class doesn't throw an error. alternatively a single solver script
        # could be used and the scheme type could be passed in json parameters.
        self.velocity_relaxation_factor = custom_settings["velocity_relaxation_factor"].GetDouble()
        self.pressure_relaxation_factor = custom_settings["pressure_relaxation_factor"].GetDouble()
        base_settings = custom_settings
        base_settings.RemoveValue("velocity_relaxation_factor")
        base_settings.RemoveValue("pressure_relaxation_factor")

        # call base class constructor with remaining parameters
        super().__init__(model, base_settings)

        if self.settings["consider_periodic_conditions"].GetBool() == True:
            raise ValueError("consider_periodic_conditions not supported yet.")

        if self.settings["move_mesh_strategy"].GetInt() != 0:
            raise ValueError("move_mesh_strategy not supported yet.")

    def Initialize(self):

        self.computing_model_part = self.GetComputingModelPart()

        MoveMeshFlag = False

        # TODO: TURBULENCE MODELS ARE NOT ADDED YET
        #~ # Turbulence model
        #~ if self.use_spalart_allmaras:
            #~ self.activate_spalart_allmaras()

        # creating the solution strategy
        self.conv_criteria = KratosCFD.VelPrCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                     self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                     self.settings["relative_pressure_tolerance"].GetDouble(),
                                                     self.settings["absolute_pressure_tolerance"].GetDouble())

        (self.conv_criteria).SetEchoLevel(self.settings["echo_level"].GetInt())

        self.time_scheme = KratosCFD.ResidualBasedSimpleSteadyScheme(self.velocity_relaxation_factor, self.pressure_relaxation_factor, self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        # TODO: TURBULENCE MODELS ARE NOT ADDED YET
        #~ if self.turbulence_model is None:
            #~ if self.periodic == True:
                #~ self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent\
                #~ (self.alpha, self.move_mesh_strategy, self.domain_size, PATCH_INDEX)
            #~ else:
                #~ self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent\
                #~ (self.alpha, self.move_mesh_strategy, self.domain_size)
        #~ else:
            #~ self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent\
                #~ (self.alpha, self.move_mesh_strategy, self.domain_size, self.turbulence_model)

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

        if self.settings["stabilization"].Has("dynamic_tau"):
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["stabilization"]["dynamic_tau"].GetDouble())
        if self.settings["stabilization"].Has("oss_switch"):
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH, self.settings["stabilization"]["oss_switch"].GetInt())

        print ("Monolithic solver initialization finished.")
