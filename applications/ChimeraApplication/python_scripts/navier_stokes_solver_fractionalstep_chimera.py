from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ChimeraApplication as KratosChimera
from KratosMultiphysics.ChimeraApplication import chimera_setup_utils

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_solver_fractionalstep import NavierStokesSolverFractionalStep


def CreateSolver(model, custom_settings):
    return NavierStokesSolverFractionalStepForChimera(model, custom_settings)

class NavierStokesSolverFractionalStepForChimera(NavierStokesSolverFractionalStep):


    def __init__(self, model, custom_settings):
        [self.chimera_settings, self.chimera_internal_parts, custom_settings] = chimera_setup_utils.SeparateAndValidateChimeraSettings(custom_settings)
        super(NavierStokesSolverFractionalStepForChimera,self).__init__(model,custom_settings)
        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverFractionalStepForChimera", "Construction of NavierStokesSolverFractionalStepForChimera finished.")

    def AddVariables(self):
        super(NavierStokesSolverFractionalStepForChimera,self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.CHIMERA_DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATIONAL_ANGLE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATIONAL_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATION_MESH_DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosChimera.ROTATION_MESH_VELOCITY)


        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverFractionalStepForChimera", "Fluid solver variables added correctly.")


    def ImportModelPart(self):
        if(self.settings["model_import_settings"]["input_type"].GetString() == "chimera"):
            chimera_mp_import_settings = []
            for chimera_part_levels in self.chimera_settings["chimera_parts"]:
                for chimera_part_settings in chimera_part_levels:
                    chimera_mp_import_settings.append( chimera_part_settings["model_import_settings"].Clone() )

            material_file_name = self.settings["material_import_settings"]["materials_filename"].GetString()
            import KratosMultiphysics.ChimeraApplication.chimera_modelpart_import as chim_mp_imp
            chim_mp_imp.ImportChimeraModelparts(self.main_model_part, chimera_mp_import_settings, material_file=material_file_name, parallel_type="OpenMP")
            KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverFractionalStepForChimera", " Import of all chimera modelparts completed.")
        else:# we can use the default implementation in the base class
            super(NavierStokesSolverFractionalStepForChimera,self).ImportModelPart()

    def Initialize(self):
        self.chimera_process = chimera_setup_utils.GetApplyChimeraProcess(self.model, self.chimera_settings, self.settings)
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

        chimera_setup_utils.SetChimeraInternalPartsFlag(self.model, self.chimera_internal_parts)


    def InitializeSolutionStep(self):
        self.chimera_process.ExecuteInitializeSolutionStep()
        super(NavierStokesSolverFractionalStepForChimera,self).InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(NavierStokesSolverFractionalStepForChimera,self).FinalizeSolutionStep()
        ## Depending on the setting this will clear the created constraints
        self.chimera_process.ExecuteFinalizeSolutionStep()