"""
 --- apply_mass_local_conservation_process.py ---- Fri, May 22, 2020 11:27:27 AM ----
  Altair Manufacturing Solver

  Author: ddiez --- Maintained by: ddiez
  Copyright: Altair Engineering, Inc. 2015 - 2020
 ************************************************************
"""

from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FillingApplication as Filling
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.python_linear_solver_factory as python_linear_solver_factory

if KratosMultiphysics.DataCommunicator.GetDefault().IsDistributed():
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos
    import KratosMultiphysics.TrilinosApplication.trilinos_linear_solver_factory as trilinos_linear_solver_factory

def Factory(settings, Model):
    if( not isinstance(settings, KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyLocalMassConservationCheckProcess( Model, settings["Parameters"] )

class ApplyLocalMassConservationCheckProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "model_part_name"                    : "fluid_model_part",
            "mass_correction_settings"           : {
                "echo_level" : 1
            },
            "log_file_name"                      : "mass_conservation.log",
            "convector_settings"                 : {},
            "check_volume_loss_at_the_end"       : false,
            "min_dt_factor"                      : 1e-5,
            "max_dt_factor"                      : 1.0,
            "correct_backwards"                  : true,
            "maximum_iterations"                 : 5
        }""" )

        settings.ValidateAndAssignDefaults(default_parameters)
        convector_settings = KratosMultiphysics.Parameters( """{
            "linear_solver_settings": {
                "solver_type": "amgcl"
            },
            "max_cfl" : 1.0,
            "cross_wind_stabilization_factor" : 0.7,
            "max_substeps" : 5
        }""")
        settings["convector_settings"].ValidateAndAssignDefaults(convector_settings)
        self.settings = settings

        if ( self.settings["model_part_name"].GetString() == "" ):
            KratosMultiphysics.Logger.PrintInfo("ApplyMassConservationCheckProcess","The value (string) of the parameter 'model_part_name' must not be empty.")

        self._fluid_model_part = Model[self.settings["model_part_name"].GetString()]
        self._my_log_file = self.settings["log_file_name"].GetString()
        self._check_volume_loss_at_the_end = self.settings["check_volume_loss_at_the_end"].GetBool()
        self._correct_backwards = self.settings["correct_backwards"].GetBool()

        self._min_dt_factor = self.settings["min_dt_factor"].GetDouble()
        self._max_dt_factor = self.settings["max_dt_factor"].GetDouble()

        self._maximum_iterations = self.settings["maximum_iterations"].GetInt()

        self.mass_conservation_utility = KratosFluid.MassConservationUtility(self._fluid_model_part, self.settings["mass_correction_settings"])

        KratosMultiphysics.Logger.PrintInfo("ApplyMassConservationCheckProcess","Construction finished.")


    def ExecuteInitialize(self):
        self.forward_convection_process = self._set_levelset_convection_process()
        self.mass_conservation_utility.Initialize()
        KratosMultiphysics.Logger.PrintInfo("ApplyMassConservationCheckProcess","Initialization finished (initial volumes calculated).")

    def Execute(self):
        self.ExecuteFinalizeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        KratosMultiphysics.VariableUtils().SaveScalarVar(
            KratosMultiphysics.DISTANCE,
            KratosFluid.AUX_DISTANCE,
            self._fluid_model_part.Nodes)
        self.mass_conservation_utility.ComputeBalancedVolume()
        initial_volume_error = self.mass_conservation_utility.GetVolumeError()

        if not self._correct_backwards and initial_volume_error > 0.0:
            KratosMultiphysics.Logger.PrintInfo("Mass gained, but volume is not being corrected backwards in local corrector.")
            return

        flow_into_air = self.mass_conservation_utility.OrthogonalFlowIntoAir()
        prev_dt = self._fluid_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        dt_max = prev_dt*self._max_dt_factor if initial_volume_error < 0.0 else -prev_dt*self._max_dt_factor
        if flow_into_air < 0.0:
            dt_max = -dt_max

        self._ConvectAuxiliaryDistance(dt_max)
        self.mass_conservation_utility.ReCheckTheMassConservation()
        new_volume_error = self.mass_conservation_utility.GetVolumeError()

        if dt_max > 0.0:
            dt_low = 0.0
            dt_high = dt_max
            volume_error_low = initial_volume_error
            volume_error_high = new_volume_error
        else:
            dt_low = dt_max
            dt_high = 0.0
            volume_error_low = new_volume_error
            volume_error_high = initial_volume_error
        self.iteration_is_on_min_side = abs(volume_error_high) > abs(volume_error_low)
        self.old_iteration_was_on_min_side = not self.iteration_is_on_min_side
        self.gamma = 1.0


        if volume_error_low*volume_error_high > 0.0:
            if abs(volume_error_low) < abs(volume_error_high):
                self.mass_conservation_utility.RestoreDistanceValues(KratosFluid.AUX_DISTANCE)
            KratosMultiphysics.Logger.PrintInfo("Volume cannot be corrected, moving to global corrector")
            return

        iterations = 0

        while iterations < self._maximum_iterations:
            self.mass_conservation_utility.RestoreDistanceValues(KratosFluid.AUX_DISTANCE)
            dt = self._PredictNewDt(dt_low, dt_high, volume_error_low, volume_error_high)
            self._ConvectAuxiliaryDistance(dt)
            self.mass_conservation_utility.ReCheckTheMassConservation()
            new_volume_error = self.mass_conservation_utility.GetVolumeError()
            self.old_iteration_was_on_min_side = self.iteration_is_on_min_side
            if new_volume_error*volume_error_low > 0.0:
                volume_error_low = new_volume_error
                dt_low = dt
                self.iteration_is_on_min_side = True
            else:
                volume_error_high = new_volume_error
                dt_high = dt
                self.iteration_is_on_min_side = False
            if self.old_iteration_was_on_min_side is self.iteration_is_on_min_side:
                self.gamma *= 0.5
            else:
                self.gamma = 1.0
            iterations += 1

    def _ComputeTimeStepForConvection(self):
        """ REMARK: The pseudo time step dt has no meaning as a physical time step and is NOT accounted as such.
            Instead, it as purely used for an artificial convection that helps to conserve the mass.
            Accordingly, dt depends on the currently necessary mass convection.
        """
        dt = self.mass_conservation_utility.ComputeTimeStepForConvection()

        KratosMultiphysics.Logger.PrintInfo("Time step dt for convection in mass conservation: ", "{:e}".format(dt))
        return dt


    def _ConvectAuxiliaryDistance(self, dt):
        """ The distance field is artificially convected to extrapolate the interface motion into the future.
            The time step considered for this artificial "forward convection" is dt.
            The previous distance field is saved in an auxiliary variable
        """
        prev_dt = self._fluid_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        revert_velocity = False
        if dt < 0.0:
            dt *= -1.0
            revert_velocity = True
        self._fluid_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = dt
        if revert_velocity:
            self.mass_conservation_utility.RevertVelocityDirection()
        self.forward_convection_process.Execute()
        if revert_velocity:
            self.mass_conservation_utility.RevertVelocityDirection()
        self._fluid_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = prev_dt


    def _PredictNewDt(self, dt_low, dt_high, volume_error_low, volume_error_high):
        return (volume_error_high*dt_low - volume_error_low*dt_high)/(volume_error_high - volume_error_low)

    def _set_levelset_convection_process(self):
        if self._fluid_model_part.IsDistributed():
            self.EpetraCommunicator = KratosTrilinos.CreateCommunicator()
            self.trilinos_linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["convector_settings"]["linear_solver_settings"])
            level_set_convection_process = KratosTrilinos.TrilinosLevelSetConvectionProcess3D(
                self.EpetraCommunicator,
                KratosMultiphysics.DISTANCE,
                self._fluid_model_part,
                self.trilinos_linear_solver,
                self.settings["convector_settings"]["max_cfl"].GetDouble(),
                self.settings["convector_settings"]["cross_wind_stabilization_factor"].GetDouble(),
                self.settings["convector_settings"]["max_substeps"].GetInt())
        else:
            self.linear_solver = python_linear_solver_factory.ConstructSolver(self.settings["convector_settings"]["linear_solver_settings"])
            level_set_convection_process = KratosMultiphysics.LevelSetConvectionProcess3D(
                KratosMultiphysics.DISTANCE,
                self._fluid_model_part,
                self.linear_solver,
                self.settings["convector_settings"]["max_cfl"].GetDouble(),
                self.settings["convector_settings"]["cross_wind_stabilization_factor"].GetDouble(),
                self.settings["convector_settings"]["max_substeps"].GetInt())

        return level_set_convection_process