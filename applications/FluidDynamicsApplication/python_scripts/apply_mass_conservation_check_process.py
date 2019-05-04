from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if( not isinstance(settings, KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMassConservationCheckProcess( Model, settings["Parameters"] )

class ApplyMassConservationCheckProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "model_part_name"                           : "",
            "perform_local_corrections"                 : true,
            "perform_global_corrections"                : false,
            "correction_frequency_in_time_steps"        : 20,
            "write_to_log_file"                         : false,
            "log_file_name"                             : "mass_conservation.log"
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters)

        if ( settings["model_part_name"].GetString() == "" ):
            KratosMultiphysics.Logger.PrintInfo("ApplyMassConservationCheckProcess","The value (string) of the parameter 'model_part_name' must not be empty.")

        self._fluid_model_part = Model[settings["model_part_name"].GetString()]
        self._write_to_log = settings["write_to_log_file"].GetBool()
        self._my_log_file = settings["log_file_name"].GetString()

        self._perform_local_corr  = settings["perform_local_corrections"].GetBool()
        self._perform_global_corr = settings["perform_global_corrections"].GetBool()

        self._is_printing_rank = ( self._fluid_model_part.GetCommunicator().MyPID() == 0 )
        self.mass_conservation_check_process = KratosFluid.MassConservationCheckProcess(self._fluid_model_part, settings)

        KratosMultiphysics.Logger.PrintInfo("ApplyMassConservationCheckProcess","Construction finished.")


    def ExecuteInitialize(self):

        self.forward_convection_process = self._set_levelset_convection_process()

        # calling C++ initalization
        first_lines_string = self.mass_conservation_check_process.Initialize()

        # writing first line in file
        if ( self._write_to_log and self._is_printing_rank ):
            with open(self._my_log_file, "w") as logFile:
                logFile.write( first_lines_string )

        KratosMultiphysics.Logger.PrintInfo("ApplyMassConservationCheckProcess","Initialization finished (initial volumes calculated).")


    def ExecuteBeforeSolutionLoop(self):
        self.mass_conservation_check_process.ExecuteBeforeSolutionLoop()


    def ExecuteInitializeSolutionStep(self):
        self.mass_conservation_check_process.ExecuteInitializeSolutionStep()


    def ExecuteFinalizeSolutionStep(self):

        ### always necessary to keep track of the balanced volume
        log_line_string = self.mass_conservation_check_process.ComputeBalancedVolume()

        ### writing first line in file
        if ( self._write_to_log and self._is_printing_rank ):
            with open(self._my_log_file, "a+") as logFile:
                logFile.write( log_line_string )

        # (1) #
        if ( self._perform_local_corr ):
            ### perform the local conservation procedure
            # REMARK: The pseudo time step dt has no meaning as a physical time step and is NOT accounted as such.
            # Instead, it as purely used for an artificial convection that helps to conserve the mass.
            # Accordingly, dt depends on the currently necessary mass convection.
            dt = self.mass_conservation_check_process.ComputeDtForConvection()

            # REMARK: The distance field is artificially convected to extrapolate the inteface motion into the future.
            # The time step considered for this artificial "forward convection" is dt.
            # The generated distance field is saved in an auxiliary variable and does not influence the actual distance field.
            self.forward_convection_process.ConvectForward( dt, KratosFluid.AUX_DISTANCE )

            # REMARK: A comparison between the "forward_convected" distance field and the current distance field is made.
            # Depending on the current mass balance, one of the following options i chosen:
            # - taking min() of both fields = adding to fluid_1 ("more water")
            # - taking max() of both fields = adding to fluid_2 ("more air")
            self.mass_conservation_check_process.ApplyLocalCorrection( KratosFluid.AUX_DISTANCE )

        # (2) #
        if ( self._perform_local_corr and self._perform_global_corr ):
            ### check how much work has already been done by the local correction process
            self.mass_conservation_check_process.ReCheckTheMassConservation()

        # (3) #
        if ( self._perform_global_corr ):
            ### perform the global correction
            # without any consideration of the velocity field, volume is added by a global slight shift of the distance field
            self.mass_conservation_check_process.ApplyGlobalCorrection()


    def _set_levelset_convection_process(self):
        ### for serial and OpenMP
        serial_settings = KratosMultiphysics.Parameters("""{
            "linear_solver_settings"   : {
                "solver_type" : "amgcl"
            }
        }""")
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(serial_settings["linear_solver_settings"])

        if self._fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            level_set_convection_process = KratosMultiphysics.LevelSetForwardConvectionProcess2D(
                KratosMultiphysics.DISTANCE,
                self._fluid_model_part,
                self.linear_solver)
        else:
            level_set_convection_process = KratosMultiphysics.LevelSetForwardConvectionProcess3D(
                KratosMultiphysics.DISTANCE,
                self._fluid_model_part,
                self.linear_solver)

        return level_set_convection_process