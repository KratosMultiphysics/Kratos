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

        self._my_log_file = "Mass_Conservation.log"

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "model_part_name"                        : "default_model_part_name",
            "mass_computation_frequency"             : 20,
            "compare_to_initial_values"              : true,
            "write_to_log_file"                      : true
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters)

        self._fluid_model_part = Model[settings["model_part_name"].GetString()]
        self._write_to_log = settings["write_to_log_file"].GetBool()
        self._is_printing_rank = ( self._fluid_model_part.GetCommunicator().MyPID() == 0 )

        self.mass_conservation_check_process = KratosFluid.MassConservationCheckProcess(self._fluid_model_part, settings)

        # writing first line in file
        if ( self._write_to_log and self._is_printing_rank ):
            with open(self._my_log_file, "w") as logFile:
                logFile.write( "positiveVolume" + "\t" + "negativeVolume" + "\n" )

        if self._is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("ApplyMassConservationCheckProcess","Construction of Mass Conservation Check Process finished.")


    def Execute(self):

        updated = int( self.mass_conservation_check_process.GetUpdateStatus() )

        # force a pause until the last process has its status
        self._fluid_model_part.GetCommunicator().Barrier()

        # would take the value "0" if one is not updated
        updated = self._fluid_model_part.GetCommunicator().MinAll( updated )

        if ( updated > 0 ):

            posVol = self.mass_conservation_check_process.GetPositiveVolume()
            negVol = self.mass_conservation_check_process.GetNegativeVolume()
            initPosVol = self.mass_conservation_check_process.GetInitialPositiveVolume()
            initNegVol = self.mass_conservation_check_process.GetInitialNegativeVolume()

            self._fluid_model_part.GetCommunicator().Barrier()
            posVol = self._fluid_model_part.GetCommunicator().SumAll( posVol )
            negVol = self._fluid_model_part.GetCommunicator().SumAll( negVol )
            initPosVol = self._fluid_model_part.GetCommunicator().SumAll( initPosVol )
            initNegVol = self._fluid_model_part.GetCommunicator().SumAll( initNegVol )

            # find the printing rank
            if ( self._is_printing_rank ):
                # managing the output to the console
                KratosMultiphysics.Logger.PrintInfo("ApplyMassConservationCheckProcess", "Positive Volume = " + str(posVol) + "  ( initially " + str(initPosVol) + ")" )
                KratosMultiphysics.Logger.PrintInfo("ApplyMassConservationCheckProcess", "Negative Volume = " + str(negVol) + "  ( initially " + str(initNegVol) + ")" )
                KratosMultiphysics.Logger.Flush()

                # adds additional lines to the log file
                if ( self._write_to_log ):
                    with open(self._my_log_file, "a+") as logFile:
                        logFile.write( str(posVol) + "\t" + str(negVol) + "\n" )


    def ExecuteInitialize(self):
        self.mass_conservation_check_process.ExecuteInitialize()


    def ExecuteBeforeSolutionLoop(self):
        self.mass_conservation_check_process.ExecuteBeforeSolutionLoop()


    def ExecuteInitializeSolutionStep(self):
        self.mass_conservation_check_process.ExecuteInitializeSolutionStep()


    def ExecuteFinalizeSolutionStep(self):

        self.Execute()